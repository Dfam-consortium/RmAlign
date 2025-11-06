// main.nf
nextflow.enable.dsl = 2

import groovy.json.JsonSlurper
import groovy.json.JsonOutput

// --------------------
// Params
// --------------------
params.sequence     = params.sequence ?: null
params.project_dir  = params.project_dir ?: null
params.lib_file     = params.lib_file   ?: null

// --------------------
// Helpers
// --------------------
def fail(msg) {
    log.error msg
    System.exit(1)
}

def ensureProjectDir(p) {
    def d = new File(p as String)
    if (!d.exists()) {
        if (!d.mkdirs())
            fail("Could not create --project_dir: ${d}")
        log.info "Created project directory: ${d.absolutePath}"
    } else if (!d.isDirectory()) {
        fail("--project_dir exists but is not a directory: ${d}")
    }
    return d
}

def mapToCliFlags(Map rec) {
    rec.collect { k, v ->
        def flag = "--${k.toString().replace('_','-')}"
        switch (v) {
            case null:            return null
            case Boolean:         return v ? flag : null
            case List:            return v.collect { "${flag} ${it}" }.join(' ')
            case Map:             return "${flag} '${JsonOutput.toJson(v)}'"
            default:              return "${flag} ${v}"
        }
    }.findAll { it }.join(' ')
}

// --------------------
// Early validation
// --------------------
if (!params.sequence)     fail("Missing required parameter: --sequence <file>")
if (!params.project_dir)  fail("Missing required parameter: --project_dir <directory>")
if (!params.lib_file)     fail("Missing required parameter: --lib_file <string>")

def seqPath = file(params.sequence)
if (!seqPath.exists())    fail("The --sequence path does not exist: ${seqPath}")
if (!seqPath.isFile())    fail("The --sequence path is not a file: ${seqPath}")

def projDir = ensureProjectDir(params.project_dir)

log.info "Using sequence    : ${seqPath.toAbsolutePath()}"
log.info "Using project_dir : ${projDir.absolutePath}"
log.info "Using lib_file    : ${params.lib_file}"

// Singletons
Channel.value(seqPath.toAbsolutePath().toString()).set { CH_SEQ }
Channel.value(projDir.absolutePath).set             { CH_PROJ }
Channel.value(params.lib_file as String).set        { CH_LIB }

// Absolute paths to helper scripts in the repo
def GEN_ALIGN_SCRIPT = "${projectDir}/generate_align_tasks.py"
def ALIGN_SCRIPT     = "${projectDir}/align_task.py"
def GEN_ADJ_SCRIPT   = "${projectDir}/generate_adjudication_tasks.py"
def ADJ_SCRIPT       = "${projectDir}/adjudicate_task.py"
def POST_SCRIPT      = "${projectDir}/postprocess_task.py"

// --------------------
// Workflow
// --------------------
workflow {

    // Stage 1: generate_align_tasks.py -> produce align_tasks.json and project token
    GEN_ALIGN_TASKS(CH_SEQ, CH_PROJ, CH_LIB)

    // Parse align_tasks.json and fan out ALIGN tasks
    //GEN_ALIGN_TASKS.out.align_json
    //    .map { f ->
    //        def list = new JsonSlurper().parse(new File(f.toString()))
    //        if (!(list instanceof List))
    //            throw new IllegalStateException("align_tasks.json must be a JSON array")
    //        list
    //    }
    //    .flatten()
    //    .combine( GEN_ALIGN_TASKS.out.proj_after_align )
    //    .combine( CH_SEQ )
    //    .combine( CH_LIB )
    //    .map { rec, proj, seq, lib_file -> [ rec as Map, seq as String, proj as String, lib_file as String ] }
    //    .set { CH_ALIGN_TUPLES }

    // Parse align_tasks.jsonl (JSON Lines) and fan out ALIGN tasks
    GEN_ALIGN_TASKS.out.align_json
        .map { f ->
            def slurper = new JsonSlurper()
            def file = (f instanceof java.nio.file.Path) ? f.toFile() : new File(f.toString())
            def items = []

            file.withReader { rdr ->
                rdr.eachLine { line ->
                    def s = line?.trim()
                    if (!s) return  // skip blank lines
                    def obj = slurper.parseText(s)
                    if (!(obj instanceof Map))
                        throw new IllegalStateException("Each line in align_tasks.jsonl must be a JSON object")
                    items << (obj as Map)
                }
            }
            items
        }
        .flatten()
        .combine( GEN_ALIGN_TASKS.out.proj_after_align )
        .combine( CH_SEQ )
        .combine( CH_LIB )
        .map { rec, proj, seq, lib_file -> [ rec as Map, seq as String, proj as String, lib_file as String ] }
        .set { CH_ALIGN_TUPLES }


    ALIGN_TASK(CH_ALIGN_TUPLES)

    // Barrier: wait for all align tasks to finish
    ALIGN_TASK.out.done_tokens.collect().set { CH_ALIGN_DONE_LIST }

    // Stage 3: generate_adjudication_tasks.py (after aligns)
    GEN_ADJ_TASKS(CH_SEQ, GEN_ALIGN_TASKS.out.proj_after_align, CH_LIB, CH_ALIGN_DONE_LIST)

    // Parse adjudication_tasks.json and fan out ADJUDICATE tasks
    GEN_ADJ_TASKS.out.adj_json
        .map { f ->
            def list = new JsonSlurper().parse(new File(f.toString()))
            if (!(list instanceof List))
                throw new IllegalStateException("adjudication_tasks.json must be a JSON array")
            list
        }
        .flatten()
        .combine( GEN_ADJ_TASKS.out.proj_after_adj )
        .combine( CH_SEQ )
        .combine( CH_LIB )
        .map { rec, proj, seq, species -> [ rec as Map, seq as String, proj as String, species as String ] }
        .set { CH_ADJ_TUPLES }

    ADJUDICATE_TASK(CH_ADJ_TUPLES)

    // Barrier: wait for all adjudications to finish
    ADJUDICATE_TASK.out.done_tokens.collect().set { CH_ADJ_DONE_LIST }

    // Stage 5: postprocess_task.py (once, after adjudications)
    POSTPROCESS_TASK(CH_SEQ, GEN_ADJ_TASKS.out.proj_after_adj, CH_LIB, CH_ADJ_DONE_LIST)
}

// --------------------
// Processes
// --------------------

// Stage 1: generate and expose align_tasks.json
process GEN_ALIGN_TASKS {
    tag { "generate_align_tasks" }

    input:
    val seq
    val proj
    val lib_file

    output:
    path "align_tasks.json", emit: align_json
    val(proj),             emit: proj_after_align

    script:
    """
    set -euo pipefail
    python3 "${GEN_ALIGN_SCRIPT}" \\
      --sequence ${seq} \\
      --project_dir ${proj} \\
      --emit-batches-gcbin \\
      --library ${lib_file}
    # Expose the JSON file to Nextflow as a process output
    ln -sf "${proj}/jobs/align_tasks.jsonl" align_tasks.json
    """
}

// Stage 2: align fan-out
process ALIGN_TASK {

    tag { "align_${rec.family ?: 'x'}__${rec.sequence ?: 'x'}" }

    input:
    tuple val(rec), val(seq), val(proj), val(species)

    output:
    path "*.align.done", emit: done_tokens

    when:
    rec != null

    script:
    def flags = mapToCliFlags(rec)
    def fam   = (rec.family ?: 'x').toString()
    def label = (rec.sequence ?: 'x').toString()
    def token = "align_${fam}__${label}.align.done"
    """
    python3 "${ALIGN_SCRIPT}" \\
      --output ${proj}/results/ALIGN/${fam}__${label}.out \\
      --sequence ${proj}/sequence/${rec.sequence} \\
      --family ${fam} \\
      --div ${rec.div} \\
      --gc  ${rec.gc} \\
      --threshold ${rec.threshold} \\
      --library /u3/home/rhubley/projects/RepeatMasker5/data/lib.fa \\
      --threads 4
    echo done > ${token}
    """
}

// Stage 3: generate and expose adjudication_tasks.json
process GEN_ADJ_TASKS {
    tag { "generate_adjudication_tasks" }

    input:
    val seq
    val proj
    val species
    val align_done_list    // barrier list; not used directly

    output:
    path "adjudication_tasks.json", emit: adj_json
    val(proj),                    emit: proj_after_adj

    script:
    """
    set -euo pipefail
    # barrier consumed via 'align_done_list'
    python3 "${GEN_ADJ_SCRIPT}" \\
      --sequence ${seq} \\
      --project_dir ${proj} \\
      --species ${species}
    ln -sf "${proj}/adjudication_tasks.json" adjudication_tasks.json
    """
}

// Stage 4: adjudicate fan-out
process ADJUDICATE_TASK {

    tag { "adjudicate_${rec.batch ?: 'x'}" }

    input:
    tuple val(rec), val(seq), val(proj), val(species)

    output:
    path "*.adj.done", emit: done_tokens

    when:
    rec != null

    script:
    def flags = mapToCliFlags(rec)
    def batch = (rec.batch ?: 'x').toString()
    def token = "adjudicate_${batch}.adj.done"
    """
    set -euo pipefail
    python3 "${ADJ_SCRIPT}" \\
      --sequence ${seq} \\
      --project_dir ${proj} \\
      --species ${species} \\
      ${flags}
    echo done > ${token}
    """
}

// Stage 5: postprocess (single)
process POSTPROCESS_TASK {
    tag { "postprocess" }

    input:
    val seq
    val proj
    val species
    val adj_done_list      // barrier list; not used directly

    output:
    stdout

    script:
    """
    set -euo pipefail
    python3 "${POST_SCRIPT}" \\
      --sequence ${seq} \\
      --project_dir ${proj} \\
      --species ${species}
    """
}
