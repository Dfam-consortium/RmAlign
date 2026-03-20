# RmAlign

Family-wise transposable element alignment framework — the computational
core of a RepeatMasker-like pipeline.  A lightweight Python wrapper
(`RepeatMasker`) handles argument parsing and environment checks, then
delegates all orchestration to Nextflow.

---

## Contents

- [How it works](#how-it-works)
- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Quickstart](#quickstart)
- [RepeatMasker CLI](#repeatmasker-cli)
- [generate_align_tasks.py](#generate_align_taskspy)
- [Execution profiles](#execution-profiles)
- [Output structure](#output-structure)
- [Utilities](#utilities)

---

## How it works

The pipeline runs in five stages.  Stages 3–5 are defined but not yet
wired into the active workflow.

| # | Nextflow process | Script | Status |
|---|---|---|---|
| 1 | `GEN_ALIGN_TASKS` | `generate_align_tasks.py` | active |
| 2 | `ALIGN_TASK` (fan-out) | `align_task.py` | active |
| 3 | `GEN_ADJ_TASKS` | `generate_adjudication_tasks.py` | placeholder |
| 4 | `ADJUDICATE_TASK` (fan-out) | `adjudicate_task.py` | placeholder |
| 5 | `POSTPROCESS_TASK` | `postprocess_task.py` | placeholder |

**Stage 1 — generate_align_tasks.py**

Validates the input sequence and repeat library, chunks the sequence
(optionally), groups chunks into batch FASTA files, and writes
`jobs/align_tasks.jsonl` — one JSON record per alignment job.

Chunking and batching are independent options:

| `--chunk-size` | `--batch-method` | Result |
|---|---|---|
| `60000` (default) | `gc_bin` (default) | 60 kb chunks grouped into `{bin-35gc,…,bin-53gc}.fa` by GC% — for **rmblastn** |
| `60000` | `per_seq` | one `seq#.fa` per chunk |
| `0` (no chunking) | `all` | single `all.fa` containing all contigs — for **nhmmer** |

Each run maintains a `manifest.json` that fingerprints the sequence and
library.  Re-runs skip batching/job-emission automatically when nothing
has changed.  `--aligner` and `--batch-method` are *breaking options*:
changing them on an existing project requires `--force-update`.

**Stage 2 — align_task.py**

Runs once per record in `align_tasks.jsonl` (all jobs run in parallel).
For each job it:

1. Extracts the target family from the library into a temporary FASTA.
2. Builds a single-family BLAST database (`makeblastdb`).
3. Selects the scoring matrix from `Matrices/ncbi/nt/` (e.g. `25p41g.matrix`);
   snaps to the nearest available divergence level if an exact match is absent.
4. Runs `rmblastn` via `telib.runners.rmblast`.
5. Annotates hits with corrected E-values.
6. Writes alignments to `results/ALIGN/<family>_bin-<gc>gc.bpaf` (binary BPAF format).

---

## Repository layout

```
RmAlign/
├── RepeatMasker              # user-facing entry-point (Python wrapper)
├── main.nf                   # Nextflow DSL2 workflow
├── nextflow.config           # manifest + includeConfig directives
├── _version.py               # single source of truth for version string
├── pyproject.toml            # Python package metadata (requires Python >= 3.10)
│
├── conf/
│   ├── base.config           # per-process CPU / memory / time limits
│   └── profiles.config       # executor profiles (standard, rmcluster)
│
├── bin/                      # scripts auto-added to $PATH by Nextflow
│   ├── generate_align_tasks.py
│   ├── align_task.py
│   ├── generate_adjudication_tasks.py
│   ├── adjudicate_task.py
│   └── postprocess_task.py
│
├── utils/
│   └── bpaf_converter.py     # standalone BPAF -> crossmatch/caf/rmblast/... converter
│
├── telib/                    # bundled alignment library
│   ├── formats/              # codecs: bpaf, crossmatch, caf, rmblast, dfam_tsv
│   ├── models/               # PairwiseAlignment, GapStats, ...
│   ├── runners/              # rmblast.py -- RmblastOptions, run_rmblast()
│   ├── scoring/              # evalue.py, kimura.py, nishimaki.py
│   ├── sequences/            # FASTA / 2bit readers and writer
│   └── tests/
│
├── Matrices/ncbi/nt/         # RMBlast scoring matrices ({div}p{gc}g.matrix)
└── data/                     # small example files for testing
```

---

## Requirements

| Component | Minimum version |
|---|---|
| Python | 3.10 |
| Nextflow | 25.04.6 |
| Java | 24.0.2 |
| rmblastn | any recent RMBlast build |
| makeblastdb | bundled with RMBlast / NCBI BLAST+ |

Both `Nextflow` and `Java` versions are checked at startup by the
`RepeatMasker` wrapper.

---

## Quickstart

```bash
# Local execution (all tasks run on the current machine)
./RepeatMasker \
    --sequence    data/example-1mb.fa \
    --project_dir my_project \
    --library     data/lib.fa

# SLURM cluster (rmcluster partition)
./RepeatMasker \
    --sequence    data/example-1mb.fa \
    --project_dir my_project \
    --library     data/lib.fa \
    --nf_profile  rmcluster \
    --cpus        8

# Re-run after updating the library (only changed families re-aligned)
./RepeatMasker \
    --sequence    data/example-1mb.fa \
    --project_dir my_project \
    --library     data/lib.new.fa
```

On a re-run, `generate_align_tasks.py` compares BLAKE2b signatures of
every library sequence against the stored manifest.  Only families that
were added or modified produce new alignment jobs; unchanged families
are skipped automatically.

---

## RepeatMasker CLI

```
./RepeatMasker --sequence <FASTA> --project_dir <DIR>
               --library <FA>
               [--cpus N] [--nf_profile PROFILE]
               [--search-threshold INT] [--final-threshold INT]
               [--force-update] [--debug] [--version]
```

| Flag | Default | Description |
|---|---|---|
| `--sequence` | — | Input FASTA (required) |
| `--project_dir` | — | Project work directory; created if absent (required) |
| `--library` | — | Repeat library FASTA (required) |
| `--cpus N` | `4` | Threads per alignment task |
| `--nf_profile PROFILE` | `standard` | Nextflow executor profile (`standard` or `rmcluster`) |
| `--search-threshold INT` | `180` | Minimum rmblastn raw gapped score |
| `--no-family-threshold` | off | Disable per-family score thresholds from library headers |
| `--force-update` | off | Accept sequence / library / breaking-option changes |
| `--debug` | off | Enable Nextflow `-trace` output |

---

## generate_align_tasks.py

Can also be run directly, outside of Nextflow, for testing or custom
pipelines.

```bash
bin/generate_align_tasks.py \
    --sequence    data/example-1mb.fa \
    --library     data/lib.fa \
    --project_dir my_project \
    --chunk-size  60000 \
    --batch-method gc_bin \
    --aligner     rmblastn
```

### Key options

| Flag | Default | Description |
|---|---|---|
| `--sequence` | — | Input FASTA; `.gz`/`.bgz` accepted (required) |
| `--library` | — | Repeat library FASTA (required on first run) |
| `--project_dir` | — | Work directory (required) |
| `--chunk-size N` | `60000` | Bases per chunk; `0` = no chunking (one chunk per contig) |
| `--overlap N` | `0` | Overlap between consecutive chunks |
| `--batch-method` | `gc_bin` | `gc_bin` / `per_seq` / `all` — **breaking** |
| `--gc-denom` | `acgt` | GC% denominator for bin assignment: `acgt` (exclude Ns) or `all` (include Ns, matches legacy Perl behavior) — **breaking** |
| `--aligner` | `rmblastn` | `rmblastn` / `nhmmer` — stored in manifest, **breaking** |
| `--search-threshold` | `180` | Minimum rmblastn `-min_raw_gapped_score` |
| `--no-family-threshold` | off | Ignore `s_thresh=` tags in library FASTA headers |
| `--force-update` | off | Accept changes to sequence, library, or breaking options |
| `--skip-all-n` | off | Drop chunks composed entirely of N |
| `--jobs-format` | `jsonl` | `jsonl` / `tsv` |

### Library header tags

`generate_align_tasks.py` reads two optional key=value tokens from each
library FASTA description line:

```
>AluYa5 div=18 s_thresh=225
```

- `div=N` — divergence level used to select the scoring matrix.
- `s_thresh=N` — post-alignment score threshold for that family.  Families
  without this tag default to 225 (unless `--no-family-threshold` is set).

### Manifest and incremental updates

`<project_dir>/manifest.json` stores:

- Sequence file identity (name, size, mtime, partial and full hash).
- Per-family BLAKE2b signatures of the library.
- Breaking options (`aligner`, `batch_method`).

On re-runs the manifest is compared against the current inputs:

- **No changes** — short-circuits immediately (empty jobs file, no batching).
- **Library changes** — only added/modified families get new jobs.
- **Sequence changed / version upgrade / breaking option changed** — full
  re-alignment; requires `--force-update`.

---

## Execution profiles

Defined in `conf/profiles.config`:

| Profile | Executor | Notes |
|---|---|---|
| `standard` | local | Default; all processes run on the current machine |
| `rmcluster` | SLURM | Submits to the `rmcluster` partition; queue size 500, rate 20/min |

Resource limits per process are in `conf/base.config`:

| Process | CPUs | Memory | Time |
|---|---|---|---|
| `GEN_ALIGN_TASKS` | 1 | 4 GB | 1 h |
| `ALIGN_TASK` | `params.threads` | 16 GB | 4 h |
| `GEN_ADJ_TASKS` | 1 | 4 GB | 1 h |
| `ADJUDICATE_TASK` | 1 | 8 GB | 4 h |
| `POSTPROCESS_TASK` | 1 | 8 GB | 2 h |

---

## Output structure

After a successful run the project directory contains:

```
my_project/
├── manifest.json               # sequence + library fingerprints, breaking options
├── sequence/
│   ├── lib-cons.fa             # project-local copy of the repeat library
│   ├── bins.tsv                # chunk table (id, coords, GC%, matrix bin, ...)
│   ├── bin-35gc.fa ... bin-53gc.fa  # GC-binned batch FASTAs  (gc_bin mode)
│   ├── seq1.fa  seq2.fa ...    # per-chunk FASTAs         (per_seq mode)
│   └── all.fa                  # single combined FASTA    (all mode)
├── jobs/
│   └── align_tasks.jsonl       # one JSON record per alignment job
└── results/
    └── ALIGN/
        └── <family>_bin-<gc>gc.bpaf   # binary alignment output per job
```

### align_tasks.jsonl record fields

| Field | Modes | Description |
|---|---|---|
| `family` | all | Library family ID |
| `sequence` | all | Batch FASTA filename (relative to `sequence/`) |
| `gc` | `gc_bin` only | Integer GC% of the bin |
| `div` | all | Divergence level for matrix selection |
| `search_threshold` | all | Minimum raw score for rmblastn |
| `final_threshold` | all | Post-alignment score filter (0 = keep all) |
| `bin_bases` | all | Non-N bases in the batch FASTA (for E-value calculation) |
| `full_seq_bases` | all | Total non-N bases across the entire sequence |

---

## Utilities

### utils/bpaf_converter.py

Converts BPAF files produced by `align_task.py` to human-readable formats.

```bash
# Project-aware (resolves files from manifest automatically)
utils/bpaf_converter.py --project my_project --family AluYa5
utils/bpaf_converter.py --project my_project --family AluYa5 --gc 43 --out caf

# Direct file mode
utils/bpaf_converter.py \
    --input   my_project/results/ALIGN/AluYa5_43g.bpaf \
    --query   my_project/sequence/43g.fa \
    --library data/lib.fa \
    --out     crossmatch \
    --include-alignment
```

Supported output formats: `crossmatch` (default), `caf`, `bpaf`, `rmblast`, `cigar`, `aligned`.
