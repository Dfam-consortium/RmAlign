# RmAlign

This is a project to work on the architecture of a RepeatMasker-like 
family-wise alignment framework.  It has a main entry point script
('RepeatMasker') which sets up a project directory (if it doesn't
already exist) and turns orchestration over to Nextflow.

## Modular Workflow Pattern

Single workflow-spec (cluster and server-based) in Nextflow, with a 
lightweight wrapper (BASH) script to "rebrand" the command *and*
run Nextflow inside of the project-directory.

PROS:
- Single workflow description in DSL2
- Nextflow logs can be directed to the project directory without user-intervention to 
facilitate easier user support/debugging.
- Single ("rebranded") entry script simply checks that it can run Nextflow, change into
the project-directory and passes all the options to a "nextflow run" command.
- Purely a Nextflow based application that should run in any nf-core compatible environment

CONS:
- Restricted to Nextflow logging, and error reporting mechanisms -- which can be cryptic.
- Users must install/configure Nextflow for any type of run.
- Builds a larger barier to use of future parallel execution engines (e.g. Snakemake, Cromwell)

---

## Contents

- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Quickstart (demo)](#quickstart-demo)

---

## What it does

Stages (same for both Python + Nextflow):

1. **Generate align tasks** — `generate_align_tasks.py` writes `align_tasks.json` at the top of `--project_dir`.
2. **Align fan-out** — `align_task.py` runs once per record in `align_tasks.json` (in parallel).
3. **Generate adjudication tasks** — `generate_adjudication_tasks.py` writes `adjudication_tasks.json` at the top of `--project_dir`.
4. **Adjudicate fan-out** — `adjudicate_task.py` runs once per record in `adjudication_tasks.json` (in parallel).
5. **Post-process** — `postprocess_task.py` writes a small summary under `results/SUMMARIZE/`.

All scripts accept the **same base CLI**:

```
--sequence <file>   --project_dir <dir>   --species <string>
```

`--project_dir` is created if it doesn't exist.

---

## Repository layout

```
.
├─ main.nf                               # Nextflow DSL2 pipeline
├─ nextflow_schema.json                  # Param schema (validates --sequence/--project_dir/--species)
├─ RepeatMasker                          # Python - main entry point for workflow
├─ generate_align_tasks.py               # writes project_dir/align_tasks.json
├─ align_task.py                         # consumes one align task
├─ generate_adjudication_tasks.py        # writes project_dir/adjudication_tasks.json
├─ adjudicate_task.py                    # consumes one adjudication task
├─ postprocess_task.py                   # final summary step
└─ data/example.fa                       # tiny demo input
```

> Tip: You may also place helper scripts in a `bin/` folder and make them executable; Nextflow adds `bin/` to `PATH` automatically for all processes.

---

## Requirements

- **Python** 3.8+
- **Nextflow** 25.04.x (tested). Other recent versions likely fine; pin if needed.
- Bash / POSIX shell

Optional (recommended):

- A Python virtual environment for running helper scripts locally.

---

## Quickstart (demo)

Create a project directory and run both paths end-to-end.

```bash
# (1) Python orchestrator demo
./RepeatMasker \
  --sequence data/example.fa \
  --project_dir my_project \
  --species human
```
