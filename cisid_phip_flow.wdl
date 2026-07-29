version 1.0

workflow phippery_flow {
    input {
        File    input_sample_table
        File    input_peptide_table
        String  output_directory
        Int     read_length               = 125
        Int     oligo_tile_length         = 117
        Int     n_mismatches              = 2
        String  output_prefix             = "data"
        Boolean replicate_sequence_counts = true
        Boolean run_beer                  = false
        Boolean run_cpm_enrichment        = true
        Boolean run_z_score               = true
        String  zone                      = "us-central1-d"
        String  memory                    = "32G"
        Int     num_cpu                   = 8
        Int     preemptible               = 2
        Int     disk_space                = 250
        String  docker_registry
    }

    # ── Step 1: parse sample table into arrays for scatter ──────────────────
    call parse_sample_table {
        input:
            sample_table    = input_sample_table,
            docker_registry = docker_registry
    }

    # ── Step 2: download + merge each sample in parallel ────────────────────
    scatter (i in range(length(parse_sample_table.sample_ids))) {
        call download_and_merge {
            input:
                sample_id       = parse_sample_table.sample_ids[i],
                r1_gcs          = parse_sample_table.r1_paths[i],
                r2_gcs          = parse_sample_table.r2_paths[i],
                docker_registry = docker_registry,
                preemptible     = preemptible,
                zone            = zone
        }
    }

    # ── Step 3: gather all merged FASTQs, run full phip-flow cohort ─────────
    call run_phippery_flow {
        input:
            merged_fastqs         = download_and_merge.merged_fastq,
            sample_ids_file       = write_lines(parse_sample_table.sample_ids),
            control_statuses_file = write_lines(parse_sample_table.control_statuses),
            input_peptide_table   = input_peptide_table,
            output_directory      = output_directory,
            read_length           = read_length,
            oligo_tile_length     = oligo_tile_length,
            n_mismatches          = n_mismatches,
            output_prefix         = output_prefix,
            replicate_sequence_counts = replicate_sequence_counts,
            run_beer              = run_beer,
            run_cpm_enrichment    = run_cpm_enrichment,
            run_z_score           = run_z_score,
            zone                  = zone,
            memory                = memory,
            num_cpu               = num_cpu,
            preemptible           = preemptible,
            disk_space            = disk_space,
            docker_registry       = docker_registry
    }

    output {
        File outs = run_phippery_flow.outs_files
    }
}


# ════════════════════════════════════════════════════════════════════════════
# Task 1 — Parse sample table CSV → parallel arrays for scatter
# Small VM; only reads one CSV file.
# ════════════════════════════════════════════════════════════════════════════
task parse_sample_table {
    input {
        File   sample_table
        String docker_registry
    }

    command <<<
        python3 <<CODE
import csv, sys

with open("~{sample_table}") as f:
    rows = list(csv.DictReader(f))

required = {"sample_ID", "R1_cloud_filepath", "R2_cloud_filepath", "control_status"}
missing  = required - set(rows[0].keys()) if rows else required
if missing:
    print(f"ERROR: sample table is missing columns: {missing}", file=sys.stderr)
    sys.exit(1)

for col, fname in [
    ("sample_ID",         "sample_ids.txt"),
    ("R1_cloud_filepath", "r1_paths.txt"),
    ("R2_cloud_filepath", "r2_paths.txt"),
    ("control_status",    "control_statuses.txt"),
]:
    with open(fname, "w") as out:
        out.write("\n".join(r[col] for r in rows))

print(f"Parsed {len(rows)} samples.")
CODE
    >>>

    output {
        Array[String] sample_ids       = read_lines("sample_ids.txt")
        Array[String] r1_paths         = read_lines("r1_paths.txt")
        Array[String] r2_paths         = read_lines("r2_paths.txt")
        Array[String] control_statuses = read_lines("control_statuses.txt")
    }

    runtime {
        docker:      docker_registry
        cpu:         1
        memory:      "2G"
        disks:       "local-disk 10 HDD"
        preemptible: 3
        zones:       "us-central1-d"
    }
}


# ════════════════════════════════════════════════════════════════════════════
# Task 2 — Per-sample download + BBMerge (runs in parallel across samples)
# Small VM per sample; disk sized for one sample's R1+R2+output.
# ════════════════════════════════════════════════════════════════════════════
task download_and_merge {
    input {
        String sample_id
        String r1_gcs
        String r2_gcs
        String docker_registry
        Int    preemptible
        String zone
    }

    command <<<
        set -euo pipefail

        mkdir -p /tmp/fastqs

        # ── Download with 3-attempt retry ──────────────────────────────────
        echo "[~{sample_id}] Downloading R1..."
        gcloud storage cp "~{r1_gcs}" /tmp/fastqs/R1.fastq.gz || \
        gcloud storage cp "~{r1_gcs}" /tmp/fastqs/R1.fastq.gz || \
        gcloud storage cp "~{r1_gcs}" /tmp/fastqs/R1.fastq.gz

        echo "[~{sample_id}] Downloading R2..."
        gcloud storage cp "~{r2_gcs}" /tmp/fastqs/R2.fastq.gz || \
        gcloud storage cp "~{r2_gcs}" /tmp/fastqs/R2.fastq.gz || \
        gcloud storage cp "~{r2_gcs}" /tmp/fastqs/R2.fastq.gz

        # ── BBMerge: keep merged + unmerged reads ──────────────────────────
        # Overlapping reads (short inserts) → merged consensus
        # Non-overlapping reads (long inserts) → both R1 and R2 kept as SE
        echo "[~{sample_id}] Running BBMerge..."

        set +e
        bbmerge.sh \
            in=/tmp/fastqs/R1.fastq.gz    \
            in2=/tmp/fastqs/R2.fastq.gz   \
            out=/tmp/fastqs/merged.fastq.gz   \
            outu1=/tmp/fastqs/R1u.fastq.gz    \
            outu2=/tmp/fastqs/R2u.fastq.gz    \
            2>"~{sample_id}_bbmerge.log"
        BB_EXIT=$?
        set -e

        if [ "$BB_EXIT" -ne 0 ]; then
            echo "  BBMerge failed (exit $BB_EXIT) — falling back to raw R1+R2 concatenation"
            cat /tmp/fastqs/R1.fastq.gz /tmp/fastqs/R2.fastq.gz > "~{sample_id}.fastq.gz"
        else
            cat /tmp/fastqs/merged.fastq.gz \
                /tmp/fastqs/R1u.fastq.gz    \
                /tmp/fastqs/R2u.fastq.gz    > "~{sample_id}.fastq.gz"
        fi

        # ── Sanity check ───────────────────────────────────────────────────
        N_READS=$(( $(zcat "~{sample_id}.fastq.gz" | wc -l) / 4 ))
        echo "[~{sample_id}] $N_READS reads written to ~{sample_id}.fastq.gz"

        if [ "$N_READS" -eq 0 ]; then
            echo "ERROR: ~{sample_id} produced 0 reads — aborting" >&2
            exit 1
        fi
    >>>

    output {
        File merged_fastq = "~{sample_id}.fastq.gz"
        File bbmerge_log  = "~{sample_id}_bbmerge.log"
    }

    runtime {
        docker:      docker_registry
        cpu:         2
        memory:      "8G"
        disks:       "local-disk 50 HDD"
        preemptible: preemptible
        zones:       zone
    }
}


# ════════════════════════════════════════════════════════════════════════════
# Task 3 — Gather all merged FASTQs, run phip-flow over the full cohort
# Larger VM; phip-flow alignment + edgeR/BEER runs here.
# ════════════════════════════════════════════════════════════════════════════
task run_phippery_flow {
    input {
        Array[File] merged_fastqs
        File        sample_ids_file
        File        control_statuses_file
        File        input_peptide_table
        String      output_directory
        Int         read_length               = 125
        Int         oligo_tile_length         = 117
        Int         n_mismatches              = 2
        String      output_prefix             = "data"
        Boolean     replicate_sequence_counts = true
        Boolean     run_beer                  = false
        Boolean     run_cpm_enrichment        = true
        Boolean     run_z_score               = true
        String      zone
        String      memory
        Int         num_cpu
        Int         preemptible
        Int         disk_space
        String      docker_registry
    }

    output {
        File outs_files = "${output_prefix}_outs.txt"
    }

    command <<<
        set -euo pipefail

        cd /phipflow/

        mkdir -p /phipflow/data/seq

        cp ~{input_peptide_table} /phipflow/data/peptide_table.csv

        # ── Localize all merged FASTQs into the seq directory ──────────────
        # Cromwell stages Array[File] inputs but puts them at arbitrary paths.
        # Copy each one to /phipflow/data/seq/<sample_id>.fastq.gz using the
        # order-matched sample_ids array (WDL scatter preserves order).
        python3 <<CODE
import shutil

fastq_paths  = "~{sep=',' merged_fastqs}".split(",")
sample_ids   = open("~{sample_ids_file}").read().splitlines()

for path, sid in zip(fastq_paths, sample_ids):
    dest = f"/phipflow/data/seq/{sid}.fastq.gz"
    shutil.copy(path.strip(), dest)
    print(f"  Staged {sid}.fastq.gz")
CODE

        # ── Rebuild sample table with fastq_filepath for phip-flow ─────────
        python3 <<CODE
sample_ids       = open("~{sample_ids_file}").read().splitlines()
control_statuses = open("~{control_statuses_file}").read().splitlines()

with open("/phipflow/data/sample_table.csv", "w") as f:
    f.write("control_status,sample_ID,fastq_filepath\n")
    for sid, cs in zip(sample_ids, control_statuses):
        f.write(f"{cs},{sid},{sid}.fastq.gz\n")

print(f"Wrote sample_table.csv with {len(sample_ids)} samples.")
CODE

        df -h

        # ── Run phip-flow ──────────────────────────────────────────────────
        # Nextflow work dir goes on the mounted data disk, not the boot disk
        CMD="nextflow run main.nf"
        CMD="$CMD -w /mnt/disks/cromwell_root/nf_work"
        CMD="$CMD --ansi-log false"
        CMD="$CMD --sample_table /phipflow/data/sample_table.csv"
        CMD="$CMD --peptide_table /phipflow/data/peptide_table.csv"
        CMD="$CMD --reads_prefix /phipflow/data/seq"

        if [[ ~{read_length} -ne 125 ]]; then
            CMD="$CMD --read_length ~{read_length}"
        fi
        if [[ ~{oligo_tile_length} -ne 117 ]]; then
            CMD="$CMD --oligo_tile_length ~{oligo_tile_length}"
        fi
        if [[ ~{n_mismatches} -ne 2 ]]; then
            CMD="$CMD --n_mismatches ~{n_mismatches}"
        fi

        [[ ~{replicate_sequence_counts} == true ]] && CMD="$CMD --replicate_sequence_counts"
        [[ ~{run_beer} == true ]]                  && CMD="$CMD --run_BEER"
        [[ ~{run_cpm_enrichment} == true ]]        && CMD="$CMD --run_cpm_enr_workflow"
        [[ ~{run_z_score} == true ]]               && CMD="$CMD --run_z_score_fit_predict"

        echo "Running: $CMD"
        eval $CMD

        # ── Decompress result files before export ─────────────────────────
        find /phipflow/results/ -name "*.gz" -exec gunzip {} \; 2>/dev/null || true

        # ── Upload results + logs ──────────────────────────────────────────
        gcloud storage cp -r /phipflow/results/ ~{output_directory}/
        gcloud storage cp /phipflow/data/seq/*_bbmerge.log ~{output_directory}/bbmerge_logs/ || true

        # Sentinel file created locally so Cromwell tracks it as a real output
        touch "~{output_prefix}_outs.txt"
        gcloud storage cp "~{output_prefix}_outs.txt" ~{output_directory}/
    >>>

    runtime {
        bootDiskSizeGb: 50
        disks:          "local-disk ~{disk_space} HDD"
        docker:         "~{docker_registry}"
        cpu:            num_cpu
        zones:          zone
        memory:         memory
        preemptible:    preemptible
    }
}
