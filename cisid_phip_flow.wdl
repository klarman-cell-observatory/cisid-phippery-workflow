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

    call run_phippery_flow {
        input:
            input_sample_table        = input_sample_table,
            input_peptide_table       = input_peptide_table,
            output_directory          = output_directory,
            read_length               = read_length,
            oligo_tile_length         = oligo_tile_length,
            n_mismatches              = n_mismatches,
            output_prefix             = output_prefix,
            replicate_sequence_counts = replicate_sequence_counts,
            run_beer                  = run_beer,
            run_cpm_enrichment        = run_cpm_enrichment,
            run_z_score               = run_z_score,
            zone                      = zone,
            memory                    = memory,
            num_cpu                   = num_cpu,
            preemptible               = preemptible,
            disk_space                = disk_space,
            docker_registry           = docker_registry
    }

    output {
        File outs = run_phippery_flow.outs_files
    }
}


task run_phippery_flow {
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
        String  zone
        String  memory
        Int     num_cpu
        Int     preemptible
        Int     disk_space
        String  docker_registry
    }

    output {
        File outs_files = "~{output_prefix}_outs.txt"
    }

    command <<<
        set -euo pipefail

        cd /phipflow/

        cp ~{input_sample_table}  /phipflow/data/sample_table.csv
        cp ~{input_peptide_table} /phipflow/data/peptide_table.csv

        mkdir -p /phipflow/data/seq
        mkdir -p /mnt/disks/cromwell_root/tmp_fastqs

        # ── Download + merge all samples in parallel (16 workers) ─────────────
        python3 <<CODE
import csv, os, sys, subprocess, threading
from concurrent.futures import ThreadPoolExecutor, as_completed

seq_dir = "/phipflow/data/seq"
tmp_dir = "/mnt/disks/cromwell_root/tmp_fastqs"

_lock = threading.Lock()
def log(msg):
    with _lock:
        print(msg, flush=True)

with open("/phipflow/data/sample_table.csv") as f:
    rows = list(csv.DictReader(f))

required = {"sample_ID", "R1_cloud_filepath", "R2_cloud_filepath", "control_status"}
missing  = required - set(rows[0].keys()) if rows else required
if missing:
    print(f"ERROR: sample table missing columns: {missing}", file=sys.stderr)
    sys.exit(1)

def download(gcs_path, local_path, label, sample_id):
    for attempt in range(1, 4):
        r = subprocess.run(
            ["gcloud", "storage", "cp", gcs_path, local_path],
            capture_output=True, text=True
        )
        if r.returncode == 0:
            return
        log(f"[{sample_id}] {label} attempt {attempt} failed: {r.stderr.strip()}")
    raise RuntimeError(f"Failed to download {label} for {sample_id} after 3 attempts")

def process_sample(row):
    sample_id = str(row["sample_ID"])
    r1_gcs    = row["R1_cloud_filepath"]
    r2_gcs    = row["R2_cloud_filepath"]

    r1_local    = f"{tmp_dir}/{sample_id}_R1.fastq.gz"
    r2_local    = f"{tmp_dir}/{sample_id}_R2.fastq.gz"
    merged      = f"{seq_dir}/{sample_id}.fastq.gz"
    merged_tmp  = f"{tmp_dir}/{sample_id}_merged.fastq.gz"
    unmerged_r1 = f"{tmp_dir}/{sample_id}_R1u.fastq.gz"
    unmerged_r2 = f"{tmp_dir}/{sample_id}_R2u.fastq.gz"
    log_file    = f"{seq_dir}/{sample_id}_bbmerge.log"

    log(f"[{sample_id}] Downloading R1...")
    download(r1_gcs, r1_local, "R1", sample_id)
    log(f"[{sample_id}] Downloading R2...")
    download(r2_gcs, r2_local, "R2", sample_id)

    log(f"[{sample_id}] Running BBMerge...")
    result = subprocess.run(
        f"bbmerge.sh in={r1_local} in2={r2_local} "
        f"out={merged_tmp} outu1={unmerged_r1} outu2={unmerged_r2} "
        f"2>{log_file}",
        shell=True
    )

    if result.returncode != 0:
        log(f"[{sample_id}] BBMerge failed — falling back to raw R1+R2 concatenation")
        subprocess.run(f"cat {r1_local} {r2_local} > {merged}", shell=True, check=True)
    else:
        subprocess.run(
            f"cat {merged_tmp} {unmerged_r1} {unmerged_r2} > {merged}",
            shell=True, check=True
        )
        for f in [merged_tmp, unmerged_r1, unmerged_r2]:
            if os.path.exists(f):
                os.remove(f)

    os.remove(r1_local)
    os.remove(r2_local)

    n_reads = int(subprocess.check_output(f"zcat {merged} | wc -l", shell=True)) // 4
    log(f"[{sample_id}] {n_reads} reads written")

    if n_reads == 0:
        raise RuntimeError(f"{sample_id} produced 0 reads — check {log_file}")

    return sample_id, n_reads

with ThreadPoolExecutor(max_workers=16) as pool:
    futures = {pool.submit(process_sample, row): row["sample_ID"] for row in rows}
    failed  = []
    for future in as_completed(futures):
        sid = futures[future]
        try:
            future.result()
        except Exception as e:
            failed.append(sid)
            log(f"ERROR [{sid}]: {e}")

if failed:
    print(f"FAILED samples ({len(failed)}): {failed}", file=sys.stderr)
    sys.exit(1)

# Rewrite sample table with fastq_filepath column for phip-flow
with open("/phipflow/data/sample_table.csv", "w") as f:
    f.write("control_status,sample_ID,fastq_filepath\n")
    for row in rows:
        f.write(f"{row['control_status']},{row['sample_ID']},{row['sample_ID']}.fastq.gz\n")

print(f"All {len(rows)} samples ready.")
CODE

        df -h

        # ── Run phip-flow over the full cohort ─────────────────────────────────
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

        # ── Decompress results, upload, emit sentinel ──────────────────────────
        find /phipflow/results/ -name "*.gz" -exec gunzip {} \; 2>/dev/null || true

        gcloud storage cp -r /phipflow/results/ ~{output_directory}/
        gcloud storage cp /phipflow/data/seq/*_bbmerge.log ~{output_directory}/bbmerge_logs/ || true

        touch "~{output_prefix}_outs.txt"
        gcloud storage cp "~{output_prefix}_outs.txt" ~{output_directory}/
    >>>

    runtime {
        bootDiskSizeGb: 50
        disks:          "local-disk ~{disk_space} HDD"
        docker:         "~{docker_registry}"
        cpu:            num_cpu
        zone:           zone
        memory:         memory
        preemptible:    preemptible
    }
}
