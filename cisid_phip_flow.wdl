version 1.0

workflow phippery_flow{
    input {
        File input_sample_table
        File input_peptide_table
        String output_directory
        Int read_length = 125
        Int oligo_tile_length = 117
        Int n_mismatches = 2
        String output_prefix = "data"
        Boolean replicate_sequence_counts = true
        Boolean run_beer
        Boolean run_cpm_enrichment
        Boolean run_z_score
        String zone = "us-central1-d"
        String memory = "32G"
        Int num_cpu = 8
        Int preemptible = 2
        Int disk_space = 250
        String docker_registry
    }

    call run_phippery_flow {
        input:
            input_sample_table = input_sample_table,
            input_peptide_table = input_peptide_table,
            output_directory = output_directory,
            read_length = read_length,
            oligo_tile_length = oligo_tile_length,
            n_mismatches = n_mismatches,
            output_prefix = output_prefix,
            replicate_sequence_counts = replicate_sequence_counts,
            run_beer = run_beer,
            run_cpm_enrichment = run_cpm_enrichment,
            run_z_score = run_z_score,
            zone = zone,
            memory = memory,
            num_cpu = num_cpu,
            preemptible = preemptible,
            disk_space = disk_space,
            docker_registry = docker_registry
    }

    output {
        String outs = run_phippery_flow.outs_files
    }
}

task run_phippery_flow{
    input {
        File input_sample_table
        File input_peptide_table
        String output_directory
        Int read_length = 125
        Int oligo_tile_length = 117
        Int n_mismatches = 2
        String output_prefix = "data"
        Boolean replicate_sequence_counts = true
        Boolean run_beer = false
        Boolean run_cpm_enrichment = true
        Boolean run_z_score = true
        String zone
        String memory
        Int num_cpu
        Int preemptible
        Int disk_space
        String docker_registry
    }

    output {
        String outs_files = "${output_directory}/${output_prefix}_outs.txt"
    }

    command <<<
        #!/bin/bash
        set -euo pipefail

        cd /phipflow/

        cp ~{input_sample_table} /phipflow/data/sample_table.csv
        cp ~{input_peptide_table} /phipflow/data/peptide_table.csv

        mkdir -p /phipflow/data/seq /phipflow/data/tmp

        python3 <<CODE
import pandas as pd
import os
import subprocess
import sys

df = pd.read_csv("/phipflow/data/sample_table.csv")

required_cols = {"sample_ID", "R1_cloud_filepath", "R2_cloud_filepath"}
missing = required_cols - set(df.columns)
if missing:
    print(f"ERROR: sample table is missing columns: {missing}", file=sys.stderr)
    sys.exit(1)

seq_dir = "/phipflow/data/seq"
tmp_dir = "/phipflow/data/tmp"

for _, row in df.iterrows():
    sample_id = str(row["sample_ID"])
    r1_gcs    = row["R1_cloud_filepath"]
    r2_gcs    = row["R2_cloud_filepath"]

    r1_local  = f"{tmp_dir}/{sample_id}_R1.fastq.gz"
    r2_local  = f"{tmp_dir}/{sample_id}_R2.fastq.gz"
    merged    = f"{seq_dir}/{sample_id}.fastq.gz"
    log_file  = f"{seq_dir}/{sample_id}_bbmerge.log"

    print(f"[{sample_id}] Downloading R1...")
    subprocess.run(["gcloud", "storage", "cp", r1_gcs, r1_local], check=True)

    print(f"[{sample_id}] Downloading R2...")
    subprocess.run(["gcloud", "storage", "cp", r2_gcs, r2_local], check=True)

    print(f"[{sample_id}] Merging R1+R2 with BBMerge...")
    result = subprocess.run(
        f"bbmerge.sh in={r1_local} in2={r2_local} out={merged} outu=/dev/null 2>{log_file}",
        shell=True
    )
    if result.returncode != 0:
        print(f"  WARNING: BBMerge returned non-zero for {sample_id} — check {log_file}")

    n_merged = int(subprocess.check_output(f"zcat {merged} | wc -l", shell=True)) // 4
    print(f"[{sample_id}] {n_merged} merged reads written to {merged}")

    os.remove(r1_local)
    os.remove(r2_local)

print("All samples downloaded and merged.")
CODE

        df -h

        CMD="nextflow run main.nf --ansi-log false"
        CMD="$CMD --sample_table /phipflow/data/sample_table.csv"
        CMD="$CMD --peptide_table /phipflow/data/peptide_table.csv"
        CMD="$CMD --reads_prefix /phipflow"

        # Read length (only pass if non-default)
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

        gcloud storage cp -r /phipflow/results/ ~{output_directory}

        touch "~{output_directory}/~{output_prefix}_outs.txt"
    >>>

    runtime {
        bootDiskSizeGb: 250
        disks: "local-disk ${disk_space} HDD"
        docker: "${docker_registry}"
        cpu: num_cpu
        zone: zone
        memory: memory
        preemptible: preemptible
    }
}
