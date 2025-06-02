version 1.0

workflow phippery_flow{
    input {
        File input_sample_table = "some/path/in/cloud"
        File input_peptide_table = "some/path/in/cloud"
        String input_fastq_dir = "some/path/in/cloud"
        String output_directory = "some/path/in/cloud"
        Int read_length = 125
        Int oligo_tile_length = 117
        Int n_mismatches = 2
        String output_prefix = "root/added/to/outs"
        Boolean replicate_sequence_counts = true
        Boolean run_beer = false
        Boolean run_cpm_enrichment = true
        Boolean run_z_score = true
        String zone = "us-central1-d"
        String memory = "32G"
        Int num_cpu = 8
        Int preemptible = 2
        Int disk_space = 250
        String docker_registry = "whereverimhosting"
    }

    call run_phippery_flow {
        input:
            input_sample_table = input_sample_table,
            input_peptide_table = input_peptide_table,
            input_fastq_dir = input_fastq_dir,
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
        File input_sample_table = "some/path/in/cloud"
        File input_peptide_table = "some/path/in/cloud"
        String input_fastq_dir = "some/path/in/cloud"
        String output_directory = "some/path/in/cloud"
        Int read_length = 125
        Int oligo_tile_length = 117
        Int n_mismatches = 2
        String output_prefix = "root/added/to/outs"
        Boolean replicate_sequence_counts = true
        Boolean run_beer = false
        Boolean run_cpm_enrichment = true
        Boolean run_z_score = true
        String zone
        String memory
        String num_cpu
        Int preemptible
        Int disk_space
        String docker_registry  
    }

    output {
        String outs_files = "${output_directory}/${output_prefix}_outs.txt"
    }

    command <<<
        #!/bin/bash
        set -e

        cp ~{input_sample_table} /phipflow/data/sample_table.csv
        cp ~{input_peptide_table} /phipflow/data/peptide_table.csv

        python <<CODE

        import pandas as pd
        import os
        import subprocess

        #Moving required R1 Fastqs to the specified run folder
        df = pd.read_csv("/phipflow/data/sample_table.csv")
        print(df['cloud_filepath'])

        for i, j in zip(list(df['sample_ID']), list(df['cloud_filepath'])):
            print(f"Moving {i} to /phipflow/data/seq/")
            subprocess.run(f"gcloud storage cp {j} /phipflow/data/seq/", shell=True)
        
        #Validation for error logs
        print(os.listdir("/phipflow/data/seq/"))

        CODE

        CMD = "nextflow run main.nf"

        # Read Length
        if [[read_length -ne 125 ]]; then
            CMD = "$CMD --read_length ~{read_length}"
        fi
        #Oligo_tile_length
        if [[oligo_tile_length -ne 117 ]]; then
            CMD = "$CMD --oligo_tile_length ~{oligo_tile_length}"
        fi
        #Number of Mismatches
        if [[n_mismatches -ne 2 ]]; then
            CMD = "$CMD --n_mismatches ~{n_mismatches}"
        fi
        #Optional Run parameters
        [[replicate_sequence_counts == true ]] && CMD= "$CMD --replicate_sequence_counts"
        [[run_beer == true ]] && CMD= "$CMD --run_BEER" 
        [[run_cpm_enrichment == true ]] && CMD="$CMD --run_cpm_enr_workflow"
        [[run_z_score == true ]] && CMD="$CMD --run_z_score_fit_predict"

        echo "Running: $CMD"
        eval $CMD

        #Transfer results back into the bucket

        gcloud storage cp -r /phipflow/results/ ~{output_directory}

        if [[ -f "${output_directory}/${output_prefix}_outs.txt" ]]; then
            echo "File exists"
        else
            touch "${output_directory}/${output_prefix}_outs.txt"  # Create empty fallback
        fi
    >>>

    runtime {
        preemptible : preemptible
        bootDiskSizeGb: 10
        disks: "local-disk ${disk_space} HDD"
        docker: "${docker_registry}"
        cpu: num_cpu
        zone: zone
        memory: memory
    }
}