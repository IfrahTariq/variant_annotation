version 1.0


workflow run_vcf_to_csv {
    input {
        String sample_id
        File input_vep_vcf
        File input_funcotator_maf
        String version=""
        String docker_image="itariq/variant_annotation:sha256:975132c4fd7676e02a85c3c183249107ff142d0d88c161dfbc4cfbe15f5ce02e"
    }

    call vcf_to_csv {
        input:
            input_vep_vcf=input_vep_vcf,
            input_funcotator_maf=input_funcotator_maf,
            sample_id=sample_id,
            version=version,
            docker_image=docker_image,
    }

    output {
        File csvfile=vcf_to_csv.csvfile
    }
}

# transforms a vcf file (if possible annotated with opencravat) into:
# 1. a csv file that does some basic cleanup from vcf
task vcf_to_csv {
    input {
        File input_vep_vcf
        File input_funcotator_maf
        String sample_id
        String version=""

        String docker_image="itariq/variant_annotation:sha256:975132c4fd7676e02a85c3c183249107ff142d0d88c161dfbc4cfbe15f5ce02e"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=40
        Int cpu = 4
        Int mem = 32
    }

    command {
        python -u  /process_vep_func/vcf_to_csv.py \
              ~{input_vep_vcf} \
              ~{input_funcotator_maf} \
              ~{sample_id} \
    }

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File csvfile = "~{sample_id}_vep_funcotator.csv"
    }
}
