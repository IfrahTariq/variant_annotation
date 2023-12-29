version 1.0


workflow HgvsWorkflow {
    input {
        File input_vcf
        String sample_id
        File pLi = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/pLI_values.txt"
        File LoF = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/LoFtool_scores.txt"

        File clinvar_data = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/clinvar_20231217.vcf.gz"
        File clinvar_data_tbi = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/clinvar_20231217.vcf.gz.tbi"
        File dbsnp_data = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/00-All.vcf.gz"
        File dbsnp_data_tbi = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/00-All.vcf.gz.tbi"

        File snpeff = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/snpEff.jar"
        File snpeff_config = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/snpEff.config"
        File snpsift = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/SnpSift.jar"

        File vep_data = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/homo_sapiens_vep_110_GRCh38.tar.gz"
        File fasta = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/Homo_sapiens_assembly38.fasta.gz"
        File fai = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/Homo_sapiens_assembly38.fasta.gz.fai"
        File gzi = "gs://fc-f28a7948-a3c6-48bb-a978-56732d4aa44d/Homo_sapiens_assembly38.fasta.gz.gzi"
        Int boot_disk_size=1000
        Int disk_space=500
        Int cpu = 60
        Int mem = 1000
    }

    call annotate_hgvs_task {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
            fasta=fasta,
            fai=fai,
            gzi=gzi,
            vep_data=vep_data,
            snpeff=snpeff,
            snpsift=snpsift,
            snpeff_config=snpeff_config,
            dbsnp_data=dbsnp_data,
            dbsnp_data_tbi=dbsnp_data_tbi,
            clinvar_data=clinvar_data,
            clinvar_data_tbi=clinvar_data_tbi,
            pLi=pLi,
            LoF=LoF,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu = cpu,
            mem = mem,
    }

    output {
        File maf=annotate_hgvs_task.output_maf   
        File vep=annotate_hgvs_task.output_vep_vcf
    }
}

# Standard interface to run vcf to maf
task annotate_hgvs_task {
    input {
        File input_vcf
        File fasta
        File fai
        File gzi
        File pLi
        File LoF
        File vep_data
        File snpeff
        File snpeff_config
        File snpsift
        File dbsnp_data
        File dbsnp_data_tbi
        File clinvar_data
        File clinvar_data_tbi
        String sample_id
        String docker_image="itariq/variant_annotation:sha256:061302d61c8ec1d316befb979774778fa44db6a3e7446fb98b7759196ef1468d"
        String assembly="GRCh38"
        Int preemptible=2
        Int boot_disk_size=1000
        Int disk_space=500
        Int cpu = 60
        Int mem = 1000
    }

    command {

        bcftools norm -m- -w 10000 -f ~{fasta} -O z -o ~{sample_id}.norm.vcf.gz ~{input_vcf}

        cp ~{snpeff_config} ~{snpeff} .

        java -Xmx8g -jar ~{snpeff} GRCh38.mane.1.0.ensembl ~{sample_id}.norm.vcf.gz > ~{sample_id}.norm.snpeff.vcf

        mkdir -p ./db/GRCh38/clinvar/
        cp ~{clinvar_data} ~{clinvar_data_tbi} ./db/GRCh38/clinvar/
        java -jar ~{snpsift} Annotate -clinvar -db ~{clinvar_data} ~{sample_id}.norm.snpeff.vcf > ~{sample_id}.norm.snpeff.clinvar.vcf

        tar -C /tmp -xvzf ~{vep_data} 
        ls /tmp
        chmod 777 /tmp/homo_sapiens
        #ls /tmp/homo_sapiens
        cp ~{fai} /tmp
        cp ~{fasta} /tmp
        cp ~{gzi} /tmp
        du -sh /tmp/Homo_sapiens_assembly38.fasta.gz*
       
        cp ~{pLi} ~{LoF} /tmp

        perl -e '(print -e "/tmp/homo_sapiens") && print "\n".$^V."\n"' && perl vep -id rs699 --cache --dir

        vep --species homo_sapiens --cache --assembly ~{assembly} --no_progress --no_stats --everything --dir /tmp --input_file ~{sample_id}.norm.snpeff.clinvar.vcf \
            --output_file ~{sample_id}.norm.snpeff.clinvar.vep.vcf \
            --plugin pLI,/tmp/pLI_values.txt --plugin LoFtool,/tmp/LoFtool_scores.txt \
            --force_overwrite --offline --fasta /tmp/Homo_sapiens_assembly38.fasta.gz --fork ~{cpu} --vcf \
            --pick 

        perl ../../vcf2maf/vcf2maf.pl \
            --input-vcf ~{sample_id}.norm.snpeff.clinvar.vep.vcf \
            --output-maf ~{sample_id}.maf \
            --ref /tmp/Homo_sapiens_assembly38.fasta.gz \
            --vep-path /opt/conda/envs/vep/bin/ \
            --vep-data /tmp/ \
            --ncbi-build ~{assembly} --inhibit-vep
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
        File output_maf = "~{sample_id}.maf"        
        File output_vep_vcf = "~{sample_id}.norm.snpeff.clinvar.vep.vcf"
    }
}

