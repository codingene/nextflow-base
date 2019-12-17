#!/usr/bin/env nextflow

/*
========================================================================================
                         nextflow-base
========================================================================================
 Writtrn by - Sangram keshari Sahu
 Basic nextflow exicution for a fastq files
----------------------------------------------------------------------------------------
*/


params.reads = "$baseDir/test_data/FASTQ/SRR*_{1,2}.fastq.gz"
params.transcriptome = "$baseDir/test_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.cdna.all.1.1.10M.fa"
params.outdir = "$baseDir/results"

transcriptome_file = file(params.transcriptome)
res_dir = file(params.outdir)

Channel
    .fromFilePairs( params.reads, checkExists:true )
    .into { read_pairs_ch; read_pairs2_ch }

process fastp {
    tag "FASTP on $pair_id"
    publishDir params.outdir, mode: 'copy'

    input:
    set pair_id, file(reads) from read_pairs_ch

    output:
    file("filtred_${pair_id}/${pair_id}_{1,2}.fastq.gz") into fastp_ch
    
    script:
    """
    mkdir -p filtred_${pair_id}
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o filtred_${pair_id}/${reads[0]} \
        -O filtred_${pair_id}/${reads[1]} \
        --json filtred_${pair_id}/${pair_id}_fastp.json \
		--html filtred_${pair_id}/${pair_id}_fastp.html \
        --detect_adapter_for_pe \
        --disable_length_filtering \
        --correction
    """ 
}

process fastqc {
    tag "FASTQC on $sample_id"

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """ 
}

process index {
    tag "$transcriptome.simpleName"

    input:
    file transcriptome from transcriptome_file

    output:
    file 'index' into index_ch

    script:
    """
    kallisto index -i index $transcriptome 
    """
}

process quant {
    tag "Kallisto quant on $filtred_pair_id"
    publishDir params.outdir, mode: 'copy'

    input:
    file index from index_ch
    set filtred_pair_id, file(pair_id) from fastp_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    kallisto quant -i $index ${filtred_pair_id[0]} ${filtred_pair_id[1]} -o $pair_id
    """
}

process multiqc {

    input:
    file('*') from quant_ch.mix(fastqc_ch).collect()
    
    output:
    file('multiqc_report.html') optional true

    script:
    """
    multiqc . -f -o $res_dir
    """
    }

workflow.onComplete {
    
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $res_dir/multiqc_report.html\n" : "Oops .. something went wrong" )
}