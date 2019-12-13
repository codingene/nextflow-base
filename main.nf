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
params.outdir = "$baseDir/results"

Channel
    .fromFilePairs( params.reads, checkExists:true )
    .into { read_pairs_ch; read_pairs2_ch }

process fastqc {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

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

process multiqc {
    
    output:
    file('multiqc_report.html') optional true

    script:
    """
    multiqc ${params.outdir} -o ${params.outdir}
    """
}
