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
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'

    input:
    file index from index_ch
    set pair_id, file(reads) from read_pairs_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    kallisto quant -i $index ${reads[0]} ${reads[1]} -o $pair_id
    """
}

process multiqc {

    input:
    file('*') from fastqc_ch.collect()
    
    output:
    file('multiqc_report.html') optional true

    script:
    """
    multiqc . --force -o $res_dir
    """
    }

workflow.onComplete {
    
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $res_dir/multiqc_report.html\n" : "Oops .. something went wrong" )
}