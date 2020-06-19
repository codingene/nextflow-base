#!/usr/bin/env nextflow

/*
========================================================================================
                         nextflow-base
========================================================================================
 Writtrn by - Sangram keshari Sahu
 Basic Nextflow structure template
----------------------------------------------------------------------------------------
*/

def workflow_description() {
    return """
    ----------------------------------------------------------------------
    Workflow Name: ${workflow.manifest.name}-v${workflow.manifest.version}
    Short Description: ${workflow.manifest.description}
    ----------------------------------------------------------------------
    """.stripIndent()
}

def helpMessage() {
    log.info workflow_description()
    log.info"""
    Usage:
    nextflow run ${workflow.manifest.name}-v${workflow.manifest.version} [arguments]
    Workflow arguments:

      --reads                 Path to fastq files directory
      --cdna
      --outdir  

    Tool arguments:

      --fastp.length_required           (default: '${params.fastp.length_required}')
      --fastp.length_limit              (default: '${params.fastp.length_limit}')
      --fastp.qualified_quality_phred   (default: '${params.fastp.qualified_quality_phred}')

    System arguments:

      --max_cpus               Maximum CPU threads to be used (default: '${params.max_cpus}')
      --max_memory             Maximum memroy to be used (default: '${params.max_memory}')   
      --email                  Email ID for notofications (default: '${params.email}')
      --help                   This help menu

    Nextflow arguments:

      -profile                 Run execution settings [docker,conda,test]
      -resume                  Resume from where you left

Note: Please take care of single (-) and double (--) before an argument

Check documentation for more details on each argument.
----------------------------------------------------------------------
    """
}

// Help

if (params.help){
    helpMessage()
    exit 0
} else {
    log.info workflow_description()
}

// ######################### params summary ######################
def summary = [:]

summary['User']             = workflow.userName

summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir

summary['Input FASTA dir']  = params.reads
summary['Output dir']       = params.outdir

summary['Config Profile']   = workflow.profile
summary['Docker Container'] = workflow.container
summary['Maximum CPU']      = params.max_memory
summary['Maximum Memory']   = params.max_cpus

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")

def line() {
    return """
----------------------------------------------------------------------
""".stripIndent()
}

log.info line()

// ######################### Parameters check ######################

// Input
if(!params.reads){
    exit 1, "Please give a directory path containing pair-end fastq files."
} else {
    reads = "${params.reads}/*_{1,2}.f*.gz"
    Channel
    .fromFilePairs( reads, checkExists:true )
    .ifEmpty { error "Cannot find any reads matching: ${reads}" }
    .into { read_pairs_ch; read_pairs2_ch }
}

// Output directory
if(params.outdir) results_dir = params.outdir

// reference
if(params.cdna) cdna_file = file(params.cdna)

// ######################### Workflow Starting here ######################

process fastp {
    tag "filtering on $sample_id"
    publishDir "${results_dir}/${sample_id}", mode: 'copy'
    memory params.max_memory
    container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'

    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    set sample_id, file("fastp_filtred_reads/${reads[0]}"), file("fastp_filtred_reads/${reads[1]}") into fastp_ch, fastp_ch2
    set sample_id, file("fastp_filtred_reads/${sample_id}_fastp.json"), file("fastp_filtred_reads/${sample_id}_fastp.html") into fastp_reports

    script:
    """
    mkdir fastp_filtred_reads
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o fastp_filtred_reads/${reads[0]} \
        -O fastp_filtred_reads/${reads[1]} \
        --json fastp_filtred_reads/${sample_id}_fastp.json \
		--html fastp_filtred_reads/${sample_id}_fastp.html \
        --detect_adapter_for_pe \
        --disable_length_filtering \
        --correction \
        --thread ${task.cpus}
    """ 
}


process fastqc {
    tag "FASTQC on $sample_id"
    publishDir "${results_dir}/${sample_id}", mode: 'copy'
    cpus params.max_cpus
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    input:
    set sample_id, file(fq1), file(fq2) from fastp_ch

    output:
    file("fastqc_report_${sample_id}") into fastqc_report

    script:
    """
    mkdir fastqc_report_${sample_id}
    fastqc -f fastq -q ${fq1} ${fq2} \
        -o fastqc_report_${sample_id} \
        -f fastq -q ${fq1} ${fq2} \
        --threads ${task.cpus}
    """ 
}

process index {
    tag "$cdna.simpleName"
    container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'

    input:
    file cdna from cdna_file

    output:
    file 'index' into index_ch

    script:
    """
    kallisto index -i index $cdna 
    """
}

process quant {
    tag "Kallisto quant on $sample_id"
    publishDir "${results_dir}/${sample_id}", mode: 'copy'
    container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'

    input:
    file index from index_ch
    set sample_id, file(fq1), file(fq2) from fastp_ch2

    output:
    file("kallisto_quant_${sample_id}") into quant_ch

    script:
    """
    mkdir kallisto_quant_${sample_id}
    kallisto quant -i $index ${fq1} ${fq2} -o kallisto_quant_${sample_id} \
        &> kallisto_quant_${sample_id}/stdout.log
    """
}

process multiqc {
    container 'quay.io/biocontainers/multiqc:1.8--py_2'
    publishDir "${results_dir}", mode: 'copy'

    input:
    file('*') from quant_ch.mix(fastqc_report).collect()
    
    output:
    file('multiqc_report.html') optional true

    when: params.multiqc_report

    script:
    """
    multiqc . 
    """
    }

workflow.onComplete {
    
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $results_dir/multiqc_report.html\n" : "Oops .. something went wrong" )
}