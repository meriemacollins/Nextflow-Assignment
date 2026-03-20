#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/*_{1,2}.fq'
params.outdir = 'outputs/'
params.adapters = 'adapters.fa'
params.genome = 'extras/LG12.fasta'

log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
Genome           : ${params.genome}
"""

// Create read channel
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect { it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)
genome_ch = Channel.fromPath(params.genome)

// Define fastqc process
process fastqc {
    publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    fastqc ${reads}
    """
}

// Process trimmomatic
process trimmomatic {
    publishDir "${params.outdir}/trimmed-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path adapters_file

    output:
    tuple val("${sample}"), path("${sample}*.trimmed.fq.gz"), emit: trimmed_fq
    tuple val("${sample}"), path("${sample}*.discarded.fq.gz"), emit: discarded_fq

    script:
    """
    trimmomatic PE -phred33 ${reads[0]} ${reads[1]} ${sample}_1.trimmed.fq.gz ${sample}_1.discarded.fq.gz ${sample}_2.trimmed.fq.gz ${sample}_2.discarded.fq.gz ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

// Process 
process bwa_mem2 {
    publishDir "${params.outdir}/alignment-${sample}/", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), path(reads)
    path genome

    output:
    tuple val(sample), path("${sample}.sorted.bam")

   script:
    """
    bwa_mem2 mem -t $task.cpus ${genome} ${reads[0]} ${reads[1]} | samtools sort --threads $task.cpus -o ${sample}.sorted.bam - 
    """
}

// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
    bwa_mem2(trimmomatic.out.trimmed_fq, genome_ch)
}

