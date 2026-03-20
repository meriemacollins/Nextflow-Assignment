#!/usr/bin/env nextflow

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
