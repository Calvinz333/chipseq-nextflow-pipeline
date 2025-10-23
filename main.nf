#!/usr/bin/env nextflow
/*
 * ChIP-seq pipeline: (Local Files Version)
 * FastQC → BWA → Samtools → MACS2 → deepTools
 */

nextflow.enable.dsl=2

// -------------------------------
// PARAMETERS
// -------------------------------
params.reads = params.reads ?: "$baseDir/dataa/*_{1,2}.fastq.gz"
params.genome_fasta = params.genome_fasta ?: "$baseDir/genomee/hg38.fa"
params.genome_size = params.genome_size ?: 'hs'
params.outdir = params.outdir ?: 'results'

// -------------------------------
// CHANNELS
// -------------------------------

// Create a 6-part genome channel from the fasta param
// This makes the tuple (fa, amb, ann, bwt, pac, sa) that ALIGN expects
genome_ch = Channel.fromPath(params.genome_fasta)
                 .map { fa -> tuple(fa, 
                              file(fa.toString() + ".amb"), 
                              file(fa.toString() + ".ann"), 
                              file(fa.toString() + ".bwt"), 
                              file(fa.toString() + ".pac"), 
                              file(fa.toString() + ".sa")) 
                      }

// Create a 3-part read channel from the reads param
// This makes the tuple (sample_id, read1_path, read2_path)
read_ch = Channel.fromFilePairs(params.reads, size: 2)
                 .map { sample_id, files -> tuple(sample_id, files[0], files[1]) }

// -------------------------------
// PROCESSES
// -------------------------------

// 1️⃣ FASTQC
process FASTQC {
    tag "$accession"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(accession), path(read1), path(read2)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")

    script:
    """
    fastqc ${read1} ${read2} --outdir .
    """
}

// 2️⃣ ALIGNMENT — BWA MEM
process ALIGN {
    tag "$accession"

    input:
    tuple val(accession), path(read1), path(read2)
    tuple path(genome_fasta), path(amb), path(ann), path(bwt), path(pac), path(sa)

    output:
    tuple val(accession), path("${accession}.bam")

    script:
    """
    bwa mem -t ${task.cpus} ${genome_fasta} ${read1} ${read2} \
        | samtools view -Sb - \
        | samtools sort -o ${accession}.bam
    samtools index ${accession}.bam
    """
}

// 3️⃣ FILTERING — remove low-quality reads
process FILTER {
    tag "$accession"

    input:
    tuple val(accession), path(bam)

    output:
    tuple val(accession), path("${accession}.filtered.bam"), path("${accession}.filtered.bam.bai")
    script:
    """
    samtools view -b -q 30 ${bam} > ${accession}.filtered.bam
    samtools index ${accession}.filtered.bam
    """
}

// 4️⃣ PEAK CALLING — MACS2
process PEAKCALLING {
    tag "$accession"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(accession), path(filtered_bam), path(index)

    output:
    path("${accession}_peaks.narrowPeak")
    path("${accession}_summits.bed")

    script:
    """
    macs2 callpeak -t ${filtered_bam} -n ${accession} -g ${params.genome_size} \
        --outdir . --nomodel --shift -100 --extsize 200
    """
}

// 5️⃣ VISUALIZATION — deepTools
process VISUALIZATION {
    tag "$accession"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(accession), path(filtered_bam), path(index)

    output:
    path("${accession}.bw")
    path("profile_${accession}.png")

    script:
    """
    bamCoverage -b ${filtered_bam} -o ${accession}.bw --normalizeUsing CPM
    plotFingerprint -b ${filtered_bam} -plot profile_${accession}.png
    """
}

// -------------------------------
// WORKFLOW
// -------------------------------
workflow {
    FASTQC(read_ch)
    
    aligned = ALIGN(read_ch, genome_ch)
    filtered = FILTER(aligned)
    peaks = PEAKCALLING(filtered)
    VISUALIZATION(filtered)
}