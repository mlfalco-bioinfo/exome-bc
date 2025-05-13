# Breast Cancer Exome Pipeline (Tumor-Only)

## 1. Quality Control (Falco and MultiQC)

### Step 1: Run Falco for FASTQ Quality Control

```bash
mkdir -p results/qc
ls *.fastq.gz | parallel -j 4 'falco {} -o results/qc'
```

* **mkdir -p results/qc:** Creates a directory for QC results (if it doesn't exist).
* **ls \*.fastq.gz:** Lists all FASTQ files in the directory.
* **parallel -j 4:** Runs 4 jobs in parallel (adjust based on CPU availability).
* **falco {} -o results/qc:** Runs Falco for QC on each FASTQ file, saving outputs in the `results/qc` directory.

  * **{}:** Placeholder for each file name.
  * **-o results/qc:** Specifies the output directory.

### Step 2: Summarize QC results using MultiQC

```bash
multiqc results/qc -o results/qc
```

* **results/qc:** The directory where individual QC reports are located.
* **-o results/qc:** Sets the output directory for the aggregated report.

## 2. Alignment (BWA-MEM, SAMtools, MarkDuplicates)

### Step 3: Prepare Reference Genome

```bash
bwa index hg38.fa
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
```

* **hg38.fa:** Path to the reference genome in FASTA format.
* **bwa index:** Builds an index for the reference genome.
* **samtools faidx:** Creates a FASTA index file.
* **gatk CreateSequenceDictionary:** Creates a dictionary file for GATK compatibility.

### Step 4: Align Reads to Reference

```bash
bwa mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" hg38.fa \
    ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz | samtools view -@4 -b - > results/align/${SAMPLE}.bam
```

* **bwa mem:** Aligns paired-end reads to the reference genome.
* **-t 8:** Uses 8 threads.
* **-R:** Read group information for GATK compatibility.

### Step 5: Sort and Index BAM File

```bash
samtools sort -@4 -o results/align/${SAMPLE}.sorted.bam results/align/${SAMPLE}.bam
samtools index results/align/${SAMPLE}.sorted.bam
```

* **samtools sort:** Sorts the BAM file.
* **samtools index:** Creates an index for the sorted BAM file.

### Step 6: Mark Duplicates

```bash
gatk MarkDuplicates -I results/align/${SAMPLE}.sorted.bam \
    -O results/align/${SAMPLE}.dedup.bam -M results/align/${SAMPLE}.metrics.txt
samtools index results/align/${SAMPLE}.dedup.bam
```

## 3. Variant Calling (GATK Mutect2)

### Step 7: Call Somatic Variants

```bash
gatk Mutect2 -R hg38.fa -I results/align/${SAMPLE}.dedup.bam \
    --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.hg38.vcf.gz \
    -O results/variants/${SAMPLE}.unfiltered.vcf.gz
```

### Step 8: Filter Variants

```bash
gatk FilterMutectCalls -R hg38.fa \
    -V results/variants/${SAMPLE}.unfiltered.vcf.gz \
    -O results/variants/${SAMPLE}.filtered.vcf.gz
```

## 4. Annotation (VEP)

### Step 9: Annotate Variants

```bash
vep -i results/variants/${SAMPLE}.filtered.vcf.gz \
    -o results/annotated/${SAMPLE}.vep.vcf \
    --cache --assembly GRCh38 --vcf --fork 4 \
    --dir_cache /path/to/vep_cache --everything
```

## 5. Copy Number Variant (CNV) Detection (CNVkit)

### Step 10: Call CNVs

```bash
cnvkit.py batch results/align/${SAMPLE}.dedup.bam \
    -n -t targets.bed -f hg38.fa -d results/cnv/
```

## 6. Parallelization (GNU Parallel)

### Step 11: Parallel Variant Calling

```bash
cat samples.txt | parallel -j 4 'gatk Mutect2 -R hg38.fa -I results/align/{}.dedup.bam --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.hg38.vcf.gz -O results/variants/{}.vcf.gz'
```
