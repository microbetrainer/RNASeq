ALL_SAMPLES = glob_wildcards("data/raw_reads/{sample}.fastq.gz").sample

# output rule
rule target:
    input:
        expand("data/sense_counts/{sample}_noadapt_qtrim_mapped_sense.tab",sample=ALL_SAMPLES),expand("data/antisense_counts/{sample}_noadapt_qtrim_mapped_antisense.tab",sample=ALL_SAMPLES)


# trim adapters and linker sequences
rule run_adapttrim:
    input:
        reads = "data/raw_reads/{sample}.fastq.gz",
        adapter = "data/genomes/adapters.fa"
    output:
        "tmp/noadapt/{sample}_noadapt.fastq.gz"
    conda:
        "envs/bbtools.yaml"
    shell:
        "bbduk.sh in={input.reads} out={output} ref={input.adapter} ktrim=r k=23 mink=11 hdist=1" 

# quality trim
rule run_qualitytrim:
    input:
        "tmp/noadapt/{sample}.fastq.gz"
    output:
        "tmp/qualtrim/{sample}_qtrim.fastq.gz"
    conda:
        "envs/bbtools.yaml"
    shell:
        "bbduk.sh in={input} out={output} qtrim=r trimq=10"
   
# align reads to genome with BWA
rule bwa_map:
    input:
        "data/genomes/NATL2A_MIT1002.fna",
        "tmp/qualtrim/{sample}.fastq.gz"
    output:
        "tmp/mapped_reads/{sample}_mapped.bam"
    conda:
        "envs/bwa.yaml"
    shell: 
        "bwa mem {input} | samtools view -S -b > {output}"
 
# count sense reads with htseq-count
rule run_htseq:
    input:
        "tmp/mapped_reads/{sample}.bam",
        "data/genomes/NATL2A_MIT1002_mod.gff"
    output:
        "data/sense_counts/{sample}_sense.tab"
    conda:
        "envs/htseq.yaml"
    shell:
        "htseq-count -f bam -i ID -t CDS -s yes {input} > {output}"

# count antisense reads with htseq-count    
rule run_htseq_as:    
    input:
        "tmp/mapped_reads/{sample}.bam",
        "data/genomes/NATL2A_MIT1002_mod.gff"
    output:
        "data/antisense_counts/{sample}_antisense.tab"
    conda:
        "envs/htseq.yaml"
    shell:
        "htseq-count -f bam -i ID -t CDS -s reverse {input} > {output}"

