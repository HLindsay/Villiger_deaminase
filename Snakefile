MOUSE_DIR = "annotation/Ensembl_GRCm38.90"
MOUSE_IDX = expand("{mdir}/bwaIndex/Mus_musculus.GRCm38.dna.primary_assembly",
                   mdir = MOUSE_DIR)

CONDITIONS, SAMPLES = glob_wildcards("fastq/{cond}/{sample}_L001_R1_001.fastq.gz")
ALL_MERGED, = glob_wildcards("merged_fastq/{file}.fastq")


rule all:
  input:
    expand("bam/{cond}/{sample}.bam",
           zip, cond = CONDITIONS, sample = SAMPLES)


rule merge_pairs:
  input:
    fwd = "fastq/{cond}/{sample}_L001_R1_001.fastq.gz",
    rev = "fastq/{cond}/{sample}_L001_R2_001.fastq.gz",
  params:
    outnm = "merged_fastq/{cond}/{sample}"
  threads: 20
  output: "merged_fastq/{cond}/{sample}.assembled.fastq.gz"
  shell:  
    """ pear -f {input.fwd} -r {input.rev} -o {params.outnm} -j {threads};
        gzip {params.outnm}*"""


# Note: Only mapping the successfully assembled pairs
rule map_merged:
  input:
    fq = "merged_fastq/{cond}/{sample}.assembled.fastq.gz"
  threads: 20
  params:
    idx = MOUSE_IDX,
    tempnm = "{cond}_{sample}",
  output: "bam/{cond}/{sample}.bam"
  shell: "bwa mem -t {threads} {params.idx} {input.fq} | samtools view -Sb - > {params.tempnm}_temp.bam; "
         "samtools sort {params.tempnm}_temp.bam -o {output} && rm {params.tempnm}_temp.bam"
