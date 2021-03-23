# Smart-seq-total
This repository contains a Snakemake pipeline and R notebooks used to analyze Smart-seq-total v1 and v2 data.

![github](https://user-images.githubusercontent.com/40051862/112074504-abfcda80-8b33-11eb-86d7-45215fd1291a.png)


## Prerequisites

snakemake (https://snakemake.readthedocs.io/en/stable/)

cutadapt (https://cutadapt.readthedocs.io/en/stable/)

STAR (https://github.com/alexdobin/STAR)

htseq (https://htseq.readthedocs.io/en/master/) 

samtools (http://www.htslib.org) 

featureCounts (http://subread.sourceforge.net) 

UMI-tools (https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)

R (https://www.r-project.org)

*all can be installed through conda


# Smart-seq-total.v1

To analyse the data generated through Smart-seq-total.v1 (i.e. no UMIs), for each cell we follow a simple trim -> map -> count procedure.

We use Snakemake workflow manager to iterate through all the cells. Snakefile.v1 contains all the rules needed to generate the count files.
To launch the workflow on a Slurm operated cluster:

           snakemake all --snakefile Snakefile.v1 --cluster "sbatch --ntasks=1 --time=2:00:00 --job-name={params.name} \
                         --cpus-per-task={threads} --partition="normal" --mem={params.mem} -o {params.name}.%j.log" \
                         --keep-target-files -j 100 -w 100 -k -p -r --rerun-incomplete

## Smart-seq-total.v1 pipeline step-by-step:


STEP 1. Trim A-tails (>=10bp) from the 3'prime of the read

           cutadapt -j 4 -a AAAAAAAAAA -m 18 -o {output.trimmed.R1} {input.R1}


STEP 2. Map reads to the genome

           STAR        --genomeDir {starreference.genome} \
                       --outFileNamePrefix {cell.name.dir}. \
                       --readFilesIn {output.trimmed.R1}\
                       --runThreadN 4 \
                       --readFilesCommand zcat \
                       --outReadsUnmapped Fastx \
                       --outSAMtype BAM Unsorted SortedByCoordinate \
                       --outSAMattributes All \
                       --outFilterMismatchNoverLmax 0.05 \
                       --outFilterMatchNmin 16 \
                       --outFilterScoreMinOverLread 0 \
                       --outFilterMatchNminOverLread 0 \
                       --outMultimapperOrder Random
 
STEP 3. Count genes using either htseq or featureCounts

##### htseq ########
        samtools view -h {input.bam} | htseq-count \
                      --mode intersection-strict \
                      --stranded no \
                      --nonunique all \
                      --additional-attr gene_name \
                      --idattr Name - {input.gtf.file} > {output.counts}
  

##### featureCounts #####
         featureCounts -a {input.gtf.file} \
                       -M --primary \       ### this is otpional, if you want to also count the primary alignment of multi-mapping reads.
                       -o {output.counts} \
                       -T 4 \
                       {input.bam}





# Smart-seq-total.v2

To analyse the data generated through Smart-seq-total.v2 (with UMIs), we use UMI-tools to extract UMIs from Read1 and append this information to Read2. Read2 is then trimmed and mapped to the genome. We them count the number of reads per gene based on the mapping co-ordinate and the UMI attached to the read.
Snakefile.v2 contains all the rules needed to generate the count files.

To launch the workflow on a Slurm operated cluster:

    snakemake all --snakefile Snakefile.v2 --cluster "sbatch --ntasks=1 --time=2:00:00 --job-name={params.name} \
                     --cpus-per-task={threads} --partition="normal" --mem={params.mem} -o {params.name}.%j.log" \
                     --keep-target-files -j 100 -w 100 -k -p -r --rerun-incomplete
                         

## Smart-seq-total.v2 pipeline step-by-step:

STEP 1: Extract UMIs from Read1 and append them to Read2 headers 
    
    umi_tools extract --stdin {input.R1} --bc-pattern=NNNNNNNNNNNNNNNN \
                      --log {output.log}  --read2-in {input.R2} --read2-out {output}

  ** R1 is usually a mix of: 1) reads containing a UMI followed by a stretch of Ts (~70%) and 2) gene body reads corresponding to tagmented molecules (~30%). For precise molecule count, one should filter R1 files (and R2 respectively) and keep only UMI-containing pairs (based on the presence of TTTTTTT at the end of each R1) (type 1 read above). Otherwise, assuming that each mRNA molecule gets no more than two cuts during tagmentation, the gene body part of R1 (type 2 read) can be also treated as a UMI (since two mRNA molecules are unlikely to be cut at the exact same location).
  
STEP 2: Trim A-tails (>=10bp) from the 3'prime of the Read2

     cutadapt -j 4 -u 6 -a AAAAAAAAAA -m 18 -o {output.trimmed.R2} {input.R2}

STEP 3: Trim A-tails (>=10bp) from the 3'prime of the Read2

      STAR        --genomeDir {star.reference.genome} \
                  --outFileNamePrefix {cell.name.dir}. \
                  --readFilesIn {trimmed.R2}\
                  --runThreadN 4 \
                  --readFilesCommand zcat \
                  --outReadsUnmapped Fastx \
                  --outSAMtype BAM Unsorted SortedByCoordinate \
                  --outSAMattributes All \
                  --outFilterMismatchNoverLmax 0.05 \
                  --outFilterMatchNmin 16 \
                  --outFilterScoreMinOverLread 0 \
                  --outFilterMatchNminOverLread 0 \
                  --outMultimapperOrder Random


STEP 4: Assign reads to features
   
    featureCounts     -a {input.gtf} \
                      -o {output.fc} \
                      -R BAM {input.bam} \
                      -T 4
   
    samtools sort {fc.bam} -o {output.sorted.bam}
    samtools index {output.sorted.bam}      
                      
                      
STEP 5: Count features (genes)                     
       
    umi_tools count --per-gene \
                    --gene-tag=XT \
                    --assigned-status-tag=XS -I {sorted.indexed.bam} -S {output.counts}

