# Whole Genome Analysis

For whole-genome sequencing studies it is highly recommended to split the analysis by chromosome and then merge the results.
Splitting the work by chromosome will reduce overall runtime requirements to analyze whole-genome data.

### SLURM

```bash
NUM_CORES=64

for chrom in $(head -24 GRCh38.fasta.fai | cut -f1 | tr '\n' ' ')
do
sbatch --job-name="Lancet2_${chrom}" --cpus-per-task="${NUM_CORES}" --wrap \
"Lancet2 pipeline --num-threads ${NUM_CORES} \
    --normal normal.bam --tumor tumor.bam \
    --reference GRCh38.fasta --region ${chrom} \
    --out-vcfgz output.${chrom}.vcf.gz"
done
```

### SGE (Sun Grid Engine)

```bash
NUM_CORES=64

for chrom in $(head -24 GRCh38.fasta.fai | cut -f1 | tr '\n' ' ')
do
qsub -N "Lancet2_${chrom}" -cwd -pe smp "${NUM_CORES}" -j y -b y \
"Lancet2 pipeline --num-threads ${NUM_CORES} \
    --normal normal.bam --tumor tumor.bam \
    --reference GRCh38.fasta --region ${chrom} \
    --out-vcfgz output.${chrom}.vcf.gz"
done
```

### Merging results

After all per-chromosome jobs complete, merge the VCF files:

```bash
bcftools concat -Oz -o merged.vcf.gz output.chr*.vcf.gz
bcftools index -t merged.vcf.gz
```
