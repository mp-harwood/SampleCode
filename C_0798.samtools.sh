#!/bin/sh
#SBATCH --mem=10000M
#SBATCH --time=00-01:00

module load samtools
        
samtools view star_outfiles/C_0798/C_0798Aligned.out.bam |awk -F\t '$0 ~"vW:i:1"' > star_outfiles/C_0798/C_0798Aligned_WASP_filter_PASS.sam
samtools view -H star_outfiles/C_0798/C_0798Aligned.out.bam > header
cat header star_outfiles/C_0798/C_0798Aligned_WASP_filter_PASS.sam > star_outfiles/C_0798/C_0798Aligned_WASP_filter_PASS_header.sam
samtools view -S -b star_outfiles/C_0798/C_0798Aligned_WASP_filter_PASS_header.sam | samtools sort -l1 -o star_outfiles/C_0798/C_0798Aligned_WASP_filter_PASS_sorted.bam
