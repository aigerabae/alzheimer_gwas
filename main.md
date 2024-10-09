# Alzheimer idats to vcf:
dragena genotype call \
    --idat-sample-sheet /home/user/biostar/gwas/alzheimer/sample_sheet/SampleSheet_alz.csv \
    --bpm-manifest /home/user/biostar/gwas/alzheimer/manifest_bpm/InfiniumImmunoArray-24v2-0_A.bpm    \
    --cluster-file /home/user/biostar/gwas/alzheimer/cluster/InfiniumImmunoArray-24v2-0_A_ClusterFile.egt \    
    --output-folder /home/user/biostar/gwas/alzheimer/gtc/

dragena genotype gtc-to-vcf \
    --bpm-manifest /home/user/biostar/gwas/alzheimer/manifest_bpm/InfiniumImmunoArray-24v2-0_A.bpm \
    --genome-fasta-file /home/user/biostar/gwas/redo_october/making_vcf/GRCh37_genome/GRCh37_genome.fa \
    --gtc-sample-sheet /home/user/biostar/gwas/alzheimer/sample_sheet/SampleSheet_alz.csv \
    --csv-manifest /home/user/biostar/gwas/alzheimer/manifest_csv/InfiniumImmunoArray-24v2-0_A.csv \
    --output-folder /home/user/biostar/gwas/alzheimer/vcf/
