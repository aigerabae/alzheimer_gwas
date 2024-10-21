```bash
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

bcftools merge /home/user/biostar/gwas/alzheimer/vcf/*.vcf.gz -o alz.vcf
bcftools reheader -s sample_name_conversion.tsv -o alz1.vcf alz.vcf
plink -vcf alz1.vcf --pheno phenotypes.tsv --update-sex sex.tsv --remove non-kz.tsv --make-bed --out alz2

plink --bfile alz2 --geno 0.02 --make-bed --out alz3
plink --bfile alz3 --mind 0.02 --make-bed --out alz4
plink --bfile alz4 --maf 0.001 --make-bed --out alz5
plink --bfile alz5 --genome --min 0.2 --out pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile alz10 --missing --out missing_report
echo "D180    D180
AK015    AK015
C132    C132
C124    C124
C077    C077
C067    C067
C170    C170" > relatives_to_remove.tsv
plink --bfile alz5 --remove relatives_to_remove.tsv --allow-no-sex --make-bed --out alz6
plink2 --bfile alz6 --pca 10 --out alz_pca
```

```bash
plink --bfile alz6 --covar alz_pca.eigenvec --logistic --hide-covar --thread-num 8 --out simple_logistic
plink --bfile alz6 --covar alz_pca.eigenvec --logistic --dominant --hide-covar --out dominant_results
plink --bfile alz6 --covar alz_pca.eigenvec --logistic --recessive --hide-covar --out recessive_results
plink --bfile alz6 --assoc --out assoc_results
plink --bfile alz6 --fisher
plink --bfile alz6 --model
plink2 --bfile alz6 --glm --covar alz_pca.eigenvec
```

```bash
cat dominant_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat recessive_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat simple_logistic.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat assoc_results.assoc | awk '$9 != "NA"' | sort -gk 9,9 | head
cat plink.assoc.fisher | awk '$8 != "NA"' | sort -gk 8,8 | head
cat plink.model | awk '$10 != "NA"' | sort -gk 10,10 | head
cat plink2.PHENO1.glm.logistic.hybrid | awk '$15 != "NA"' | sort -gk 15,15 | head
```

```bash

 ./manhattan_plot.py -i assoc_results.assoc 

./manhattan_plot.py -i dominant_results.assoc.logistic 

./manhattan_plot.py -i plink.assoc.fisher

./manhattan_plot.py -i plink2.PHENO1.glm.logistic.hybrid

./manhattan_plot.py -i recessive_results.assoc.logistic

./manhattan_plot.py -i simple_logistic.assoc.logistic

```

