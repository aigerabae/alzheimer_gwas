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

awk -F'\t' '!seen[$1]++' InfiniumImmunoArray-24v2-0_A_b138_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > IIA-dictionary.txt
plink --bfile alz2 --update-name IIA-dictionary.txt --make-bed --out alz3
awk '$2 !~ /^rs/' alz3.bim | sort -k2,2 > custom_non_rs_SNP.txt
plink --bfile alz3 --exclude custom_non_rs_SNP.txt --make-bed --out alz4

plink --bfile alz4 --geno 0.02 --make-bed --out alz5
plink --bfile alz5 --mind 0.02 --make-bed --out alz6
plink --bfile alz6 --maf 0.001 --make-bed --out alz7
plink --bfile alz7 --genome --min 0.2 --out pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile alz7 --missing --out missing_report
echo "D180    D180
AK015    AK015
C132    C132
C124    C124
C077    C077
C067    C067
C170    C170" > relatives_to_remove.tsv
plink --bfile alz7 --remove relatives_to_remove.tsv --allow-no-sex --make-bed --out alz8
plink2 --bfile alz8 --pca 10 --out alz_pca
```

```bash
plink --bfile alz8 --covar alz_pca.eigenvec --logistic --hide-covar --thread-num 8 --out simple_logistic
plink --bfile alz8 --covar alz_pca.eigenvec --logistic --dominant --hide-covar --out dominant_results
plink --bfile alz8 --covar alz_pca.eigenvec --logistic --recessive --hide-covar --out recessive_results
plink --bfile alz8 --assoc --out assoc_results
plink --bfile alz8 --fisher --out fisher

plink --bfile alz8 --covar alz_pca.eigenvec  --model --out model

plink2 --bfile alz8 --glm --covar alz_pca.eigenvec --out glm
```

```bash
cat dominant_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat recessive_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat simple_logistic.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat assoc_results.assoc | awk '$9 != "NA"' | sort -gk 9,9 | head
cat fisher.assoc.fisher | awk '$8 != "NA"' | sort -gk 8,8 | head
cat model.model | awk '$10 != "NA"' | sort -gk 10,10 | head
cat glm.PHENO1.glm.logistic.hybrid | awk '$15 != "NA"' | sort -gk 15,15 | head
```

```bash
./manhattan_plot.py -i assoc_results.assoc 
./manhattan_plot.py -i dominant_results.assoc.logistic 
./manhattan_plot.py -i fisher.assoc.fisher
./manhattan_plot.py -i glm.PHENO1.glm.logistic.hybrid
./manhattan_plot.py -i model.model
./manhattan_plot.py -i recessive_results.assoc.logistic
./manhattan_plot.py -i simple_logistic.assoc.logistic
```

```bash
cat simple_logistic.assoc.logistic | awk '$9 != "NA" && $9 < 1e-3' | sort -gk 9,9
cat dominant_results.assoc.logistic | awk '$9 != "NA" && $9 < 1e-3' | sort -gk 9,9
cat recessive_results.assoc.logistic | awk '$9 != "NA" && $9 < 1e-3' | sort -gk 9,9
cat assoc_results.assoc | awk '$9 != "NA" && $9 < 1e-3' | sort -gk 9,9
cat fisher.assoc.fisher | awk '$9 != "NA" && $9 < 1e-3' | sort -gk 9,9
# too many values
cat model.model | awk '$10 != "NA" && $10 < 1e-4' | sort -gk 9,9

I WILL DO IT TOMORROW MONING
```

plink2:
```bash
plink2 --bfile alz8 --glm --covar alz_pca.eigenvec --pfilter 0.0001 --out glm_additive
plink2 --bfile alz8 --glm dominant --covar alz_pca.eigenvec --pfilter 0.0001 --out glm_dominant
plink2 --bfile alz8 --glm recessive --covar alz_pca.eigenvec --pfilter 0.0001 --out glm_recessive
plink2 --bfile alz8 --glm firth --covar alz_pca.eigenvec --pfilter 0.0001 --out glm_firth

plink2 --bfile alz8 --glm --covar alz_pca.eigenvec -out glm_additive
plink2 --bfile alz8 --glm dominant --covar alz_pca.eigenvec --out glm_dominant
plink2 --bfile alz8 --glm recessive --covar alz_pca.eigenvec  --out glm_recessive
plink2 --bfile alz8 --glm firth --covar alz_pca.eigenvec --out glm_firth

./manhattan_plot.py -i glm_additive.PHENO1.glm.logistic.hybrid
./manhattan_plot.py -i glm_dominant.PHENO1.glm.logistic.hybrid
./manhattan_plot.py -i glm_recessive.PHENO1.glm.logistic.hybrid
./manhattan_plot.py -i glm_firth.PHENO1.glm.firth
