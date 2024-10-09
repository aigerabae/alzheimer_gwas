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
```

Work I did on GenomeStudio generated FinalReport:
# Alzheimer idats to plink
```bash
awk 'NR > 1 { print $3, $2, 0, $4 }' SNP_Map.txt > alz.map
tail -n +11 "alz_FinalReport_2.txt" > alz1.txt
awk 'BEGIN {OFS="\t"} {$3=$3""$4; $4=""; print}' alz1.txt > alz2.txt
sed $'s/\r/\t/g' alz2.txt > alz3.txt
awk '{printf "%s%s", $3, (NR%253702 ? OFS : ORS)}' alz3.txt > pivoted_alz.txt
cat alz3.txt | cut -f 2 | uniq > pivoted_header.txt
awk 'BEGIN { OFS="\t" } { print 0, $0, 0, 0, 0, 0 }' pivoted_header.txt > temp_columns.txt
paste temp_columns.txt pivoted_alz.txt > alz.ped
```

Work I did on plink files that I created from FinalReport:
```bash
# Changing snp names to standard rsID names, excluding SNPs on 0 chromosome, adding sex and phenotypes:
plink --file alz --make-bed --out alz
awk -F'\t' '!seen[$1]++' InfiniumImmunoArray-24v2-0_A_b138_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > IA-dictionary.txt
plink --bfile alz --update-name IA-dictionary.txt --make-bed --out alz1
awk '$2 !~ /^rs/' alz1.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile alz1 --exclude non_rs_SNP.txt --make-bed --out alz2
awk '$1 == 0 {print $2}' alz2.bim > exclude_snps.txt
plink --bfile alz2 --exclude exclude_snps.txt --make-bed --out alz3
```

I used metadata (from different zapusks) to identify samples IDs and their matching Sentrix_ID positions; 
i wanted to manually replace 207851060016 with 207859430016 since they were the only samples that faild to replace names and they had an identical numebr of samples so i figured it could be a misspelling... but then i decided to ust keep them as they are since even if i did i would have 13 more samples with no phenotype data so i just excluded 29 samples with no phenotype data (among them 24 have no sample ID and 5 had a sample ID but were not present in the metadata).

I didn't remove other ethnicities 

```bash
plink --bfile alz3 --update-ids 'alzheimer_metadata_selected - true_dict.tsv' --make-bed --out alz4
plink --bfile alz4 --update-sex 'alzheimer_metadata_selected - sex.tsv' --make-bed --out alz5
plink --bfile alz5 --pheno 'alzheimer_metadata_selected - pheno.tsv' --make-bed --out alz6
awk '$6 == -9' alz6.fam | awk '{print $1"\t" $2}' > missing_phenotype.tsv
plink --bfile alz6 --remove missing_phenotype.tsv --allow-no-sex --make-bed --out alz7
plink --bfile alz7 --geno 0.02 --make-bed --out alz8
plink --bfile alz8 --mind 0.02 --make-bed --out alz9
plink --bfile alz9 --maf 0.001 --make-bed --out alz10
plink --bfile alz10 --genome --min 0.2 --out pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile alz10 --missing --out missing_report
echo "0    D180
0    C020
0    AK015
0    C132
0    C124
0    C077
0    C067
0    C170" > relatives_to_remove.tsv
plink --bfile alz10 --remove relatives_to_remove.tsv --allow-no-sex --make-bed --out alz11
plink --bfile alz11 --genome --min 0.2 --out 2pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' 2pihat_min0.2.genome > 2related_pairs.txt
plink2 --bfile alz11 --pca 10 --out alz_pca
```

```bash
plink --bfile alz11 --covar alz_pca.eigenvec --logistic --hide-covar --thread-num 8 --out simple_logistic
plink --bfile alz11 --covar alz_pca.eigenvec --logistic --dominant --hide-covar --out dominant_results
plink --bfile alz11 --covar alz_pca.eigenvec --logistic --recessive --hide-covar --out recessive_results
plink --bfile alz11 --allow-no-sex --assoc --out assoc_results
```

```bash
cat dominant_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat recessive_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat simple_logistic.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat assoc_results.assoc | awk '$9 != "NA"' | sort -gk 9,9 | head 
```
No association so far! MAX = 10e-7, 10-20 10e-5
