# GWAS-MSL-Dog
Script repository for "Role of canid hosts for the genetic remodeling of Leishmania infantum strain" publication.
Authors: Luís Fábio Batista and João Luís Reis-Cunha.

Log - Editing datasets to run PCA of dogs based in MSL (DEL X Non-DEL)

 #To do a merge.map file

   ./plink --bfile merge --recode --tab --out merge --autosome-num 38 --allow-extra-chr

#Convert merge dataset in Plink format to VCF format

  ./plink --file merge --recode vcf --out mergepca --autosome-num 38 –allow-extra-chr

#Remove chromosome ‘0’ from the vcf file

vcftools --vcf mergepca.vcf --not-chr 0 --recode --recode-INFO-all --out mergepca_outchr0

#Genomic filtering
--maf: Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value. Allele frequency is defined as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site.
--mac: Include only sites with Minor Allele Count greater than or equal to the "--mac" value
--max-missing: Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
--missing-indv
Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
--remove-indv
Removed 3 individuals with imiss > 0.1
Obs.: A amostra BAU08 tem 34% de Missiness mas foi mantida após filtragem pq é um dos cães cujo parasito foi filtrado para MSL.

vcftools --vcf mergepca.vcf --not-chr 0 --max-missing 0.95 --mac 1 --maf 0.01 --remove-indv LBD_PRU02 --remove-indv LBD_PRU04 --remove-indv LBD_PRU07 --recode --recode-INFO-all --out mergepca_filtered.recode
#It returned 234 from 237 individuals and 149456 SNPs from 173,662 of the original set.

#make plink files
vcftools --plink --vcf mergepca_filtered.recode.vcf --out mergepca_filtered.recode

#make a bed file
./plink --file mergepca_filtered.recode --make-bed --dog --out mergepca_filtered.recode

#To do Exact test for Hard-Weinberg Equilibrium and do the H-W list:
./plink --bfile mergepca_filtered.recode --hardy --dog 
#To sort based in Exact test for Hard-Weinberg Equilibrium p-value:
 sort -k9,9n plink.hwe > plink_sorted.hwe
	more plink_sorted.hwe
  	tail plink_sorted.hwe

 #To do a list of SNP with p-value < 10e-5 in Exact test for Hard-Weinberg Equilibrium:
awk '{if($9 <= 10e-5) print $0}' plink_sorted.hwe > plink_sorted_filt.hwe
 more plink_sorted_filt.hwe
 wc plink_sorted_filt.hwe
Which comprises 31586  SNPs

#To remove SNPs out of Hard-Weimberg Equilibrium (< 10e-5):
./plink --bfile mergepca_filtered.recode --exclude plink_sorted_filt.hwe --make-bed --out   mergepca_filtered.recode_hwe.filt --dog
#Which returned 117870 in mergepca_filtered.recode_hwe.filt from 149456  SNPs in   mergepca_filtered.recode


### GWAS for MSL:

#lets create a kinship matrix:
./ldak5.2.linux --calc-kins-direct mergepca_filtered.recode_hwe.filt.kinship --bfile mergepca_filtered.recode_hwe.filt --ignore-weights YES --power -0.25


#19/09/2023 - Luís Fábio Batista

#Testing GWAS Linear model for MSL using PHENO_DEF.txt as phenotype file. After change all of traits to number code. Legend for all traits are in file "Códigos dos #fenótipos_BEPE" no path #"D:\Documentos\DOCUMENTOS EM 18-9-2017\Documentos\POS DOC\BEPE\Analise_Projeto\Códigos dos fenótipos_BEPE.txt" 

#The path for this GWAS for MSL using bfile (genotypes) filtered for call rate 95% and Exact Hard-Weimberg Test p-value > 1X10e-5 is "D:\Documentos\DOCUMENTOS EM 18-9-2017\Documentos\POS DOC\BEPE\Analise_Projeto\GWAS_MSL_H-W_test\"

#To teste GWAS Linear model for MSL I first used:
./ldak5.2.linux --linear GWASbySNP_MSL --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41
#Linear regression results: GWASbySNP_MSL.assoc

#23/09/2023 - Luís Fábio Batista

#To test GWAS Linear model for MSL with ORIGIN as covariate
cut -f 1,2,37 PHENO_DEF.txt > origincovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_origincovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar origincovar.txt

#To test GWAS Linear model for MSL with SEX as covariate
cut -f 1,2,38 PHENO_DEF.txt > sexcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_sexcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar sexcovar.txt

#To test GWAS Linear model for MSL with AGE as covariate
cut -f 1,2,39 PHENO_DEF.txt > agecovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_agecovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar agecovar.txt

#To test GWAS Linear model for MSL with REPELLENT COLLAr as covariate
cut -f 1,2,40 PHENO_DEF.txt > collarcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_collarcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar collarcovar.txt

#To test GWAS Linear model for MSL with VACCINE as covariate
cut -f 1,2,41 PHENO_DEF.txt > vaccinecovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_vaccinecovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar vaccinecovar.txt

#To test GWAS Linear model for MSL with TREATMENT as covariate
cut -f 1,2,42 PHENO_DEF.txt > treatmentcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_treatmentcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar treatmentcovar.txt

#To test GWAS Linear model for MSL with SEX, AGE, ORIGIN, REPELLENT COLLAR, VACCINE and TREATMENT as covariates
#cut -f 1,2,37-42 PHENO_DEF.txt > allcovars.txt
./ldak5.2.linux --linear GWASbySNP_MSL_allcovars --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar allcovars.txt
#Linear regression results: GWASbySNP_MSL_allcovars.assoc

#To test GWAS Linear model for MSL with CLINIC as covariate
cut -f 1,2,44 PHENO_DEF.txt > cliniccovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_cliniccovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar cliniccovar.txt

#To test GWAS Linear model for MSL with STAGING as covariate
cut -f 1,2,45 PHENO_DEF.txt > stagingcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_stagingcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar stagingcovar.txt

#To test GWAS Linear model for MSL with Parasite load in limph node as covariate
cut -f 1,2,46 PHENO_DEF.txt > PL.LFNDcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_PL.LFNDcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar LOG.PL.LFNDcovar.txt --max-iter
#Error, eigen decomp failed; please tell Doug (info 2, length 3)

#To test GWAS Linear model for MSL with Log10 of Parasite load in limph node as covariate
cut -f 1,2,47 PHENO_DEF.txt > LOG.PL.LFNDcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_LOG.PL.LFNDcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar LOG.PL.LFNDcovar.txt --max-iter 500 --tolerance 0.01

#To test GWAS Linear model for MSL with IgG as covariate
cut -f 1,2,48 PHENO_DEF.txt > IgGcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_IgGcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgGcovar.txt


#24/09/2023 - Luís Fábio Batista

#To test GWAS Linear model for MSL with IgA as covariate
cut -f 1,2,49 PHENO_DEF.txt > IgAcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_IgAcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgAcovar.txt

#To test GWAS Linear model for MSL with IgM as covariate
cut -f 1,2,50 PHENO_DEF.txt > IgMcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_IgMcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgMcovar.txt

#To test GWAS Linear model for MSL with IgE as covariate
cut -f 1,2,51 PHENO_DEF.txt > IgEcovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_IgEcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgEcovar.txt

#To test GWAS Linear model for MSL with IgG anti-SALIVA as covariate
cut -f 1,2,52 PHENO_DEF.txt > IgGsalivacovar.txt
./ldak5.2.linux --linear GWASbySNP_MSL_IgGsalivacovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgGsalivacovar.txt
#27/09/2023 - Luís Fábio Batista

#To test GWAS Linear model for MSL with SEX, AGE, ORIGIN, REPELLENT COLLAR, and VACCINE as covariates
cut -f 1,2,4-8 PHENO_DEF.txt > allcovars2.txt
./ldak5.2.linux --linear GWASbySNP_MSL_allcovars2 --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar allcovars2.txt
#Linear regression results: GWASbySNP_MSL_allcovars2.assoc

28/9/2023 - Luís Fábio Batista

## Permutation

#To run 100 permutations to GWAS MSL with no covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 7.8824e-09


#To run 100 permutations to GWAS MSL using ORIGIN as covariate:
#cut -f 1,2,4 FENOTIPOS_MERGE_FINAL2.txt > origincovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_origincovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar origincovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_origincovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_origincovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_origincovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 2.5156e-09


#To run 100 permutations to GWAS MSL using SEX as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_sexcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar sexcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_sexcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_sexcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_sexcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.1110e-09


#To run 100 permutations to GWAS MSL using AGE as covariate:
#cut -f 1,2,6 FENOTIPOS_MERGE_FINAL2.txt > agecovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_agecovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar agecovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_agecovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_agecovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_agecovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 8.6074e-09


#To run 100 permutations to GWAS MSL using REPELLENT COLLAR as covariate:
#cut -f 1,2,7 FENOTIPOS_MERGE_FINAL2.txt > collarcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_collarcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar collarcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_collarcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_collarcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_collarcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.1520e-08


#To run 100 permutations to GWAS MSL using VACCINE as covariate:
#cut -f 1,2,8 FENOTIPOS_MERGE_FINAL2.txt > vaccinecovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_vaccinecovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar vaccinecovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_vaccinecovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_vaccinecovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_vaccinecovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 2.6067e-09


#To run 100 permutations to GWAS MSL using TREATMENT as covariate:
#cut -f 1,2,9 FENOTIPOS_MERGE_FINAL2.txt > treatmentcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_treatmentcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar treatmentcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_treatmentcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_treatmentcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_treatmentcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 5.0998e-09


#To run 100 permutations to GWAS MSL using CLINIC as covariate:
#cut -f 1,2,11 FENOTIPOS_MERGE_FINAL2.txt > cliniccovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_cliniccovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar cliniccovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_cliniccovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_cliniccovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_cliniccovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 2.4330e-08


#To run 100 permutations to GWAS MSL using STAGING as covariate:
#cut -f 1,2,12 FENOTIPOS_MERGE_FINAL2.txt > stagingcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_stagingcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar stagingcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_stagingcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_stagingcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_stagingcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 1.3472e-09


#To run 100 permutations to GWAS MSL using PL.LFND (LOG) as covariate:
#cut -f 1,2,15 FENOTIPOS_MERGE_FINAL2.txt > pl.lfnd.covar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_PL.LFNDcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar LOG.PL.LFNDcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_PL.LFNDcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_PL.LFNDcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_PL.LFNDcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 3.8152e-09


To run 100 permutations to GWAS MSL using Anti-Leish IgG as covariate:
#cut -f 1,2,18 FENOTIPOS_MERGE_FINAL2.txt > IgGcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgGcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgGcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgGcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgGcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgGcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 5.3635e-09


To run 100 permutations to GWAS MSL using Anti-Leish IgA as covariate:
#cut -f 1,2,19 FENOTIPOS_MERGE_FINAL2.txt > IgAcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgAcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgAcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgAcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgAcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgAcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 1.0545e-08


To run 100 permutations to GWAS MSL using Anti-Leish IgM as covariate:
#cut -f 1,2,20 FENOTIPOS_MERGE_FINAL2.txt > IgMcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgMcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgMcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgMcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgMcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgMcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 2.9937e-09


To run 100 permutations to GWAS MSL using Anti-Leish IgE as covariate:
#cut -f 1,2,21 FENOTIPOS_MERGE_FINAL2.txt > IgEcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgEcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgEcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgEcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgEcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgEcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 2.6580e-08


To run 100 permutations to GWAS MSL using IgGSALIVA as covariate:
#cut -f 1,2,22 FENOTIPOS_MERGE_FINAL2.txt > IgGSALIVAcovar.txt
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgGSALIVAcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgGsalivacovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgGSALIVAcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgGSALIVAcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgGSALIVAcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 2.6703e-09


#To run 100 permutations to GWAS MSL using allcovars as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_allcovars_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar allcovars.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_allcovars_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_allcovars_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_allcovars_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.3236e-09



#To run 100 permutations to GWAS MSL using SEX as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_sexcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno pheno PHENO_DEF.txt --mpheno 41 --covar sexcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_sexcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_sexcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_sexcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 1.5442e-09


#To run 100 permutations to GWAS MSL using IgM as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgMcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno pheno PHENO_DEF.txt --mpheno 41 --covar IgMcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgMcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgMcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgMcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.3631e-09

#To run 100 permutations to GWAS MSL using SEX as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_sexcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar sexcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_sexcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_sexcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_sexcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 1.5442e-09


#To run 100 permutations to GWAS MSL using IgM as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_MSL_IgMcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgMcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_MSL_IgMcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_MSL_IgMcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_MSL_IgMcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.3631e-09

#### Heritability_Analysis_Luís Fábio da Silva Batista 25-42-2023

#Heritability_Analysis by RMEL model, using MSL as phenotype, without cvariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT --pheno PHENO_DEF.txt --mpheno 41 --grm mergepca_filtered.recode_hwe.filt.kinship 
#Result: Her_All -0.180701 0.404583 in file 'Heritability_MSL_H-W_FILT.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and ORIGIN as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_origincovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar origincovar.txt  
#Result: Her_All -0.169170 0.421952 e Covar_Heritability 0.0048 no arquivo 'Heritability_MSL_H-W_FILT_origincovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and SEX as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_sexcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar sexcovar.txt  
#Result: Her_All 0.030998 0.564111 e Covar_Heritability 0.0276 no arquivo 'Heritability_MSL_H-W_FILT_sexcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and AGE as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_agecovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar agecovar.txt  
#Result: Her_All -0.202214 0.386031 and Covar_Heritability 0.0210 at the file 'Heritability_MSL_H-W_FILT_agecovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and REPELLENT COLLAR as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_collarcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar collarcovar.txt  
#Result: Her_All 0.006466 0.516458 56720 and Covar_Heritability 0.0268 at the file 'Heritability_MSL_H-W_FILT_collarcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and VACCINE as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_vaccinecovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar vaccinecovar.txt  
#Result: Her_All -0.263675 0.333654 56720.42 e Covar_Heritability 0.0465 at the file 'Heritability_MSL_H-W_FILT_vaccinecovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and TREATMENT as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_treatmentcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar treatmentcovar.txt  
#Result: Her_All -0.219410 0.406471 and Covar_Heritability 0.0929 at the file 'Heritability_MSL_H-W_FILT_treatmentcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and CLINIC as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_cliniccovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar cliniccovar.txt  
#Result: Her_All -0.474932 0.209258 and Covar_Heritability 0.0503 at the file 'Heritability_MSL_H-W_FILT_cliniccovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and CLINICAL STAGING as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_stagingcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar stagingcovar.txt  
#Result: Her_All -0.480847 0.149248 Covar_Heritability 0.0793 at the file 'Heritability_MSL_H-W_FILT_stagingcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and Log of Parasite load in limph node as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_LOG.PL.LFNDcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar LOG.PL.LFNDcovar.txt  
#Result: Her_All -0.204834 0.405873 and Covar_Heritability 0.0043 at the file 'Heritability_MSL_H-W_FILT_LOG.PL.LFNDcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and IgG as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_IgGcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgGcovar.txt  
#Result: Her_All -0.329654 0.325828 and Covar_Heritability 0.0076 at the file 'Heritability_MSL_H-W_FILT_IgGcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and IgA as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_IgAcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgAcovar.txt  
#Result: Her_All 0.120521 0.690596 and Covar_Heritability 0.1243 at the file 'Heritability_MSL_H-W_FILT_IgAcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and IgM as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_IgMcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgMcovar.txt  
#Result: Her_All -0.194706 0.433181 and Covar_Heritability 0.0420 at the file 'Heritability_MSL_H-W_FILT_IgMcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and IgE as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_IgEcovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgEcovar.txt  
#Result: Her_All -0.182199 0.406200 and Covar_Heritability 0.0033 at the file 'Heritability_MSL_H-W_FILT_IgEcovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and IgG anti-SALIVA as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_IgGsalivacovar --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar IgGsalivacovar.txt  
#Result: Her_All -0.203886 0.392497 and Covar_Heritability 0.0028 at the file 'Heritability_MSL_H-W_FILT_IgGsalivacovar.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and ORIGIN, AGE, SEX, REPELLENT COLLAR, VACCINE, TREATMENT as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_allcovars --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar allcovars.txt  
#Result: Her_All -0.799148 0.125574 and Covar_Heritability 0.1226 at the file 'Heritability_MSL_H-W_FILT_allcovars.reml'

#Heritability_Analysis by RMEL model, using MSL as phenotype, and ORIGIN, AGE, SEX, REPELLENT COLLAR, VACCINE as covariate:
./ldak5.2.linux --reml Heritability_MSL_H-W_FILT_allcovars2 --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 41 --covar allcovars2.txt  
#Result: Her_All 0.081079 0.609885 and Covar_Heritability 0.0926 at the file 'Heritability_MSL_H-W_FILT_allcovars2.reml'

#### Running GWAS Single-Predictor for anti-_L.infantum_ traits including MSL as covariate

#To test GWAS Linear model for IgM ant-Leish without covariates
./ldak5.2.linux --linear GWASbySNP_IgM --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48
#Linear regression results: GWASbySNP_IgM.assoc

#To test GWAS Linear model for IgM ant-Leish with MSL as covariates
./ldak5.2.linux --linear GWASbySNP_IgM_MSLcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar MSLcovar.txt
#Linear regression results: GWASbySNP_IgM_MSLcovar.assoc

#To test GWAS Linear model for IgM ant-Leish with TREATMENT as covariates
./ldak5.2.linux --linear GWASbySNP_IgM_treatmentcovar --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar treatmentcovar.txt
#Linear regression results: GWASbySNP_IgM_treatmentcovar.assoc

#To test GWAS Linear model for IgM ant-Leish CLINICAL STAGING with SEX, AGE, ORIGIN, REPELLENT COLLAR, VACCINE, and TREATMENT as covariates
./ldak5.2.linux --linear GWASbySNP_IgM_allcovars --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar allcovars.txt
#Linear regression results: GWASbySNP_IgM_allcovars.assoc

#To test GWAS Linear model for IgM ant-Leish with SEX, AGE, ORIGIN, REPELLENT COLLAR, VACCINE, TREATMENT and MSL as covariates
./ldak5.2.linux --linear GWASbySNP_IgM_allcovarsMSL --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar allcovarsMSL.txt
#Linear regression results: GWASbySNP_IgM_allcovarsMSL.assoc

## Permutation

#To run 100 permutations to GWAS anti-Leish IgM with no covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_IgM_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_IgM_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_IgM_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_IgM_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 1.5335e-07


#To run 100 permutations to GWAS anti-Leish IgM using MSLcovar as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_IgM_MSLcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar MSLcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_IgM_MSLcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_IgM_MSLcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_IgM_MSLcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.0047e-07


#To run 100 permutations to GWAS anti-Leish IgM using TREATMENT as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_IgM_treatmentcovar_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar treatmentcovar.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_IgM_treatmentcovar_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_IgM_treatmentcovar_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_IgM_treatmentcovar_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 6.2016e-07


#To run 100 permutations to GWAS anti-Leish IgM using allcovars as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_IgM_allcovars_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar allcovars.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_IgM_allcovars_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_IgM_allcovars_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_IgM_allcovars_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 5.3477e-07


#To run 100 permutations to GWAS anti-Leish IgM using allcovarsMSL as covariate:
for ((i=1;i<=100;i++)); do ./ldak5.2.linux --linear simulation_GWASbySNP_IgM_allcovarsMSL_permutation.${i} --bfile mergepca_filtered.recode_hwe.filt --grm mergepca_filtered.recode_hwe.filt.kinship --pheno PHENO_DEF.txt --mpheno 48 --covar allcovarsMSL.txt --max-iter 500 --tolerance 0.01 --permute YES; done

#Obtaining the cutoff:
#To get the 5% percentile of the lowest Wald_P values for the 100 replicates:
for i in $(ls simulation_GWASbySNP_IgM_allcovarsMSL_permutation.*.assoc); do sort -g -k7,7 $i | awk '$7!="NA"{print $7}' | head -n 2 | tail -n 1 >> Z_GWASbySNP_IgM_allcovarsMSL_assoc_all_lowest; done

#As we have 100 replicates, the 5% percentile will be the position 5
sort -g Z_GWASbySNP_IgM_allcovarsMSL_assoc_all_lowest | head -n 5 | tail -n1
#Which returne the cutoff value of 4.6089e-07


### Linear Regression Model to test association between MSL genotype frequence and kernel Density of canine visceral leishmaniasiss and Kernel concentration of human VL cases in Municipality of Bauru - São Paulo.

setwd("/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/MANUSCRIPT_3/ARQUIVOS_REVISÃO_29-9-2025")
> pheno<-read.table(file = "/media/fmusp/TOSHIBA EXT/Documentos/DOCUMENTOS EM 18-9-2017/Documentos/POS DOC/BEPE/MANUSCRIPT_3/ARQUIVOS_REVISÃO_29-9-2025/MSLxKERNEL_ONLY_BAURU_OUT_MIX_tab_Manuscript_MSL-Dog_06-10-2025.txt", header=T, sep="\t", stringsAsFactors = F)

# Kernel Density of canine VL cases:
Call:
lm(formula = RISK_CLASS ~ SEX + AGE + AREA + MSL, data = pheno)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.2819 -0.4121 -0.1361  0.4968  1.4228 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)   
(Intercept)   7.0781     2.7850   2.542  0.01943 * 
SEX          -0.5943     0.4122  -1.442  0.16479   
AGE          -1.3508     0.8766  -1.541  0.13898   
AREAB         0.9615     0.4866   1.976  0.06212 . 
AREAC         1.1458     0.4116   2.784  0.01146 * 
AREAD         1.1075     0.3830   2.892  0.00902 **
MSL          -0.8896     0.4031  -2.207  0.03915 * 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.7343 on 20 degrees of freedom
  (16 observations deleted due to missingness)
Multiple R-squared:  0.4567,	Adjusted R-squared:  0.2938 
F-statistic: 2.802 on 6 and 20 DF,  p-value: 0.03816

# Kernel Concentration of human VL cases:
Call:
lm(formula = HUMAN_RISK_CLASS ~ SEX + AGE + AREA + MSL + REPELLENT_COLLAR, 
    data = pheno)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.10679 -0.49671 -0.08481  0.35605  1.91519 

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)  
(Intercept)       1.52027    3.89873   0.390   0.7009  
SEX               0.38837    0.55807   0.696   0.4949  
AGE              -0.08823    1.24415  -0.071   0.9442  
AREAB             1.38521    0.65924   2.101   0.0492 *
AREAC            -0.29360    0.55799  -0.526   0.6049  
AREAD            -0.56141    0.56232  -0.998   0.3306  
MSL               0.82922    0.54658   1.517   0.1457  
REPELLENT_COLLAR  0.46284    0.43488   1.064   0.3005  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.9942 on 19 degrees of freedom
  (16 observations deleted due to missingness)
Multiple R-squared:  0.4949,	Adjusted R-squared:  0.3088 
F-statistic:  2.66 on 7 and 19 DF,  p-value: 0.04254

