data/GTEX_phenotype.gz is the phenotype corresponding file of GTEx dataset (normal colon tissue).
data/gtex_RSEM_gene_fpkm.gz is the expression data file of GTEx dataset (expression data of normal colon tissue).
1. run GTEx.R section 1, extract GTEx data
2. run GTExSymbol.pl, convert gene ID to symbol (gene name)
3. run GTEx.R section 2, delete genes with low expression of GTEx dataset.

data/gdc_download_20220321_023306.495140.tar.gz is the expression data file of TCGA dataset (expression data of colorectal tumor tissue).
data/gdc_download_20220324_141628.148216.tar.gz is the clinical data file of TCGA dataset (colorectal tumor tissue).
data/gene.txt is a list of immune genes.
1. run TCGAmoveFiles.pl and TCGAmerge.pl,  extract TCGA expression data
2. run TCGAgetMrna.pl, extract mRNA and convert gene ID to symbol of TCGA dataset
3. run TCGAgetClinical.pl, extract TCGA clinical data
4. run TCGA.R section 1, delete genes with low expression of TCGA dataset.
5. run TCGA.R section 2, delete TCGA normal sample
6. run TCGAmergeExpTime.pl, Intersect TCGA expression data and clinical data
7. run TCGA.R section 3, Merge TCGA and GTEx
8. run TCGA.R section 4, differentially expressed genes between TCGA and GTEx
9. run TCGA.R section 5, differentially expressed immune genes between TCGA and GTEx
10. run TCGA.R section 6, get differentially immune genes expression in TCGA
11. run TCGAmergeExpTime.pl, merge TCGA expression data and clinical data
12. run TCGA.R section 7, Divided into training group and verification group

construct model
1. run model.R section 1, LASSO analysis of trainging dataset, get 14 prognosis-related immune genes.
2. run TCGAgetLassoExp.pl, get Multivariate Cox analysis input file
3. run model.R section 2, Multivariate cox analysis, get 9 prognosis-related immune genes from 14 candidate genes to calculate risk score
4. run model.R section 3, get optimal cut-off value by surv_cutpoint function, and survival analysis of training dataset
5. run model.R section 4, plot heatmap and survstat map
6. run model.R section 5, ROC curve of training dataset
7. run TCGAprepareIndep.pl, prepare Independent prognostic analysis input file
8. run model.R section 6, Multivariate independent prognostic analysis of training dataset

model validation
1. run ModelValidation.R section 1, calculate risk score of testing dataset
2. run ModelValidation.R section 2, get optimal cut-off value by surv_cutpoint function, and survival analysis of testing dataset
3. run ModelValidation.R section 3, ROC curve of testing dataset
4. run ModelValidation.R section 4, calculate risk score of entire TCGA dataset
5. run ModelValidation.R section 5, get optimal cut-off value by surv_cutpoint function, and survival analysis of entire TCGA dataset
6. run ModelValidation.R section 6, ROC curve of entire TCGA dataset