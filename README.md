<h1 align="center">
    <b>🧬 POPULATION GENETICS - PIPELINE DOCUMENTATION </b> <br>
    (LD, HWE, HWD, MAF, EM)
</h1>

<h3 align="center">
    Author: Walter Odur <br>
    Affiliation: ACE-Uganda
</h3>

---

## 📌 By the end of this tutorial, you will have gained insights into;
1.  Measuring pairwise **`Linkage Disequilibrium (LD)`** for a group of SNPs
2.  Testing for **`Hardy-Weinberg Equilibrium HWE`**
3.  Testing for **`Hardy-Weinberg Disequilibrium HWD`**
4.  Identifying **`population substructure`**
5.  Adjusting for multiple testing; **`Bonferroni Correction`**
6.  Identifying **`Minor Alleles`** and their frequencies **`(MAF)`**
7.  **`Haplotype frequency`** estimation; **`Expectation Maximization (EM)`** approach
8.  **`Testing hypotheses`** about haplotype frequencies within the EM framework  
Required packages;

-  [genetics](https://rdrr.io/cran/genetics/)
-  [MASS](https://cran.r-project.org/web/packages/MASS/index.html)
-  [combinat](https://cran.r-project.org/web/packages/combinat/index.html)
-  [gdata](https://cran.r-project.org/web/packages/gdata/index.html)
-  [gtools](https://cran.r-project.org/web/packages/gtools/index.html)
-  [mvtnorm](https://cran.r-project.org/web/packages/mvtnorm/index.html)
-  [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
-  [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap)
-  [haplo.stats](https://cran.r-project.org/web/packages/haplo.stats/index.html)

## **1: Measuring Pairwise LD for a Group of SNPs**  

We are interested in calculating pairwise LD for a group of SNPs within four of the 17 genes 
studied in the FAMuSS study. The FAMuSS study found various genetic loci associated with 
muscle size and strength at baseline and in response to resistance training; find the [article here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3885233/). 
The four genes that will be used are: 
-  ANKRD6 (Ankyrin repeat domain 6), ACTN3 (Actinin, alpha 3), LEPR (Leptin receptor) and RETN (Resistin).
  
📌 How many SNPs in total if you add up those in all our four genes: ANKRD6, ACTN3, LEPR and RETN?

-  Create an upper triangular matrix of pairwise LD measures for the SNPs. Compare the level of linkage disequilibrium between genes and that across genes.
-  The result can also be illustrated using a heatmap. Create a heatmap from the matrix of pairwise LD measures for the SNPs.
  
📌 What do you observe? 
In the above exercise, we considered measures of pairwise LD between two alleles. More generally, 
interest may lie in determining whether a group of alleles are in Linkage Disequilibrium. One 
intuitively appealing measure of LD across a region comprised of multiple SNPs is simply the 
average of all pairwise measures D'. For each of the four genes we are exploring, determine the 
average LD. Remember to specify that missing values should be removed by including 
na.rm=TRUE. Which genes are regions of high LD? (Considering high to be average LD >= 0.7). 
Hints: The str_detect function in package stringr can be used to search SNPs in the above three 
genes (ankrd6, lepr, resistin, actn3). You can use regular expressions to search for SNPs beginning 
with those letters. You then proceed by subsetting the FMS dataset for only those SNPs. You may 
then create a list object with each element of the list containing a genotype object for each SNP. 
The lapply function can be used to apply the genotype function to a list of SNPs. The result 
can be 

### **1.1 How many SNPs in total if you add up those in all our four genes: ANKRD6, ACTN3, LEPR, and RETN?**

```r
# Load necessary packages
library(genetics)
library(MASS)
library(combinat)
library(gdata)
library(gtools)
library(mvtnorm)
library(stringr)
library(pheatmap)
library(haplo.stats)

# Read in data
FMS_data <- read.table("FMS_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Select genes
selected_genes <- c("ankrd6", "actn3", "lepr", "resistin")
selected_snps <- FMS_data[, str_detect(names(FMS_data), paste(selected_genes, collapse = "|"))]

num_snps <- ncol(selected_snps)
cat("Total SNPs in the four genes:", num_snps, "\n")




# Convert character SNPs into correct genotype format
selected_snps[] <- lapply(selected_snps, function(x) {
  x <- as.character(x)  
  x <- gsub("(\\w)(\\w)", "\\1/\\2", x)  
  factor(x, levels = unique(na.omit(x)))  
})

snp_genotypes <- lapply(selected_snps, genotype)
snp_matrix <- as.data.frame(snp_genotypes)
ld_matrix <- LD(snp_matrix)$`D'`
upper_ld_matrix <- ld_matrix
upper_ld_matrix[lower.tri(upper_ld_matrix, diag = TRUE)] <- NA  















