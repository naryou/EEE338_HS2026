# GWAS/GS exercise using GAPIT
## temporary: new URL for genome browser (slide 23): plants.ensembl.org/Oryza_sativa also try (this)[https://rice.uga.edu/jb2/?session=local-FCWky9JwdGzmRw33L332K] 

## 0. start RStudio 

on your computer or at [FGCZ](https://fgcz-genomics.uzh.ch)

## 1. Setup and load packages

Clean up a workplace first.
```
# Detach all non-base packages
pkgs <- paste("package:", names(sessionInfo()$otherPkgs), sep = "")
if (length(pkgs) > 0) {
  invisible(lapply(pkgs, detach, character.only = TRUE, unload = TRUE))
}

rm(list=ls())
```

Load GAPIT source code and its dependency. Wait ca. 15 min. to install everything.  
Some packages are not installed but they are negligible.  
```
install.packages("devtools")
# BiocManager::install("snpStats") # It is installed on the system. You don't need to install it.
devtools::install_github("SFUStatgen/LDheatmap")
devtools::install_github("jiabowang/GAPIT3@078fe28",force=TRUE) # install ver Aug. 2022
library(GAPIT3) # load GAPIT3 package
```

## 2. Load phenotypes and genotypes

Download phenotype data
```
url <- "http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"
p <- read.table(url, sep="\t", header=TRUE)
nrow(p) # No. of plants
head(p)
```

Read genotype data and marker information
```
g <- read.table("RiceDiversity_44K_Genotypes_PLINK_imputed.txt.gz",
                header=TRUE, sep="\t")
gm <- read.table("RiceDiversity_44K_Genotypes_PLINK_info.txt.gz",
                 header=TRUE, sep="\t")
nrow(g[,-1]) # No. of plants
ncol(g) # No. of SNPs
head(gm) # marker info
```

## 3. Run GWAS

Compare the general linear model (GLM) and mixed linear model (MLM).  
Some warnings occur but the analysis still works. When they are finished, output files appear in the current directory.  
```
myGAPIT <- GAPIT(
  Y=p[,c("HybID","Seed.length")],
  GD=g,
  GM=gm,
  SNP.MAF=0.05, # cut-off minor alleles at 0.05
  Inter.Plot=TRUE,
  model=c("GLM", "MLM"),
  kinship.algorithm="VanRaden",
  Multiple_analysis=TRUE)
```

**Note: When you run GAPIT twice, the second run may not work. In such a case, log-out once and retry from data loading.**  

## 4. gBLUP

Calculate BLUP for the flowering time 2006 at Arkansas.  
Some warnings occur but the analysis still works.  
```
myGAPIT_BLUP <- GAPIT(
  Y=p[,c("HybID","Year06Flowering.time.at.Arkansas")],
  GD=g,
  GM=gm,
  SNP.MAF=0.05,
  model="gBLUP",
  kinship.algorithm="VanRaden",
  file.output=FALSE)
```

Load results of genomic prediction  
```
pred <- myGAPIT_BLUP$Pred
head(pred)
```

Align predicted and observed traits following the taxa name  
```
pred <- pred[order(pred$Taxa),]
y <- p[order(p$HybID),]
```

Pearson's correlation between predicted and observed flowering  
```
# calculate Pearson's correlation between predicted and observed flowering
cor.test(pred$Prediction, y$Year06Flowering.time.at.Arkansas, method = 'pearson') 
```

Predicting flowering time 2007  
```
# perform a linear regression to estimate the slope and intercept
res <- lm(y$Year07Flowering.time.at.Arkansas~pred$Prediction)

# plot the results
plot(pred$Prediction,
     y$Year07Flowering.time.at.Arkansas,
     ylab="flowering 2007", xlab="predicted",
     main=paste("r =",round(sqrt(summary(res)$r.squared),2)))
abline(res)
```

## 5. Correlation between masked-then-predicted and truly observed values, and predicting flowering time of missing accessions
```
## the section below is added this year (2026) (using Claude):
## ---------------------------------------------------------------
## Two separate things:
##   1) k-fold cross-validation -> a real accuracy estimate, using only
##      accessions where the true value is known (so we can check
##      predictions against ground truth).
##   2) predictions for the accessions that have NO observed phenotype
##      at all -> no accuracy number is possible here, only the
##      predicted values themselves.
## ---------------------------------------------------------------

set.seed(123) # reproducible fold assignment

## ---- 1) k-fold cross-validation ----

k <- 5
known_idx <- which(!is.na(p$Flowering.time.at.Aberdeen)) # only accessions we can score
folds <- sample(rep(1:k, length.out = length(known_idx))) # random fold label per known accession

cv_results <- data.frame(
  HybID = character(),
  Observed = numeric(),
  Predicted = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:k) {
  # mask this fold's phenotypes as if they were unobserved
  p_cv <- p
  mask_idx <- known_idx[folds == i]
  p_cv$Flowering.time.at.Aberdeen[mask_idx] <- NA
  
  gapit_cv <- GAPIT(
    Y = p_cv[, c("HybID", "Flowering.time.at.Aberdeen")],
    GD = g,
    GM = gm,
    SNP.MAF = 0.05,
    model = "gBLUP",
    kinship.algorithm = "VanRaden",
    file.output = FALSE
  )
  
  pred_cv <- gapit_cv$Pred
  # explicit match
  pred_cv <- pred_cv[match(p$HybID[mask_idx], pred_cv$Taxa), ]
  
  cv_results <- rbind(cv_results, data.frame(
    HybID = p$HybID[mask_idx],
    Observed = p$Flowering.time.at.Aberdeen[mask_idx],
    Predicted = pred_cv$Prediction
  ))
}

# accuracy: correlation between masked-then-predicted and truly observed values
cor.test(cv_results$Observed, cv_results$Predicted, method = "pearson")

res_cv <- lm(Observed ~ Predicted, data = cv_results)
plot(cv_results$Predicted, cv_results$Observed,
     xlab = "Predicted (5-fold cross-validation)",
     ylab = "Observed Flowering.time.at.Aberdeen",
     main = paste("Cross-validation r =", round(sqrt(summary(res_cv)$r.squared), 2))
)
abline(res_cv)

## ---- 2) predictions for accessions with NO observed phenotype ----
## (these were never masked above - they're actually missing in the raw data)

NA06 <- is.na(p$Flowering.time.at.Aberdeen)

# `pred` here is myGAPIT_BLUP$Pred from the original full-data gBLUP run
missing_predictions <- data.frame(
  Taxa = p$HybID[NA06],
  Predicted_Flowering = pred$Prediction[match(p$HybID[NA06], pred$Taxa)]
)
head(missing_predictions)
nrow(missing_predictions)

#
```

