# Clear the workspace
rm(list=ls())

setwd("/blue/mresende/share/viannam/GS")

# Load required packages
library(BGLR)
library(AGHmatrix)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

## Phenotyping information
## Pheno file should have geno ID column (same order as in the Gmatrix) + trait columns

# Read phenotypic data
BLUEs_trait <- read.csv("/blue/mresende/share/viannam/GS/newBLUES_bySeason_CAP.csv", sep = ",")

# Display the first few rows of the phenotypic data
head(BLUEs_trait)

## K-mers GS

# Assure that the genotypes are in the same order on both files
kmers_names <- read.table("/blue/mresende/share/viannam/Kmers/kmersGWAS/kmers_table.names")

# Display the first few rows of kmers_names
head(kmers_names)

# Reorder rows of BLUEs_trait based on genotype IDs
BLUEs <- BLUEs_trait[match(kmers_names$V1, BLUEs_trait[,1]), ]

# Display the first few rows of BLUEs
head(BLUEs)

## K-mers relationship matrix

# Read the K-mers relationship matrix
kmers_kinship <- read.table("/blue/mresende/share/viannam/Kmers/kmersGWAS/kmers_table.kinship", header = FALSE)
G_mat <- as.matrix(kmers_kinship)
colnames(G_mat) <-  kmers_names$V1
rownames(G_mat) <- kmers_names$V1


# GBLUP K-mers Model

# Extract phenotypic data for traits
y <- (BLUEs[, 3:22])

# Define a function to run GBLUP with cross-validation
run_GBLUP_CV <- function(data, folds, kinship, nIter, burnIn, nReps) {
  predM1_all <- data.frame()
  traitNames = colnames(data)

  for (trait in 1:ncol(data)) {
    predM1 <- data.frame()
    Y <- data[, trait]

    for (Rep in 1:nReps) {
      nFolds= sample(1:folds, size = length(Y), replace = TRUE)

      for (i in 1:max(nFolds)) {
        tst <- which( nFolds == i)
        yNA <- Y
        yNA[tst] <- NA

        # Fit GBLUP model
        fm <- BGLR(y = as.matrix(yNA),
                   ETA = list(list(K = as.matrix(kinship), model = 'RKHS')),
                   nIter = nIter,
                   burnIn = burnIn,
                   verbose = FALSE)

        # Predicted values
        PC <- cor(Y[tst], fm$yHat[tst], use = "pairwise.complete.obs")
        predM1 <- rbind(predM1, data.frame(Trait = traitNames[trait],
                                           k_Fold = i,
                                           Rep = Rep,
                                           Acc = PC))
      }
    }
    predM1_all <- rbind(predM1, predM1_all)  # Store results for each trait in a list element
  }
  return(predM1_all)
}


# Run GBLUP with cross-validation
set.seed(1234)
predM1_all <- run_GBLUP_CV(data = y,  # Extract phenotypic data
                           folds = 5,
                           kinship = G_mat,
                           nIter = 10000,
                           burnIn = 1000,
                           nReps = 10)

# Save model output for each trait

save(predM1_all, file = "/blue/mresende/share/viannam/GS/Accuracy_GBLUP_kmers_Multi-Year.RData")
