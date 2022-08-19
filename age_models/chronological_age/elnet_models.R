#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena


library("optparse")
#library("glmnet")
#library("tidyverse")
library("biglasso")
library("bigmemoryExt")
library("matrixStats")

# Note: be careful about excess storage here: /dev/shm!
# Also do lsof /dev/shm for hidden things
# To kill hidden processes: pkill -U ebernab3

# Get arguments
######################################################

option_list = list(
    make_option(c("--meth"), type="character", default=NULL, 
              help="Methylation prepped data file", metavar="character"),
    make_option(c("--methdir"), type="character", default=NULL, 
              help="Directory for methylation data file", metavar="character"),
    make_option(c("--pheno"), type="character", default=NULL, 
              help="Pheno file with fold info location", metavar="character"),
    make_option(c("--name"), type="character", default=NULL, 
              help="Run name", metavar="character"),
    make_option(c("--sex"), type="character", default=NULL, 
              help="Sex stratified option (T/F)", metavar="character"),
    make_option(c("--cpg"), type="character", default=NULL, 
              help="Scale option per CpG site (T/F)", metavar="character"),
    make_option(c("--sample"), type="character", default=NULL, 
              help="Scale option per sample (T/F)", metavar="character"),
    make_option(c("--subset"), type="character", default=NULL, 
              help="Subset total number of CpGs (T/F)", metavar="character"),
    make_option(c("--subfil"), type="character", default=NULL, 
              help="Location of file containing CpGs to subset to", metavar="character"),
    make_option(c("--squared"), type="character", default=NULL, 
              help="CpG squared option (T/F)", metavar="character"),
    make_option(c("--squaredsubset"), type="character", default=NULL, 
              help="CpG squared subset option (T/F)", metavar="character"),
    make_option(c("--squaredsubsetfil"), type="character", default=NULL, 
              help="Location of file containing CpGs^2 to subset to", metavar="character"),
     make_option(c("--squaredsubsetfil_cpg2"), type="character", default=NULL, 
              help="Location of file containing CpGs^2 to subset to (CpG^2 lm models)", metavar="character"),             
    make_option(c("--squaredcpg2"), type="character", default=NULL, 
              help="Obtain squared CpGs from CpG2 linear model", metavar="character"),
    make_option(c("--age2cpg2"), type="character", default=NULL, 
              help="Use subsets from from age^2 EWAS and from age ~ CpG^2 linear models", metavar="character"),
    make_option(c("--external"), type="character", default=NULL, 
              help="External data option (T/F)", metavar="character"),
    make_option(c("--ext_select"), type="character", default=NULL, 
              help="External select sets, separated by commas", metavar="character"),
    make_option(c("--random"), type="character", default=NULL, 
              help="Totally random option (T/F) - will ignore assigned folds and randomize data", metavar="character"),
    make_option(c("--lasso"), type="character", default=NULL, 
              help="Lasso instead of elnet option (T/F)", metavar="character"),
    make_option(c("--logage"), type="character", default=NULL, 
              help="Option to run models with log(age) instead of age (T/F)", metavar="character"),
    make_option(c("--out"), type="character", default=NULL, 
              help="Output directory", metavar="character")           
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


set.seed(1234) # Set seed to ensure fold variation minimised 
seed <- 1234


# Function for faster row scaling
# Creds: https://www.r-bloggers.com/2016/02/a-faster-scale-function/
######################################################

rowScale = function(x,
    center = TRUE,
    scale = TRUE,
    add_attr = TRUE,
    rows = NULL,
    cols = NULL) {
    if (!is.null(rows) && !is.null(cols)) {
        x <- x[rows, cols, drop = FALSE]
    } else if (!is.null(rows)) {
        x <- x[rows, , drop = FALSE]
    } else if (!is.null(cols)) {
        x <- x[, cols, drop = FALSE]
    }
    cm = rowMeans(x, na.rm = TRUE)
    if (scale) {
        csd = rowSds(x, center = cm)
    } else {
        # just divide by 1 if not
        csd = rep(1, length = length(cm))
    }
    if (!center) {
        # just subtract 0
        cm = rep(0, length = length(cm))
    }
    x = (x - cm) / csd
    if (add_attr) {
        if (center) {
            attr(x, "scaled:center") <- cm
        }
        if (scale) {
            attr(x, "scaled:scale") <- csd
        }
    }
    return(x)
}


# Altered cbindBM for big, big matrices
# OG stuff by cdeterman on GitHub!
######################################################

cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}


# x needs to be a list of big matrices, and they have to have the same rownames
cbindBM_list <- function(x, binding="right", 
                    z=NULL, type=NULL, separated=NULL,
                    backingfile=NULL, backingpath=NULL,
                    descriptorfile=NULL, binarydescriptor=FALSE,
                    shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
        cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
        cols_list <- append(cols_list, list(cols))
        total_cols <- total_cols + ncol(x[[i]])
  }    

  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
        z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
        z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)

    if (erase == TRUE) {
        cat("\nErasing chunk and liberating memory...\n\n")
        x[[i]] <- "Replacement"
        gc()
    }
  }
  return(z)
}


# For terminal testing
######################################################

#sex_s <- "F"
#cpg_s <- "F"
#ext <- "F"
#rand <- "F"
#pheno <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/random/gs_basic_folds_random.tsv", header = TRUE, row.names = 2)
#xd <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/"
#xf <- "methylation_training_noscale.rds"
#meth <- readRDS(paste0(xd, xf))
#cpg_s_subset <- "T"
#cpg_s_df <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/subsets/age2_ewas_900.tsv", row.names = 1, header = TRUE)
#scale_cpg <- "F"
#scale_sample <- "F"
#lasso <- "F"
#subset <- "T"
#subfil <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/subsets/age_ewas_20K.tsv"

# Log age?
######################################################

logage <- opt$logage


# Scale options?
######################################################

scale_cpg <- opt$cpg
scale_sample <- opt$sample


# Using LBC + GEO in training as well?
######################################################

ext <- opt$external
ext_select <- unlist(strsplit(opt$ext_select, ","))


# Sex stratified models?
######################################################

sex_s <- opt$sex


# Include CpG^2?
######################################################

cpg_s <- opt$squared # All CpGs in model will also be included as CpG^2
cpg_s_subset <- opt$squaredsubset # Only CpGs in subset file will be included as CpG^2 (and if not in linear CpGs, will be added as well)
cpg_s_cpg2 <- opt$squaredcpg2 # Taking CpGs from lm on age vs CpG^2 instead of age^2 EWAS
cpg_s_age2cpg2 <- opt$age2cpg2 # Subsets from both age^2 EWAS and from age ~ CpG^2 lm
if ((cpg_s_subset == "T") & (cpg_s_age2cpg2 == "F")) {
    cpg_s_df <- read.delim(opt$squaredsubsetfil, row.names = 1, header = TRUE)
}
if ((cpg_s_subset == "T") & (cpg_s_age2cpg2 == "T")) {
    cpg_s_df <- read.delim(opt$squaredsubsetfil, row.names = 1, header = TRUE)
    cpg_s_df_cpg2 <- read.delim(opt$squaredsubsetfil_cpg2, row.names = 1, header = TRUE)
}


# Totally random fold assignment?
######################################################

rand <- opt$random


# Lasso instead of elnet?
######################################################

lasso <- opt$lasso


# For file naming
######################################################

lasso_f <- ""
if (lasso == "T") {
    lasso_f <- "_lasso"
}

random_f <- ""
if (rand == "T") {
    random_f <- "_randomizedfolds"
}

squared_f <- ""
if (cpg_s == "T") {
    squared_f <- "_squared"
}

squared_subset_f <- ""
if (cpg_s_subset == "T") {
    squared_subset_f <- paste0("_squaredsubset", length(rownames(cpg_s_df)))
}

squared_cpg_f <- ""
if ((cpg_s_cpg2 == "T") & (cpg_s_age2cpg2 == "F")) {
    squared_subset_f <- paste0("_squaredsubset", length(rownames(cpg_s_df))/1000, "K")
    squared_cpg_f <- "_cpg2lm"
}

if ((cpg_s_cpg2 == "T") & (cpg_s_age2cpg2 == "T")) {
    squared_subset_f <- paste0("_squaredsubset", length(rownames(cpg_s_df)), "_squaredsubsetcpg2lm", length(rownames(cpg_s_df_cpg2))/1000, "K")
}


# Import methylation plus pheno data
######################################################

# Pheno
if (ext == "T") {
    pheno <- read.delim(opt$pheno, header = TRUE, row.names = 1)
    pheno <- pheno[pheno$cohort %in% c("W1", "W3", "W4", ext_select),]
} else {
    pheno <- read.delim(opt$pheno, header = TRUE, row.names = 2)
}

if (logage == "T") {
    pheno$age <- log(pheno$age)
}

if (sex_s == "T") {
    pheno_F <- pheno[pheno$sex == "F",]
    pheno_M <- pheno[pheno$sex == "M",]
}
cat("\nImported pheno + fold data!\n")
cat(paste0("Remaining individuals: ",length(rownames(pheno)),"\n"))

# Methylation
xd <- opt$methdir
xf <- opt$meth
# Import
meth <- readRDS(paste0(xd, xf))
meth <- meth[which(rownames(meth) %in% rownames(pheno)),]
pheno <- pheno[which(rownames(pheno) %in% rownames(meth)),] 
cat("\nImported methylation data!\n")

cat("\nData prep...\n")


# Filter to subset of CpGs if option has been selected
######################################################

if (opt$subset == "T") {
    cpgs <- rownames(read.delim(opt$subfil, header = TRUE, row.names = 1))
    if ((cpg_s_subset == "T") & (cpg_s_age2cpg2 == "F")) {
        cpgs <- unique(c(cpgs, rownames(cpg_s_df))) # Keep those in CpG^2 subset as well (either from age^2 EWAS or CpG^2 ~ age lm)
    }
    if ((cpg_s_subset == "T") & (cpg_s_age2cpg2 == "T")) {
        cpgs <- unique(c(cpgs, rownames(cpg_s_df), rownames(cpg_s_df_cpg2))) # Keep those in CpG^2 subset as well (for both age^2 EWAS and CpG^2 ~ age lms)
    }
    meth <- meth[,which(colnames(meth) %in% cpgs)]
    cat("\nSubset to CpGs of interest. Resulting dimensions: ")
    cat(dim(meth))
}


# Scale data (if selected)
######################################################

# Scale per sample -> This part is very memory hungry, be careful! Also slow!
if (scale_sample == "T") {
    meth <- rowScale(meth)
    cat("\n\nScaled per sample as per user preference.")    
}

cat("\n\nRAM clean up...\n\n")
gc()

# Scale per CpG, per fold
if ((scale_cpg == "T") & (sex_s == "F")) {
    for (fold in unique(pheno[["Fold"]])) {
        samples <- rownames(pheno[pheno[["Fold"]] == fold,])
        sub <- scale(meth[samples,])
        rownames(sub) <- samples
        meth[samples,] <- sub
        rm(sub)
    }
    cat("\nScaled per CpG, per fold, as per user preference.\n")
} else if ((scale_cpg == "T") & (sex_s == "T")) {
    for (fold in unique(pheno[["Fold"]])) {
        for (s in c("F", "M")) {
            samples <- rownames(pheno[(pheno$Fold == fold) & (pheno$sex == s),])
            sub <- scale(meth[samples,])
            rownames(sub) <- samples
            meth[samples,] <- sub
            rm(sub)
        }
    }
    cat("\nScaled per CpG, per fold and per sex, as per user preference.\n")
}

cat("\nRAM clean up...\n\n")
gc()


# Check order and divide if sex-stratifying
######################################################

# Check rownames are the same and in the same order
if (sex_s == "T") {
    meth_F <- meth[rownames(meth) %in% rownames(pheno_F), ]
    meth_M <- meth[rownames(meth) %in% rownames(pheno_M), ]
    rm(meth)

    if (identical(rownames(meth_F), rownames(pheno_F)) & identical(rownames(meth_M), rownames(pheno_M))) {
        cat("\nRownames match for both sexes.\n")
    } else {
        meth_F <- meth_F[match(rownames(pheno_F), rownames(meth_F)),]
        meth_M <- meth_M[match(rownames(pheno_M), rownames(meth_M)),]
        cat("\nRownames have been matched.\n")
    }
} else {
    if (identical(rownames(meth), rownames(pheno))) {
        cat("\nRownames match.\n")
    } else {
        meth <- meth[match(rownames(pheno), rownames(meth)),]
        cat("\nRownames have been matched.\n")
    }
}

cat("\nRAM clean up...\n\n")
gc()


# Add squared CpGs if that has been selected
######################################################

if (cpg_s == "T") { # All CpG^2
    if (sex_s == "T") {
        names_F <- c(colnames(meth_F),paste0(colnames(meth_F),'_2'))
        names_M <- c(colnames(meth_M),paste0(colnames(meth_M),'_2'))
        meth_F <- cbind(meth_F,meth_F^2)
        meth_M <- cbind(meth_M,meth_M^2)
        colnames(meth_F) <- names_F
        colnames(meth_M) <- names_M
        cat("\nSquared CpGs.\n")
    } else {
        names <- c(colnames(meth),paste0(colnames(meth),'_2'))
        meth <- cbind(meth,meth^2)
        colnames(meth) <- names
        cat("\nSquared CpGs.\n")
    }
}

if (cpg_s_subset == "T") { # Only those selected in subset as CpG^2
    if (sex_s == "T") {
        if (cpg_s_age2cpg2 == "T") {
            names_F <- c(colnames(meth_F), paste0(colnames(meth_F[,which(colnames(meth_F) %in% unique(c(rownames(cpg_s_df), rownames(cpg_s_df_cpg2))))], '_2')))
            names_M <- c(colnames(meth_M), paste0(colnames(meth_M[,which(colnames(meth_M) %in% unique(c(rownames(cpg_s_df), rownames(cpg_s_df_cpg2))))], '_2')))
            meth_F <- cbind(meth_F, meth_F[,which(colnames(meth_F) %in% unique(c(rownames(cpg_s_df), rownames(cpg_s_df_cpg2))))]^2)
            meth_M <- cbind(meth_M, meth_M[,which(colnames(meth_M) %in% unique(c(rownames(cpg_s_df), rownames(cpg_s_df_cpg2))))]^2)
            colnames(meth_F) <- names_F
            colnames(meth_M) <- names_M
        } else {
            names_F <- c(colnames(meth_F), paste0(colnames(meth_F[,which(colnames(meth_F) %in% rownames(cpg_s_df))], '_2')))
            names_M <- c(colnames(meth_M), paste0(colnames(meth_M[,which(colnames(meth_M) %in% rownames(cpg_s_df))], '_2')))
            meth_F <- cbind(meth_F, meth_F[,which(colnames(meth_F) %in% rownames(cpg_s_df))]^2)
            meth_M <- cbind(meth_M, meth_M[,which(colnames(meth_M) %in% rownames(cpg_s_df))]^2)
            colnames(meth_F) <- names_F
            colnames(meth_M) <- names_M            
        }
        cat("\nSquared subset of CpGs.\n")
    } else {
        names <- c(colnames(meth), paste0(colnames(meth[,which(colnames(meth) %in% rownames(cpg_s_df))]), '_2'))
        meth <- cbind(meth, meth[,which(colnames(meth) %in% rownames(cpg_s_df))]^2)
        colnames(meth) <- names
        cat("\nSquared subset of CpGs.\n")
    }
} 

cat("\nRAM clean up...\n\n")
gc()


# Create big matrix object for biglasso
######################################################

cat("\nMaking big matrix objects for big lasso...\n")
div <- 10 # Number of chunks to divide OG methylation dataframe
if (sex_s == "T") {

    por_F <- ceiling(length(colnames(meth_F))/div)
    por_M <- ceiling(length(colnames(meth_M))/div)
    
    chunk_list_F <- list()
    chunk_list_M <- list()

    # F
    for (i in 1:div) {
        cat(paste0("\nWorking on chunk (F): ", i, " of ", div))
        if (i == 1) {
            chunk <- as.big.matrix(meth_F[,1:(por_F-1)])
        } else if (i == div) {
            chunk <- as.big.matrix(meth_F[,(por_F*(i-1)):length(colnames(meth_F))])
        } else {
            chunk <- as.big.matrix(meth_F[,(por_F*(i-1)):((por_F*i)-1)])
        }
        cat("\nMade chunk. Appending to F chunk list...\n")
        chunk_list_F <- append(chunk_list_F, list(chunk))

        cat("\nRAM clean up...\n\n")
        gc()
    }
    
    # Saving names prior to chunk fusing
    names_F <- colnames(meth_F)
    rm(meth_F)
    gc()   
    
    cat("\nFusing F chunks!\n")
    x_F <- cbindBM_list(x = chunk_list_F)
    rm(chunk, chunk_list_F)


    # M
    for (i in 1:div) {
        cat(paste0("\nWorking on chunk (M): ", i, " of ", div))
        if (i == 1) {
            chunk <- as.big.matrix(meth_M[,1:(por_M-1)])
        } else if (i == div) {
            chunk <- as.big.matrix(meth_M[,(por_M*(i-1)):length(colnames(meth_M))])
        } else {
            chunk <- as.big.matrix(meth_M[,(por_M*(i-1)):((por_M*i)-1)])
        }
        cat("\nMade chunk. Appending to M chunk list...\n")
        chunk_list_M <- append(chunk_list_M, list(chunk))

        cat("\nRAM clean up...\n\n")
        gc()
    }

    # Saving names prior to chunk fusing
    names_M <- colnames(meth_M)
    rm(meth_M)
    gc()  

    cat("\nFusing M chunks!\n")
    x_M <- cbindBM_list(x = chunk_list_M)
    rm(chunk, chunk_list_M)
     
} else {

    por <- ceiling(length(colnames(meth))/div)
    chunk_list <- list()

    for (i in 1:div) {
        cat(paste0("\nWorking on chunk: ", i, " of ", div))
        if (i == 1) {
            chunk <- as.big.matrix(meth[,1:(por-1)])
        } else if (i == div) {
            chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
        } else {
            chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
        }
        cat("\nMade chunk. Appending to chunk list...\n")
        chunk_list <- append(chunk_list, list(chunk))
        gc()
    }

    # Saving names prior to chunk fusing
    names <- colnames(meth)
    rm(meth)
    
    cat("\nRAM clean up...\n\n")
    gc()

    cat("\nFusing chunks!\n\n")
    if ((cpg_s == "T")) {
        fn <- paste0(opt$out, "temp")
        fn_d <- paste0(opt$out, "temp.desc")
        if (file.exists(fn)) {
            file.remove(fn)
            file.remove(fn_d)
        }
        x <- cbindBM_list(x = chunk_list, backingfile="temp", backingpath = opt$out) # File backed big matrix to get around limits in /dev/shm
    } else {
        x <- cbindBM_list(x = chunk_list)
    }
    
    rm(chunk, chunk_list)
}

cat("\nRAM clean up...\n\n")
gc()


# Stuff for biglasso
######################################################

nfolds <- length(unique(pheno[["Fold"]]))
if (rand == "T") {
    if (length(rownames(pheno)) < 10000) { #w1w3
        nfolds <- 10
    } else if ( (length(rownames(pheno)) > 10000) & (length(rownames(pheno)) < 15000)) { #w1w3 + LBC + GEO
        nfolds <- 15
    } else if ( (length(rownames(pheno)) > 15000) & (length(rownames(pheno)) < 20000)) { #w1w3w4
        nfolds <- 20
    } else if (length(rownames(pheno)) > 20000) { #w1w3w4 + LBC + GEO
        nfolds <- 25
    }
}
cat(paste0("\nTotal number of folds: ", nfolds,"\n"))

# y
if (sex_s == "T") {
    y_F <- pheno_F$age
    y_M <- pheno_M$age
} else {
    y <- as.numeric(pheno$age)
}

# Fold
if (ext == "T") { # If external, make sure folds are ascending with no gaps after filtering for certain cohorts
    count <- 1
    if (sex_s == "T") {
        for (fold in sort(unique(pheno$Fold))) {
            pheno_F[pheno_F$Fold == fold, "Fold"] <- count
            pheno_M[pheno_M$Fold == fold, "Fold"] <- count
            count <- count + 1            
        }
    } else {
        for (fold in sort(unique(pheno$Fold))) {
            pheno[pheno$Fold == fold, "Fold"] <- count
            count <- count + 1
        }
    }
}

if (sex_s == "T") {
    folds_F <- pheno_F$Fold
    folds_M <- pheno_M$Fold
} else {
    folds <- pheno$Fold
}

if (sex_s == "F") {
    cat(paste0("x dimensions:\n"))
    cat(dim(x))
    cat("\n")
    cat(paste0("y dimensions:\n"))
    cat(length(y))
}

# Lasso/elnet parameters
penalty <- "enet"
alpha <- 0.5
if (lasso == "T") {
    penalty <- "lasso"
    alpha <- 1
}

cat("\n\nRAM clean up...\n\n")
gc()


# Get CV'd lambda and obtain effects of each CpG
######################################################

cat("\nFitting elastic net model for age...\n")
if (sex_s == "T") {
    # Obtain lambda
    if (rand == "T") {
        cat("\nMode selected: sex-stratified and randomized folds.\n")
        cvfit_F <- cv.biglasso(x_F, y_F, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, nfolds = nfolds, penalty = penalty)
        cvfit_M <- cv.biglasso(x_M, y_M, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, nfolds = nfolds, penalty = penalty)
    } else {
        cat("\nMode selected: sex-stratified.\n")
        cvfit_F <- cv.biglasso(x_F, y_F, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, cv.ind = folds_F, nfolds = nfolds, penalty = penalty)
        cvfit_M <- cv.biglasso(x_M, y_M, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, cv.ind = folds_M, nfolds = nfolds, penalty = penalty)
    }
    lambda_F <- cvfit_F$lambda.min # Get lambda which gives minimum MSE
    lambda_M <- cvfit_M$lambda.min
    # Obtain coeffs
    fit_F <- biglasso(x_F, y_F, family = "gaussian", alpha = alpha, ncores = 8, lambda = lambda_F, penalty = penalty)
    fit_M <- biglasso(x_M, y_M, family = "gaussian", alpha = alpha, ncores = 8, lambda = lambda_M, penalty = penalty)
} else {
    # Obtain lambda
    if (rand == "T") {
        cat("\nMode selected: randomized folds.\n")
        cvfit <- cv.biglasso(x, y, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, nfolds = nfolds, penalty = penalty)
    } else {
        cvfit <- cv.biglasso(x, y, family = 'gaussian', seed = seed, alpha = alpha, ncores = 8, cv.ind = folds, nfolds = nfolds, penalty = penalty)
    }
    lambda <- cvfit$lambda.min # Get lambda which gives minimum MSE
    # Obtain coeffs
    fit <- biglasso(x, y, family = "gaussian", alpha = alpha, ncores = 8, lambda = lambda, penalty = penalty)
}

cat("\nRAM clean up...\n\n")
gc()


# Get the good stuff!
######################################################

cat("\nNow extracting info of interest.\n")
if (sex_s == "T") {
    # F
    coefs_F <- coef(fit_F) # Get betas
    coefs_F <- as.data.frame(coefs_F[which(coefs_F!=0),]) # Remove coeficients that are 0 
    names(coefs_F)[1] <- "Coefficient"
    coefs_F["CpG_loc"] <- noquote(sub('V', '', rownames(coefs_F)))
    coefs_F[1, "CpG_loc"] <- NA
    coefs_F[["CpG_loc"]] <- as.numeric(coefs_F[["CpG_loc"]])
    coefs_F[2:nrow(coefs_F), "Variable"] <- names_F[coefs_F[2:nrow(coefs_F), "CpG_loc"]]
    coefs_F[1, "Variable"] <- "Intercept"
    coefs_F <- coefs_F[c("Variable", "Coefficient")]
    # M
    coefs_M <- coef(fit_M) # Get betas
    coefs_M <- as.data.frame(coefs_M[which(coefs_M!=0),]) # Remove coeficients that are 0 
    names(coefs_M)[1] <- "Coefficient"
    coefs_M["CpG_loc"] <- noquote(sub('V', '', rownames(coefs_M)))
    coefs_M[1, "CpG_loc"] <- NA
    coefs_M[["CpG_loc"]] <- as.numeric(coefs_M[["CpG_loc"]])
    coefs_M[2:nrow(coefs_M), "Variable"] <- names_M[coefs_M[2:nrow(coefs_M), "CpG_loc"]]
    coefs_M[1, "Variable"] <- "Intercept"
    coefs_M <- coefs_M[c("Variable", "Coefficient")]
} else {
    coefs <- coef(fit) # Get betas
    coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
    names(coefs)[1] <- "Coefficient"
    coefs["CpG_loc"] <- noquote(sub('V', '', rownames(coefs)))
    coefs[1, "CpG_loc"] <- NA
    coefs[["CpG_loc"]] <- as.numeric(coefs[["CpG_loc"]])
    coefs[2:nrow(coefs), "Variable"] <- names[coefs[2:nrow(coefs), "CpG_loc"]]
    coefs[1, "Variable"] <- "Intercept"
    coefs <- coefs[c("Variable", "Coefficient")]
}


# Export coeffs
######################################################

cat("\nExporting!\n\n")

# Export
if (sex_s == "T") {
    filename_F <- paste0(opt$out, "elnet_coefficients_", opt$name, squared_f, squared_subset_f, squared_cpg_f, random_f, lasso_f, "_F.tsv")
    filename_M <- paste0(opt$out, "elnet_coefficients_", opt$name, squared_f, squared_subset_f, squared_cpg_f, random_f, lasso_f, "_M.tsv")
    write.table(coefs_F, file = filename_F, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(coefs_M, file = filename_M, row.names = FALSE, sep = "\t", quote = FALSE)
} else {
    filename <- paste0(opt$out, "elnet_coefficients_", opt$name, squared_f, squared_subset_f, squared_cpg_f, random_f, lasso_f, ".tsv")
    #filename <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/elnet_train/w1w3w4/random_noscalesample_noscalecpg_subset20K/elnet_coefficients_random_noscalesample_noscalecpg_subset20K_squaredsubset900.tsv"
    write.table(coefs, file = filename, row.names = FALSE, sep = "\t", quote = FALSE)
}

# Erase temp and temp.desc files
if (cpg_s == "T") {
    fn <- paste0(opt$out, "temp")
    fn_d <- paste0(opt$out, "temp.desc")
    if (file.exists(fn)) {
        file.remove(fn)
        file.remove(fn_d)
    }
}