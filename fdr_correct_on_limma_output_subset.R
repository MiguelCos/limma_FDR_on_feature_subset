### Script for two-conditions limma analysis (i.e. Tumor vs. Normal) and 
### apply multiple-testing correction on a subset of features
## Miguel Cosenza v0.1 - last edit 11.10.2021

### 1. Please provide a meaningful name for the dataset/experiment ----

exper_code <- "My_experiment_123"

### 2. Would you like to run a robust regression? ----

robust <- TRUE # or FALSE if not.

### 3. Define the interesting feature type
# Example: "semi_Nterm"
# this will dependend on your annotation_features.txt

interesting_features <- "semi_Nterm"

## Required packages ----

packages <- c("dplyr", "here", "tidyr", "ggplot2", "rmarkdown", "knitr", "reshape", "tibble")

biopackgs <- c("limma")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
                    install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
                    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
                    
                    BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
                    
}

library(dplyr)
library(tibble)
library(stringr)
library(limma)
library(rmarkdown)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(reshape)
library(plotly)

## Load data ----
expr_dat <- read.delim("data/input_limma.txt") %>%
                    dplyr::mutate(ID = stringr::str_remove_all(ID, ";.*$"))

annot_dat <- read.delim("data/annotation.txt")

annot_feat <- read.delim("data/annotation_features.txt")

## Define the design matrix ----

## Defining the design matrix itself

groups <- as.factor(annot_dat$Group)

design <- model.matrix(~0+groups)

row.names(design) <- annot_dat$Sample_ID

colnames(design) <- colnames(design) %>% str_remove(., "groups") %>% str_trim()

# 3. DEFINE WHICH GROUPS YOU WISH TO COMPARE ----

contrast.matrix <- makeContrasts(x1_7-x1_50, levels=design) # MODIFY THIS LINE

## Prep expression data matrix ----
n_contrasts <- dim(contrast.matrix)[2]

expr_dat[1:5,1:5]
design %>% head()
row.names(design) %>% head()

tomat <- dplyr::select(expr_dat, row.names(design))

row.names(tomat) <- expr_dat$ID

mat <- as.matrix(tomat)

## Execute the linear model / limma ----

if(robust == TRUE){
                    regression <- "robust"
} else {
                    regression <- "ls"
}

fit <- lmFit(mat, design = design, method = regression)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

output_limma2 <- topTable(fit2, adjust.method = "BH", number = Inf) %>%
                    rownames_to_column("ID")

# merge limma output with feature annotation  ----

tab_limma_feature_annot <- left_join(output_limma2, annot_feat,
                                     by = "ID")

# filter limma table based on interesting feature ----

tab_limma_subsetin <- filter(tab_limma_feature_annot,
                           feature_category == interesting_features) %>%
# apply FDR correction only on the subset
                    mutate(adj.P.Val = p.adjust(p = P.Value, method = "BH"))

# create a table of non-interesting features containing the globally adjusted p-values
tablimma_subsetout <- filter(tab_limma_feature_annot,
                             feature_category != interesting_features)

output_limma3 <- bind_rows(tab_limma_subsetin,
                           tablimma_subsetout) %>% 
                    # mark/label the features' FDR/Adjusted p-values as 'global' or 'feature-specific'
                    mutate(fdr_correction = if_else(feature_category == interesting_features,
                                                    true = 'feature-specific',
                                                    false = 'global'))
                    
## Generate output ----

tab_out <- paste0("Output/tab_output_paired_DE_analysis_limma_",exper_code,
                  ".txt")

write.table(output_limma3,
            file = tab_out,
            row.names = FALSE, col.names = TRUE)

sig_hits <- dplyr::filter(output_limma3, 
                          adj.P.Val <= 0.05) %>% row.names(.)

sig_hits_sel <- dplyr::filter(output_limma3, 
                              adj.P.Val <= 0.05,
                              fdr_correction == "feature-specific") %>% row.names(.)

n_significant_total <- length(sig_hits)

n_significant_sel_features <- length(sig_hits_sel)

## Prep volcano plots for contrasts ----


tovolc_all <- output_limma3 %>% 
                    mutate(Differentially_expressed = case_when(adj.P.Val <= 0.05 ~ TRUE,
                                                                            TRUE ~ FALSE))

tovolc_featu <- output_limma3 %>% 
                    filter(fdr_correction == "feature-specific") %>%
                    mutate(Differentially_expressed = case_when(adj.P.Val <= 0.05 ~ TRUE,
                                                                TRUE ~ FALSE))
volcanoes <- ggplot(data = tovolc_all, 
                    mapping = aes(x = logFC,
                                  y = -log10(adj.P.Val),
                                  color = Differentially_expressed,
                                  label = ID)) + 
                    geom_point(alpha = 0.5) + 
                    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
                    ggtitle("Volcano plot of all limma-tested features") + 
                    theme(legend.position = "none")

volcano_features <- ggplot(data = tovolc_featu, 
                    mapping = aes(x = logFC,
                                  y = -log10(adj.P.Val),
                                  color = Differentially_expressed,
                                  label = ID)) + 
                    geom_point(alpha = 0.5) + 
                    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
                    ggtitle("Volcano plot of selected/interesting limma-tested features") + 
                    theme(legend.position = "none")


print(volcanoes)
print(volcano_features)

rep_out <- paste0("Output/simple_report_output_paired_DE_analysis_limma_",exper_code,
                  ".html")

rmarkdown::render("renderReport.R", 
                  output_file = rep_out,
                  output_dir = "Output")

