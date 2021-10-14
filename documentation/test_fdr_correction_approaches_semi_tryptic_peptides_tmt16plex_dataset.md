TMT16plex titration dataset - FragPipe search processed on Percolator -
Evaluate performance of FDR correction on subset
================
Miguel Cosenza
11 October, 2021

-   [1 Load data](#load-data)
-   [2 Questions/aims:](#questionsaims)
-   [3 Processing](#processing)
    -   [3.1 Annotate peptides](#annotate-peptides)
    -   [3.2 Label peptide with free N-terminus (no TMT at
        N-terminus)](#label-peptide-with-free-n-terminus-no-tmt-at-n-terminus)
    -   [3.3 Merge TMT-integrator report with annotated peptide
        IDs](#merge-tmt-integrator-report-with-annotated-peptide-ids)
    -   [3.4 Label peptides as E.coli or
        human](#label-peptides-as-ecoli-or-human)
-   [4 Summarized count results](#summarized-count-results)
-   [5 Distribution of E.coli quant per feature
    category](#distribution-of-ecoli-quant-per-feature-category)
-   [6 Differential Abundance analysis
    performance](#differential-abundance-analysis-performance)
    -   [6.0.1 Prep abundance matrix and design matrix and
        comparisons](#prep-abundance-matrix-and-design-matrix-and-comparisons)
    -   [6.0.2 Prep annotation table mapping peptide/index to
        specificity
        category](#prep-annotation-table-mapping-peptideindex-to-specificity-category)
    -   [6.0.3 Prep PSM count information for weighted correction and
        DEqMS](#prep-psm-count-information-for-weighted-correction-and-deqms)
    -   [6.1 Approach 1](#approach-1)
    -   [6.2 Approach 2](#approach-2)
    -   [6.3 Approach 3](#approach-3)
    -   [6.4 Approach 4](#approach-4)
    -   [6.5 Approach 5](#approach-5)
    -   [6.6 Approach 5.2](#approach-52)
    -   [6.7 Approach 6](#approach-6)
    -   [6.8 Approach 7](#approach-7)
    -   [6.9 Approach 8](#approach-8)
    -   [6.10 Approach 9](#approach-9)
    -   [6.11 Approach 10](#approach-10)
-   [7 Distribution of adjusted p-values per
    approach](#distribution-of-adjusted-p-values-per-approach)
    -   [7.1 Observations:](#observations)

# 1 Load data

``` r
tmt_75_pept <- read_tsv(here("data/tmt16plex_variable-Nter_0-75_peptide_lev_norm_fp_v16_msfrag3.3_percolator_shuffled_tag/tmt-report/abundance_peptide_MD.tsv")) %>%
                    clean_names()

tmt_75_prot <- read_tsv(here("data/tmt16plex_variable-Nter_0-75_peptide_lev_norm_fp_v16_msfrag3.3_percolator_shuffled_tag/tmt-report/abundance_protein_MD.tsv")) %>% 
          clean_names()

ids_75_pept <- read_tsv(here("data/tmt16plex_variable-Nter_0-75_peptide_lev_norm_fp_v16_msfrag3.3_percolator_shuffled_tag/MC-FragPipe/peptide.tsv")) %>% 
          clean_names() %>% 
          mutate(Peptide = peptide, Genes = protein_id)

ids_75_prot <- read_tsv(here("data/tmt16plex_variable-Nter_0-75_peptide_lev_norm_fp_v16_msfrag3.3_percolator_shuffled_tag/MC-FragPipe/protein.tsv")) %>% 
          clean_names()

psms_75 <- read_tsv(here("data/tmt16plex_variable-Nter_0-75_peptide_lev_norm_fp_v16_msfrag3.3_percolator_shuffled_tag/MC-FragPipe/psm.tsv")) %>% 
          clean_names()

fasta75 <- read.fasta(here("data/tmt16plex_variable-Nter_0-75_peptide_lev_norm_fp_v16_msfrag3.3_percolator_shuffled_tag/MC-FragPipe/protein.fas"),
  seqtype = "AA", as.string = TRUE)
```

# 2 Questions/aims:

This report is intended to explore and show how to extract (filter-in)
N-terminal peptides from FragPipe and annotate as potentially canonical
or non-canonical based on Uniprot processing information.

#### 2.0.0.1 Objectives:

-   Extract N-terminal peptides from FragPipe output.
-   Extract annotated processing information from Uniprot.  
-   Use this information to mark the identified N-terminal peptides as
    potentially canonical or not.

# 3 Processing

## 3.1 Annotate peptides

``` r
## PIF 0.75 

if(!file.exists(here("results/peptide_annotation_pif75_32search.tsv"))){
            
    peptide_annot_75 <- annotate_peptides(ids_75_pept, fasta75)
    
    write_tsv(peptide_annot_75, here("results/peptide_annotation_pif75_32search.tsv"))
            
} else {
            
    peptide_annot_75 <- read_tsv(here("results/peptide_annotation_pif75_32search.tsv"))
            
}

protein_length <- ids_75_prot %>%
  dplyr::select(protein_id, length)

pept_75_annot <- left_join(ids_75_pept, peptide_annot_75,
                           by = c("protein_id", "Peptide")) %>%
  left_join(., protein_length)
```

## 3.2 Label peptide with free N-terminus (no TMT at N-terminus)

``` r
nterannot_75 <- pept_75_annot %>% 
  mutate(nterm = case_when(str_detect(assigned_modifications, "N-term\\(304.2072\\)") ~ "TMT-labelled",
                           str_detect(assigned_modifications, "N-term\\(42.0106\\)") ~ "acetylated",
                           TRUE ~ "free")) %>%
  mutate(tmt_tag = case_when(str_detect(assigned_modifications, "N-term\\(304.2072\\)") ~ "nterm",
                           str_detect(assigned_modifications, "K\\(304.2072\\)") ~ "lysine",
                           str_detect(assigned_modifications, "K\\(304.2072\\)",
                                      negate = TRUE) & str_detect(assigned_modifications, "N-term\\(42.0106\\)") ~ "untagged_acetylated",
                           str_detect(assigned_modifications, "K\\(304.2072\\)",
                                      negate = TRUE) & nterm == "acetylated" ~ "untagged_acetylated",
                           str_detect(assigned_modifications, "K\\(304.2072\\)",
                                      negate = TRUE) & nterm == "free" ~ "untagged_free",
                           TRUE ~ "untagged")) %>%
  mutate(specificity = case_when(nterm == "acetylated" & str_detect(last_aa, "R|K") ~ "tryptic",
                                 nterm == "acetylated" & str_detect(last_aa, "R|K", negate = TRUE) ~ "unspecific",
                                 TRUE ~ specificity),
         semi_type  = case_when(nterm == "acetylated" & str_detect(last_aa, "R|K") ~ "tryptic_nterm",
                                nterm == "acetylated" & str_detect(last_aa, "R|K", negate = TRUE) ~ "unspecific_nterm",
                                 TRUE ~ semi_type)) %>% 
  dplyr::select(-c(matches("hek_"), matches("x1"))) # eliminate columns with quant information (not normalized)
```

## 3.3 Merge TMT-integrator report with annotated peptide IDs

``` r
## PIF 0.75
tmt_reprt_75_annot <- left_join(tmt_75_pept, nterannot_75)

excluded_tmt_reprt_75_annot <- anti_join(nterannot_75, tmt_75_pept)
```

## 3.4 Label peptides as E.coli or human

``` r
zero2na <- function(x){
          y <- ifelse(x == 0,
                      yes = NA,
                      no = x)
          
          return(y)
}

#zero2na(0.223) 

prequant_75 <- tmt_reprt_75_annot %>% 
          mutate(human_ecoli = ifelse(test = str_detect(entry_name, "HUMAN"),
                                     yes = "human",
                                     no = "ecoli")) %>%
          mutate(across(matches("x1_"),zero2na)) %>%
          dplyr::select(peptide, matches("^x1_"), peptide, specificity, nterm, semi_type, tmt_tag, human_ecoli) %>%
          na.omit() %>%
          rowwise() %>%
          mutate(`FC_0.06/0.02` = mean(x1_17_8, x1_17_10, x1_17_12, x1_17_14, na.rm = TRUE)/mean(x1_50_15, x1_50_17, x1_50_19, x1_50_21, na.rm = TRUE),
                 `FC_0.14/0.06` = mean(x1_7_16, x1_7_18, x1_7_20,x1_7_22, na.rm = TRUE)/mean(x1_17_8, x1_17_10, x1_17_12, x1_17_14, na.rm = TRUE),
                 `FC_0.14/0.02` = mean(x1_7_16, x1_7_18, x1_7_20,x1_7_22,na.rm = TRUE)/mean(x1_50_15, x1_50_17, x1_50_19, x1_50_21, na.rm = TRUE))


fc_quant_eval1_75 <- dplyr::select(prequant_75,
                                peptide, specificity, nterm, semi_type, tmt_tag, human_ecoli,
                                `FC_0.06/0.02`, `FC_0.14/0.06`,`FC_0.14/0.02`) %>% 
          pivot_longer(cols = c(`FC_0.06/0.02`, `FC_0.14/0.06`,`FC_0.14/0.02`),
                       values_to = "FCs",
                       names_to = "Ratios") %>%
          pivot_longer(cols = c(specificity, nterm, semi_type, tmt_tag),
                       values_to = "feature_type",
                       names_to = "category") %>%
          mutate(FCs = ifelse(is.nan(FCs),
                              yes = NA,
                              no = FCs),
                 Expected_ratio = str_remove(Ratios, "^FC_")) %>%
          mutate(Expected_FC = eval(parse(text = Expected_ratio))) %>%
          mutate(log2_FC = log2(FCs),
                 log2_Expected_FC = log2(Expected_FC))

df2 <- fc_quant_eval1_75 %>%
  group_by(feature_type, Ratios) %>%
  summarise(Expected = log2_Expected_FC)
```

``` r
# PIF 0.75
to_count_info_75 <- prequant_75 %>% 
          dplyr::select(peptide, specificity, nterm, semi_type, tmt_tag)

to_count_info_acetyl <- to_count_info_75 %>% 
  filter(nterm == "acetylated")
```

# 4 Summarized count results

``` r
## PIF 75
n_semi <- dplyr::count(to_count_info_75, specificity) %>% 
          dplyr::rename(feature_type = specificity) %>%
          dplyr::mutate(category = "specificity")

n_term <- dplyr::count(to_count_info_75, nterm) %>% 
          dplyr::rename(feature_type = nterm) %>%
          dplyr::mutate(category = "N-term")

n_semi_type <- dplyr::count(to_count_info_75, semi_type) %>% 
               dplyr::rename(feature_type = semi_type) %>%
          dplyr::mutate(category = "Semi type")

n_tmt_tag <- dplyr::count(to_count_info_75, tmt_tag) %>% 
               dplyr::rename(feature_type = tmt_tag)  %>%
          dplyr::mutate(category = "TMT location")

n_total <- tibble(feature_type = "Total",
                  n = nrow(to_count_info_75),
                  category = "Total")
          
summary_count_75 <- bind_rows(n_semi,
                              n_term,
                              n_semi_type,
                              n_tmt_tag,
                              n_total)
```

# 5 Distribution of E.coli quant per feature category

These quant values were extracted based on the list of peptides
summarized after TMT-integrator filtering.

``` r
fc_quant_eval2_75 <- mutate(fc_quant_eval1_75,
                            feature_type = factor(feature_type,
                                                  levels = c("tryptic","semi_tryptic", "unspecific",
                                                          "semi_Cterm", "semi_Nterm",
                                                          "TMT-labelled", "acetylated","free",
                                                          "lysine", "nterm")),
                            category = factor(category,
                                              levels = c("specificity", "semi_type",
                                                        "tmt_tag", "nterm")))

df2 <- fc_quant_eval2_75 %>% # this is generating a df that can be used to plot the expected log2 ratio
  group_by(feature_type, Ratios) %>%
  summarise(Expected = log2_Expected_FC)
```

``` r
ggplot(fc_quant_eval2_75,
       aes(x = log2_FC, fill = human_ecoli)) + 
  geom_density() +
  #geom_vline(data = df2, mapping = aes(xintercept = Expected), color = "red",
  #           linetype = "dashed") +
  facet_grid(category+feature_type~Ratios) +
  labs(title = "Distribution of Fold-changes",
      subtitle = "Visualized by feature type") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10),
               axis.text.y = element_text(hjust = 0.5, size = 10),
               panel.background = element_blank(),
               panel.grid.major = element_line(color = "grey"),
               panel.border = element_rect(colour = "black", fill=NA, size=1.5),
               axis.title=element_text(size=12,face="bold"))
```

![](test_fdr_correction_approaches_semi_tryptic_peptides_tmt16plex_dataset_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# 6 Differential Abundance analysis performance

This section is intended to evaluate the performance for differential
abundant semi-tryptic peptides based on different approaches.

1.  Limma Robust on whole dataset + Global FDR correction on all
    peptides.
2.  Limma Robust on whole dataset + FDR only on semi-tryptic peptides.
3.  Limma Robust on whole dataset + weighted FDR correction on all
    peptides.
4.  Limma Robust on whole dataset + weighted FDR on semi-tryptic
    peptides.
5.  Subset semi-tryptic peptides + Limma Robust (Global FDR correction,
    on the subset)
6.  DEqMS on whole dataset + Global FDR correction.
7.  DEqMS on whole dataset + FDR only on semi-tryptic peptides.
8.  DEqMS on whole dataset + weighted FDR correction on all peptides.
9.  DEqMS on whole dataset + weighted FDR on semi-tryptic peptides.
10. Subset semi-tryptic peptides + DEqMS (Global FDR correction, on the
    subset).

### 6.0.1 Prep abundance matrix and design matrix and comparisons

For this test, I will focus on the comparison that yield the biggest
expected theoretical differences: 0.16 vs 0.02 (*E. coli*:HEK ratios).

The subset of peptides will be based on N-term peptides because these
show the best separation in the density plots.

``` r
tomat <- tmt_75_pept %>%
                    dplyr::select(index, where(is.numeric)) %>%
                    dplyr::select(-c(max_pep_prob, reference_intensity)) %>%
                    dplyr::select(-starts_with("hek_only")) %>%
                    dplyr::select(index, starts_with("x1_7"), starts_with("x1_50")) %>% 
                    na.omit()

mat <- tomat %>%
                    column_to_rownames("index") %>%
                    as.matrix()

ratio <- colnames(mat) %>%
                    str_sub(end = -4)

design <- model.matrix(~ratio)
row.names(design) <- colnames(mat)
```

### 6.0.2 Prep annotation table mapping peptide/index to specificity category

An annotation table is necessary for this.

``` r
index2semitype <- tmt_reprt_75_annot %>%
                    dplyr::select(index, peptide, semi_type)
```

### 6.0.3 Prep PSM count information for weighted correction and DEqMS

``` r
peptide2spectralcount <- ids_75_pept %>%
                    dplyr::select(peptide, spectral_count)
                    
annotation_correction <- left_join(index2semitype, peptide2spectralcount)
```

#### 6.0.3.1 Distribution of spectral counts per peptide

``` r
ggplot(annotation_correction,
       aes(x = spectral_count, fill = semi_type)) + 
                    geom_histogram(binwidth = 2)
```

![](test_fdr_correction_approaches_semi_tryptic_peptides_tmt16plex_dataset_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Semi tryptic peptides tend to have lower spectral counts, represented by
the higher amount of peptides with count 1-2.

## 6.1 Approach 1

Limma Robust on whole dataset + Global FDR correction on all peptides.

``` r
fit1 <- lmFit(mat, # expression matrix
              design = design, # design matrix
              method = "robust", # set to robust regression  
              maxit = 50) 

fit1 <- eBayes(fit1)

fit_tab1 <- topTable(fit1, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "Limma",
         correction = "bh_global") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")

fit_tab1_2 <- left_join(fit_tab1, annotation_correction)
```

``` r
fit_tab1_3 <- filter(fit_tab1_2,
                     semi_type == "semi_Nterm") 
```

## 6.2 Approach 2

Limma Robust on whole dataset + FDR only on semi-tryptic peptides.

``` r
fit_tab2 <- topTable(fit1, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "Limma",
         correction = "bh_subset") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")

fit_tab2_2 <- left_join(fit_tab2, annotation_correction)
```

#### 6.2.0.1 Apply adjustment only on subset

``` r
fit_tab2_3 <- filter(fit_tab2_2,
                     semi_type == "semi_Nterm") %>%
                    mutate(adj.P.Val = p.adjust(p = P.Value, method = "BH"))
```

## 6.3 Approach 3

Limma Robust on whole dataset + weighted FDR global correction.

``` r
fit_tab3 <- topTable(fit1, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "Limma",
         correction = "ihw_global") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")

fit_tab3_2 <- left_join(fit_tab3, annotation_correction)
```

``` r
ihw_corr_3 <- ihw(P.Value ~ spectral_count,  data = fit_tab3_2, alpha = 0.05)
```

``` r
ihw_corr_3_df <- dplyr::rename(ihw_corr_3@df,
                               adj.P.Val = adj_pvalue,
                               P.Value = pvalue)

fit_tab3_21 <- fit_tab3_2 %>% dplyr::select(-adj.P.Val)

fit_tab3_22 <- left_join(fit_tab3_21, ihw_corr_3_df, by = "P.Value")

fit_tab3_3 <- filter(fit_tab3_22,
                     semi_type == "semi_Nterm")

fit_tab3_4 <- dplyr::select(fit_tab3_3,
                            names(fit_tab2_3))
```

## 6.4 Approach 4

Limma Robust on whole dataset + weighted FDR on semi-tryptic peptides.

``` r
fit_tab4_2 <- fit_tab3_2 %>%
                    mutate(correction = "ihw_subset") %>%
                    filter(semi_type == "semi_Nterm")
```

``` r
ihw_corr_4 <- ihw(P.Value ~ spectral_count,  data = fit_tab4_2, alpha = 0.05)
```

``` r
fit_tab4_3 <- dplyr::select(fit_tab4_2,
                            names(fit_tab2_3))
```

## 6.5 Approach 5

Subset semi-tryptic peptides + Limma Robust (Global FDR correction, on
the subset)

**Filter abundance matrix**

``` r
semi_peptides <- annotation_correction %>%
                    filter(semi_type == "semi_Nterm") %>%
                    pull(index)

tomat_fil <- tomat %>%
                    filter(index %in% semi_peptides)

mat_fil <- tomat_fil %>%
                    column_to_rownames("index") %>%
                    as.matrix()
```

**Limma fit**

``` r
fit5 <- lmFit(mat_fil, # expression matrix
              design = design, # design matrix
              method = "robust", # set to robust regression  
              maxit = 50) 

fit5 <- eBayes(fit5)

fit_tab5 <- topTable(fit5, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "Limma",
         correction = "prefiltered_bh") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")

fit_tab5_2 <- left_join(fit_tab5, annotation_correction)
```

## 6.6 Approach 5.2

Subset semi-tryptic peptides + Limma Robust + Weighted FDR correction
(on the subset)

``` r
fit_tab52_2 <- fit_tab5_2 %>%
                    mutate(correction = "prefiltered_ihw") %>%
                    filter(semi_type == "semi_Nterm")
```

``` r
ihw_corr_52 <- ihw(P.Value ~ spectral_count,  data = fit_tab52_2, alpha = 0.05)
```

## 6.7 Approach 6

DEqMS on whole dataset + Global FDR correction.

``` r
fit6 <- fit1

fit6_pretab <- topTable(fit6, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>%
                    rownames_to_column("index") %>%
                    left_join(., annotation_correction)

fit6$count <- fit6_pretab$spectral_count

fit6 = spectraCounteBayes(fit6)

fit_tab6 <- topTable(fit6, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "deqms",
         correction = "bh_global") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")  

fit_tab6_2 <- left_join(fit_tab6, annotation_correction)
```

``` r
VarianceBoxplot(fit6,n=30,main="Variance boxplots",xlab="PSM count")
```

![](test_fdr_correction_approaches_semi_tryptic_peptides_tmt16plex_dataset_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

DEqMS is developed for protein quant and using PSM counts to weight the
protein variance in limma [DEqMS
vignette](https://bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html)…

In terms of peptides (here), the number of spectral counts does not seem
to have a consistent impact on variance.

``` r
fit_tab6_3 <- filter(fit_tab6_2,
                     semi_type == "semi_Nterm") 
```

## 6.8 Approach 7

DEqMS on whole dataset + FDR only on semi-tryptic peptides.

``` r
fit_tab7 <- topTable(fit6, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "deqms",
         correction = "bh_subset") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")  

fit_tab7_2 <- left_join(fit_tab7, annotation_correction)
```

#### 6.8.0.1 Apply adjustment only on subset

``` r
fit_tab7_3 <- filter(fit_tab7_2,
                     semi_type == "semi_Nterm") %>%
                    mutate(adj.P.Val = p.adjust(p = P.Value, method = "BH"))
```

## 6.9 Approach 8

DEqMS on whole dataset + weighted FDR correction on all peptides.

``` r
fit_tab8 <- topTable(fit6, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "deqms",
         correction = "ihw_global") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")

fit_tab8_2 <- left_join(fit_tab8, annotation_correction)
```

``` r
ihw_corr_8 <- ihw(P.Value ~ spectral_count,  data = fit_tab8_2, alpha = 0.05)
```

``` r
ihw_corr_8_df <- dplyr::rename(ihw_corr_8@df,
                               adj.P.Val = adj_pvalue,
                               P.Value = pvalue)

fit_tab8_21 <- fit_tab8_2 %>% dplyr::select(-adj.P.Val)

fit_tab8_22 <- left_join(fit_tab8_21, ihw_corr_8_df, by = "P.Value")

fit_tab8_3 <- filter(fit_tab8_22,
                     semi_type == "semi_Nterm") 

fit_tab8_4 <- dplyr::select(fit_tab8_3,
                            names(fit_tab7_2))
```

## 6.10 Approach 9

DEqMS on whole dataset + weighted FDR on semi-tryptic peptides.

``` r
fit_tab9_2 <- fit_tab8_2 %>%
                    mutate(correction = "ihw_subset") %>%
                    filter(semi_type == "semi_Nterm")
```

``` r
ihw_corr_9 <- ihw(P.Value ~ spectral_count,  data = fit_tab9_2, alpha = 0.05)
```

``` r
fit_tab9_3 <- dplyr::select(fit_tab9_2,
                            names(fit_tab7_2))
```

## 6.11 Approach 10

Subset semi-tryptic peptides + DEqMS (Global FDR correction, on the
subset).

**DEqMS-Limma fit**

``` r
fit10 <- fit5

fit10_pretab <- topTable(fit10, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>%
                    rownames_to_column("index") %>%
                    left_join(., annotation_correction)

fit10$count <- fit10_pretab$spectral_count

fit10 = spectraCounteBayes(fit10)

fit_tab10 <- topTable(fit10, 
                            coef = "ratiox1_7", 
                            number = Inf, 
                            adjust.method = "BH") %>% # FDR correction of DE analysis happens here
  mutate(algorithm = "deqms",
         correction = "prefiltered_bh") %>% # just adding a column indicating which model was applied.
  rownames_to_column("index")  

fit_tab10_2 <- left_join(fit_tab10, annotation_correction)
```

# 7 Distribution of adjusted p-values per approach

``` r
all_semi_tables <- bind_rows(fit_tab1_3,
                             fit_tab2_3,
                             fit_tab3_4,
                             fit_tab4_3,
                             fit_tab5_2,
                             fit_tab6_3,
                             fit_tab7_3,
                             fit_tab8_4,
                             fit_tab9_3,
                             fit_tab10_2)
```

``` r
ggplot(all_semi_tables,
       aes(x = adj.P.Val, fill = correction)) + 
                    geom_histogram(binwidth = 0.03, position = "dodge") + 
                    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") + 
                    facet_wrap(~algorithm, ncol = 1)
```

![](test_fdr_correction_approaches_semi_tryptic_peptides_tmt16plex_dataset_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

## 7.1 Observations:

Including peptide spectral counts as weights for limma (DEqMS) does not
improve the sensitivity when doing peptide-level differential abundance
analysis.

The sensitivity seems to be better when applying the BH correction on a
subset of p-values.

Fitting the whole `limma` model directly on the subset is not
necessarily better than the global FDR correction on the whole set of
identified values and it is definitely worse than applying FDR
correction on the subset after `limma` on the whole data set.

Using weighted correction is similar than applying global FDR
correction.
