How to generate an `annotation_features.txt` file
================
Miguel Cosenza
14 October, 2021

-   [1 Example based on annotated semi-specific
    peptides](#example-based-on-annotated-semi-specific-peptides)
    -   [1.1 First load the file of annotated
        peptides](#first-load-the-file-of-annotated-peptides)
    -   [1.2 Then select the interesting
        feature](#then-select-the-interesting-feature)
-   [2 Example based on protein names from
    headers](#example-based-on-protein-names-from-headers)

In this document I am including a small example on how to generate a
`annotation_features.txt` file to be used as an input to do ‘selected’
FDR control after `limma` on a interestin features.

**Note**: the outputs that can be generated in this example would not
match the ones in the sample data. This document is intender to provide
a guide on how to approach such a problem of creating the
`annotation_features.txt`. You need to adapt it so the IDs in your
`annotation_features.txt` file would match the ones in your
`input_limma.txt` file.

# 1 Example based on annotated semi-specific peptides

We already have a script that allows us to annotate our peptides
regarding their specificity ([Mapping peptides to proteins in FASTA -
Annotate](https://github.com/MiguelCos/mapping_peptides_to_proteins_from_fasta_file))
and we can use that to create our `annotation_features.txt` file.

## 1.1 First load the file of annotated peptides

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.2     v dplyr   1.0.7
    ## v tidyr   1.1.3     v stringr 1.4.0
    ## v readr   2.0.1     v forcats 0.5.1

    ## Warning: package 'readr' was built under R version 4.1.1

    ## Warning: package 'dplyr' was built under R version 4.1.1

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(here)
```

    ## here() starts at C:/Users/migue/OneDrive/Documentos/R_Projects/7_scripts_workflows/limma_FDR_on_feature_subset

``` r
annot_peptides <- read_tsv(here("results/peptide_annotation_pif75_32searchv2.tsv"))
```

    ## Rows: 78566 Columns: 14

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (11): protein_id, protein_description, Peptide, last_aa, aa_after, aa_be...
    ## dbl  (3): protein_length, start_position, end_position

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

## 1.2 Then select the interesting feature

Take a look at the names of the columns

``` r
names(annot_peptides)
```

    ##  [1] "protein_id"          "protein_description" "protein_length"     
    ##  [4] "Peptide"             "start_position"      "end_position"       
    ##  [7] "last_aa"             "aa_after"            "aa_before"          
    ## [10] "following_10_resid"  "previous_10_resid"   "previous_all_resid" 
    ## [13] "semi_type"           "specificity"

In this case, the columns related to the unique ID of the Peptide +
Protein (index) and the type of specificity.

**First we create our `index` variable so it matches the one in the
`input_limma.txt` file.**

``` r
annot_peptides_2 <- mutate(annot_peptides, # annotated peptides table as input
                           index = paste(protein_id, Peptide, "_")) # new variable index is created by concatenating protein ID and peptide.
```

**Then we select the interesting columns**

In this case: index and semi_type (containing the specificity
information per peptide/feature)

``` r
annotation_features <- dplyr::select(annot_peptides_2,
                                     ID = index, # changed index to ID column name
                                     feature_category = semi_type)  # changed semi_type to feature_category column name
```

We can now save the `annotation_features.txt` file

``` r
write_delim(annotation_features,
            file = here("sample/data/annotation_features_semi_try.txt"), 
            delim = "\t")
```

# 2 Example based on protein names from headers

In the case of metaproteomics studies, we might be interested in
labelling proteins as ‘human’ or ‘non-human/bacteria’.

If we have the complete Uniprot identifier of the proteins, we can do
that annotation easily.

We load a sample dataset containing protein headers/IDs in the first
column

``` r
proteins_w_headers <- read_tsv("sample/data/sample_file_ids_header.tsv")  
```

    ## Rows: 3029 Columns: 7

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): header_id
    ## dbl (6): DOTH41b_SF_C, DOTH42b_SF_C, DOTH43b_SF_C, DOTH44b_SF_C, DOTH45b_SF_...

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

If we check the first 6 rows, we see that the first column contain the
protein headers names from the fasta file. We can use that to know if it
comes from bacteria or human.

``` r
head(proteins_w_headers)
```

    ## # A tibble: 6 x 7
    ##   header_id     DOTH41b_SF_C DOTH42b_SF_C DOTH43b_SF_C DOTH44b_SF_C DOTH45b_SF_C
    ##   <chr>                <dbl>        <dbl>        <dbl>        <dbl>        <dbl>
    ## 1 >sp|A0A075B6~       36138.        42161           NA       214431       563395
    ## 2 >sp|A0A075B6~    13406900      14760200      2048760           NA       864503
    ## 3 >sp|A0A075B6~          NA            NA           NA           NA       165199
    ## 4 >tr|A0A075B7~     1434720        823523           NA       608046       240933
    ## 5 >tr|A0A075B7~    40787000      22400700     12236400      2292910      5819930
    ## 6 >sp|A0A087WS~     2666320       2868430           NA      3204700      3568570
    ## # ... with 1 more variable: DOTH46b_SF_C <dbl>

Now we annotate the proteins as bacterial or human

``` r
proteins_w_annotation <- mutate(proteins_w_headers, # input object
                                feature_category = if_else(str_detect(header_id,"_HUMAN"), 
                                                           true = "human",
                                                           false = "bacteria")) # feature category is the name that this column should have in the `annotation_features.txt`
```

If we look at the last two columns, we see that, for each protein, we
have a feature annotation.

``` r
proteins_w_annotation[24:33,c(5:8)]
```

    ## # A tibble: 10 x 4
    ##    DOTH44b_SF_C DOTH45b_SF_C DOTH46b_SF_C feature_category
    ##           <dbl>        <dbl>        <dbl> <chr>           
    ##  1      3967610           NA           NA human           
    ##  2           NA           NA           NA human           
    ##  3      6140610      6053250      7989480 human           
    ##  4           NA           NA           NA bacteria        
    ##  5      2073590           NA           NA bacteria        
    ##  6           NA           NA           NA bacteria        
    ##  7       231245           NA           NA bacteria        
    ##  8           NA           NA           NA bacteria        
    ##  9      1639280           NA           NA bacteria        
    ## 10      7970010           NA           NA bacteria

We can now select the two interesting columns (ID and feature_category)
and save them in a file.

``` r
annotation_features2 <- dplyr::select(proteins_w_annotation,
                                     ID = header_id, # changed index to ID column name
                                     feature_category)  # changed semi_type to feature_category column name
```

``` r
write_delim(annotation_features2,
            file = here("sample/data/annotation_features_metaprot.txt"), 
            delim = "\t")
```
