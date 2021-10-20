#' ---
#' title: "Limma - Differential expression analyses; two-group comparison: script report"
#' format:
#'html:
#'                    theme: cosmo
#'toc: true
#'toc-depth: 4
#'number-sections: true
#' ---

#'## Dataset/experiment code: `r exper_code`
#'
#'## Grouping information summary: 
#+ echo = FALSE, message = FALSE
kable(group_by(annot_dat, Group = Group) %>% summarise(N = n()))%>%
          kable_styling(bootstrap_options = "striped", full_width = F, position = "left")
#'
#'
#'## Total number of features that are differentially abundant: `r n_significant_total`
#'
#'## Number of interesting/selected features that are differentially expressed: `r n_significant_sel_features`
#'
#'
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE, message = FALSE, warning = FALSE
plotly::ggplotly(volcanoes)
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE, message = FALSE, warning = FALSE
plotly::ggplotly(volcano_features)
#'
#'Note: p-values of interesting/selected features were independently corrected using BH method. 
#'
#'
#'
#'
print(sessionInfo())
