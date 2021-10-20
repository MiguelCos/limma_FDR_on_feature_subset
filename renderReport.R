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
#'## Number of features which are differentially expressed: `r n_significant`
#'
#'
#'
#+ fig.width=11.7, fig.height= 7.1, echo = FALSE, message = FALSE, warning = FALSE
plotly::ggplotly(volcanoes)
#'
#'
#'
#'
#'
print(sessionInfo())
