####################################################
##
## Project: scAgeCom
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## additional figure preparation
##
####################################################
##

## Extended Data LRI distribution (not included) ####

ExtFig1 <- ggplot(
  data = NLRI_template,
  aes(
    y = tissue,
    x = N
  )
) + geom_boxplot(
  outlier.shape = NA
) + facet_wrap(
  ~ dataset,
  ncol = 5
) + ggplot2::scale_y_discrete(
    limits = sort(
      unique(NLRI_table$tissue),
      decreasing = TRUE
    )
) + xlab(
  "Number of detected LRIs per cell-type pair"
) + geom_jitter(
  size = 0.2
) + theme(
  #text = element_text(size = 40, face = "bold"),
  #axis.text.x = element_text(size = 30),
  #axis.text.y = element_text(size = 36, face = "bold"),
  #axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  panel.spacing = unit(2, "lines")
)
ExtFig1
