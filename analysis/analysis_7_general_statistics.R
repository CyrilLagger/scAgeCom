####################################################
##
## Project: scAgeCom
##
## Last update - June 2021
##
## cyril.lagger@liverpool.ac.uk
## ursu_eugen@hotmail.com
## anais.equey@etu.univ-amu.fr
##
## general statistics
##
####################################################
##

## libraries ####
library(data.table)
library(ggplot2)

## load results ####

cci_stat_table <- readRDS("../data_scAgeCom/analysis/outputs_data/data_4_CCI_table_unprocessed.rds")
shiny_results <- readRDS("../data_scAgeCom/analysis/outputs_data/scAgeCom_shiny_data.rds")


## create regulation table for supplemental data and stat ####

regulation_table <- cci_stat_table[
  ,
  .N,
  by = c("IS_CCI_EXPRESSED_YOUNG", "IS_CCI_SPECIFIC_YOUNG", "IS_CCI_SCORE_YOUNG",
         "IS_CCI_EXPRESSED_OLD", "IS_CCI_SPECIFIC_OLD", "IS_CCI_SCORE_OLD",
         "IS_DE_LOGFC", "IS_DE_SIGNIFICANT", "DE_DIRECTION",
         "IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD", "IS_CCI_DE",
         "REGULATION")
]

## summary statistics ####

cci_stat_table[, .N]
cci_stat_table[, .N, by = "REGULATION"][, N/sum(N)*100]

regulation_table[
  REGULATION == "NSC",
  sum(N),
  by = c("IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD")
][, V1/sum(V1)*100]

regulation_table[
  REGULATION %in% c("UP", "DOWN"),
  sum(N),
  by = c("IS_CCI_DETECTED_YOUNG", "IS_CCI_DETECTED_OLD")
]

## number of LRI per ER ####

cci_stat_table[, DTER_CELLTYPES := paste(Dataset, Tissue, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, sep = "_")]
NLRI_table <- cci_stat_table[, .N, by = c("Dataset", "DTER_CELLTYPES")]

ggplot(
  data = NLRI_table
) + geom_histogram(
  aes(
    x = N,
    fill = Dataset,
    color = Dataset
  ),
  bins = 100,
  position = "identity",
  alpha = 0.3
)


distr_LRI <- c(NLRI_table$N, rep(0, sum(cci_stat_table[, uniqueN(EMITTER_CELLTYPE), by = c("Dataset", "Tissue")][, V1^2]) - nrow(NLRI_table) ))
mean(distr_LRI)
sd(distr_LRI)
hist(distr_LRI, breaks = 100)

## regulation distribution per tissue ####

regulation_distr <- cci_stat_table[, .N, by = c("Dataset", "Tissue", "REGULATION")]
regulation_distr <- dcast.data.table(
  regulation_distr,
  Dataset + Tissue ~ REGULATION,
  value.var = "N"
)
regulation_distr[is.na(regulation_distr)] <- 0

regulation_distr[
  ,
  total := UP + DOWN + FLAT + NSC
]
regulation_distr[
  ,
  c(
    "UP",
    "DOWN",
    "FLAT",
    "NSC"
  ) :=
    list(
      UP/total,
      DOWN/total,
      FLAT/total,
      NSC/total
    )
]
regulation_distr_long <- melt.data.table(
  regulation_distr,
  id.vars = c("Dataset", "Tissue"),
  measure.vars = c("UP", "DOWN", "FLAT", "NSC"),
  variable.name = "REGULATION",
  value.name = "pct"
)
regulation_distr_long[, id := paste(Dataset, Tissue, sep = "_")]

regulation_distr_long[
  ,
  Dataset := factor(
    Dataset,
    c("TMS FACS (male)", "TMS FACS (female)", "TMS Droplet (male)", "TMS Droplet (female)", "Calico Droplet (male)")
  )
]

## Figure 5 ####
library(webshot)
library(htmlwidgets)
library(shiny)

display_KEYWORD_counts_local <- function(
  ora_keyword_counts,
  category,
  regulation,
  go_aspect = NULL
) {
  ORA_CATEGORY <- ORA_REGULATION <- ASPECT <- NULL
  dt <- ora_keyword_counts[
    ORA_CATEGORY == category &
      ORA_REGULATION == regulation
  ]
  data.table::setnames(dt, old = "VALUE", new = category)
  if (category == "GO Term") {
    temp_aspect <- ifelse(
      go_aspect == "Biological Process",
      "biological_process",
      ifelse(
        go_aspect == "Molecular Function",
        "molecular_function",
        "cellular_component"
      )
    )
    dt <- dt[ASPECT == temp_aspect]
    dt <- dt[, -c(1,2,10)]
    if (go_aspect == "Biological Process") {
      category_label <- paste0("GO ", go_aspect, "es")
    } else {
      category_label <- paste0("GO ", go_aspect, "s")
    }
  } else {
    dt <- dt[, -c(1,2,10,11)]
    if(grepl("Family", category)){
      category_label <- sub("Family", "Families", category)
    } else {
      category_label <- paste0(category, "s")
    }
  }
  DT <- DT::datatable(
    data = dt,
    filter = list(
      position ="top",
      clear = FALSE,
      plain = FALSE
    ),
    class = "display compact",
    options =list(
      pageLength = 5,
      dom = '<"top"f>rt<"bottom"lip><"clear">',
      columnDefs = list(
        list(width = '300px', targets = c(1))
      )
    ),
    caption = tags$caption(
      style = paste0(
        "caption-side: top; ",
        "text-align: center; ",
        "color: black; ",
        "font-size: 120%;"
      ),
      paste0(
        "Number of tissues in which ",
        category_label,
        " are over-represented among ",
        regulation,
        "-regulated cell-cell interactions"
      )
    )
  ) %>% DT::formatStyle(
    colnames(dt)[-1],
    `text-align` = 'center'
  )
  if (category == "GO Term") {
    DT <- DT %>% DT::formatStyle(c(7), `border-right` = "solid 2px")
  }
  DT
}

table_LRI_up <- display_KEYWORD_counts_local(
  shiny_results$ORA_KEYWORD_COUNTS,
  category = "Ligand-Receptor Interaction",
  regulation = "UP"
)
table_LRI_up$width <- "600px"

table_LRI_up

table_LRI_up_html <- "table_LRI_up.html"
saveWidget(table_LRI_up, table_LRI_up_html)
webshot(
  table_LRI_up_html,
  "../../../../../table_LRI_up.png",
  zoom = 3
) 

table_LRI_down <- display_KEYWORD_counts_local(
  shiny_results$ORA_KEYWORD_COUNTS,
  category = "Ligand-Receptor Interaction",
  regulation = "DOWN"
)
table_LRI_down$width <- "600px"

table_LRI_down_html <- "table_LRI_down.html"
saveWidget(table_LRI_down, table_LRI_down_html)
webshot(
  table_LRI_down_html,
  "../../../../../table_LRI_down.png",
  zoom = 3
)

table_L_down <- display_KEYWORD_counts_local(
  shiny_results$ORA_KEYWORD_COUNTS,
  category = "Receptor",
  regulation = "DOWN"
)
table_L_down$width <- "600px"

table_LRI_down_html <- "table_LRI_down.html"
saveWidget(table_LRI_down, table_LRI_down_html)
webshot(
  table_LRI_down_html,
  "../../../../../table_LRI_down.png",
  zoom = 3
)


table_BP_up <- display_KEYWORD_counts_local(
  shiny_results$ORA_KEYWORD_COUNTS,
  category = "GO Term",
  regulation = "UP",
  go_aspect = "Biological Process"
)
table_BP_up$width <- "650px"

table_BP_up_html <- "table_BP_up.html"
saveWidget(table_BP_up, table_BP_up_html)
webshot(
  table_BP_up_html,
  "../../../../../table_BP_up.png",
  zoom = 3
)


plot_KEYWORD_summary_local <- function(
  ora_keyword_summary,
  ora_keyword_template,
  category,
  keyword
) {
  ORA_CATEGORY <- VALUE <- Regulation <- i.ORA_REGULATION <-
    Dataset <- Tissue <- NULL
  dt <- copy(ora_keyword_summary)[
    ORA_CATEGORY == category
  ]
  if (!(keyword %in% dt$VALUE)) {
    stop("`keyword` not found")
  }
  dt <- dt[
    VALUE == keyword
  ]
  dt <- copy(ora_keyword_template)[
    dt,
    on = c("Tissue", "Dataset"),
    Regulation := i.ORA_REGULATION
  ]
  dt[is.na(dt)] <- "Not Detected"
  p <- ggplot2::ggplot(dt) +
    ggplot2::geom_tile(
      ggplot2::aes(
        Dataset,
        Tissue,
        fill = Regulation,
        width = 0.9,
        height = 0.9
      ),
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c(
        #"No Data" = "transparent",
        "Not Over-represented" = "white",
        "Not Detected" = "gray",
        "UP" = "red",
        "DOWN" = "blue",
        "FLAT" = "green"#,
        #"UP:DOWN" = "yellow"
      )
    ) +
    ggplot2::ggtitle(
      stringr::str_trunc(
        paste0(
          "Over-representation of ",
          keyword
        ),
        70, 
        "right"
      )
    ) +
    ggplot2::scale_x_discrete(
      limits = c(
        "TMS FACS (male)",
        "TMS FACS (female)" ,
        "TMS Droplet (male)",
        "TMS Droplet (female)",
        "Calico Droplet (male)"
      ),
      labels = c(
        "TMS\nFACS\n(male)",
        "TMS\nFACS\n(female)",
        "TMS\nDroplet\n(male)",
        "TMS\nDroplet\n(female)",
        "Calico\nDroplet\n(male)"
      )
    ) +
    ggplot2::scale_y_discrete(
      limits = sort(
        unique(dt$Tissue),
        decreasing = TRUE
      )
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    #theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(text=ggplot2::element_text(size = 30)) +
    ggplot2::theme(
      axis.text=ggplot2::element_text(size = 30),
      legend.position = c(0.85, 0.8)
      )# +
  #ggplot2::theme(legend.position = c(0.8, 0.8))
  p
  #plotly::ggplotly(
  #  p,
  #  source = "TCA_PLOT_KEYWORD_SUMMARY",
  #  tooltip = c("Dataset", "Tissue", "Regulation")
  #) #%>% plotly::layout(
  #  legend = list(
  #    title = list(text = "")
  #  ) 
  # )
}

plot_KEYWORD_summary_local(
  ora_keyword_summary = shiny_results$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_results$ORA_KEYWORD_TEMPLATE,
  category = "GO Term",
  keyword = "T cell differentiation"
)
#2000x1400

plot_KEYWORD_summary_local(
  ora_keyword_summary = shiny_results$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_results$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "B2m:Cd3g"
)

plot_KEYWORD_summary_local(
  ora_keyword_summary = shiny_results$ORA_KEYWORD_SUMMARY,
  ora_keyword_template = shiny_results$ORA_KEYWORD_TEMPLATE,
  category = "Ligand-Receptor Interaction",
  keyword = "Gpi1:Amfr"
)

## Figure 6 ####

ggplot(
  data = regulation_distr_long,
  aes(
    y = Tissue,
    x = pct,
    fill = REGULATION
  )
) + geom_bar(
  stat = "identity",
  position = "fill",
  alpha = 0.8
) + scale_fill_manual(
  "Age-Regulation",
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "green", "NSC" = "grey")
) + facet_wrap(
  ~ Dataset,
  ncol = 5
) + ggplot2::scale_y_discrete(
    limits = sort(
      unique(regulation_distr_long$Tissue),
      decreasing = TRUE
    )
) + xlab(
  "Fraction of CCIs per regulation group"
) + theme(
  text = element_text(size = 40, face = "bold"),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 36, face = "bold"),
  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.905, 0.8)
)
#manual save 3000x1800






