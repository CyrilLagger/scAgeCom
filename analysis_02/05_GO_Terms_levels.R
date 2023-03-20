library(data.table)
library(ggplot2)
library(ggridges)

r = readRDS("data/scAgeComShiny_data_14_05_2022.rds")

# hist(
#   r$ORA_table[
#     ORA_CATEGORY == "GO_TERMS" & ORA_REGULATION %in% c("UP", "DOWN") &
#       Dataset == "TMS FACS (male)" & Tissue == "Liver"
#   ]$`GO Level`, breaks = 100)

dt = r$ORA_table[
  ORA_CATEGORY == "GO_TERMS" 
  & ORA_REGULATION %in% c("UP", "DOWN")
]

# dt = dt[
#   Dataset == "TMS FACS (female)"
# ]

dt$ID = paste0(dt$Dataset, "-", dt$Tissue)

# Ridge plot
ggplot(dt, aes(x = `GO Level`, y = ID, fill = Dataset)) +
  geom_density_ridges() +
  ggtitle("Distribution of levels for significant UP and \nDOWN-represented GO Terms across datasets") +
  theme_ridges() + 
  theme(legend.position = "none")

## Checking the average levels of overrepresented terms
dt_medians = dt[!is.na(`GO Level`), list(MEDIAN_GO_LEVEL = median(`GO Level`)), by=ID]
(ggplot(data=dt_medians, aes(x=MEDIAN_GO_LEVEL))
  + geom_histogram(binwidth = 0.1)
  + scale_x_continuous(breaks = seq(4, 8, 0.5))
  + xlab("Median GO Level")
  + ggtitle(glue::glue("Distribution of the median of the \nGO levels across all datasets.\nN={nrow(dt_medians)}\n"))
)
median(dt_medians$MEDIAN_GO_LEVEL)
