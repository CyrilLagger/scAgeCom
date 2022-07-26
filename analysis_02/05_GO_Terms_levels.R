library(ggridges)

r = readRDS("data/scAgeComShiny_data_14_05_2022.rds")

# hist(
#   r$ORA_table[
#     ORA_CATEGORY == "GO_TERMS" & ORA_REGULATION %in% c("UP", "DOWN") &
#       Dataset == "TMS FACS (male)" & Tissue == "Liver"
#   ]$`GO Level`, breaks = 100)

dt = r$ORA_table[
  ORA_CATEGORY == "GO_TERMS" & ORA_REGULATION %in% c("UP", "DOWN")
]

# dt = dt[
#   Dataset == "TMS FACS (female)"
# ]

dt$ID = paste0(dt$Dataset, "-", dt$Tissue)

ggplot(dt, aes(x = `GO Level`, y = ID, fill = Dataset)) +
  geom_density_ridges() +
  ggtitle("Distribution of levels for significant UP and DOWN-represented GO Terms across datasets") +
  theme_ridges() + 
  theme(legend.position = "none")
