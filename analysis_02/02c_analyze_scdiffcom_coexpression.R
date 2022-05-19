library(data.table)

library(gridExtra)

library(purrr)
library(furrr)
library(future)
future::plan(future::multisession, workers = 8)

library(glue)
library(ggplot2)
library(ggpubr)

dir_path = "data/datasets_detections"

dt_cci = data.table()
for (file in list.files(dir_path)) {
  
  dataset = strsplit(file, ".csv")[[1]][1]
  dataset = strsplit(dataset, "detected_")[[1]][2]
  filepath = glue("{dir_path}/{file}")
  
  dt = data.table(read.csv(filepath))
  dt$dataset = dataset
  dt_cci = rbind(dt_cci, dt)
}

# Global coexpression
ggplot(data=dt_cci, aes(x=dataset, y=coexpression)) +
  geom_boxplot()

# Per dataset
ggplot(data=dt_cci, aes(y=dataset, x=coexpression, fill=REGULATION)) +
  geom_boxplot()
