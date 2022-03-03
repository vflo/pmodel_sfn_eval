library(sapfluxnetr)
library(tidyverse)
library(lubridate)
library(tsibble)
# path <- "../0.1.5/RData/plant"
# sfn_meta <- sapfluxnetr::read_sfn_metadata(folder = path, .write_cache = TRUE)
# save(sfn_meta, file = "sfn_meta.RData")
# load(file = "sfn_meta.RData")

# 
# FRA_FON <- read_sfn_data("FRA_FON", folder = path) %>%
#   sfn_metrics(period = "1 hour", .funs = list(~ mean(., na.rm = TRUE)),
#               solar = TRUE, interval = "general") %>%
#   sapfluxnetr::metrics_tidyfier(metadata = sfn_meta,interval = "general")
  
# 
# save(FRA_FON, file = "FRA_FON.RData")
load(file = "FRA_FON.RData")
