library(tidyverse)
paths <- list.files("DATA/ISR-Yat",full.names = TRUE)

as.list(paths) %>% 
purrr::map_df(function(x){
    df <- read.table(x, sep = ",", header=T)
    # df[df == -9999]<- NA
    if(!any(colnames(df) == "GPP")){df <- df %>% mutate(GPP = c(0.1))}
    #   mutate(TIMESTAMP_START = strptime(df$TIMESTAMP_START, "%Y%m%d%H%M"),
    #          TIMESTAMP_END = strptime(df$TIMESTAMP_END, "%Y%m%d%H%M"))
df
  }
)->df

write_csv(df,file = "DATA/EC/IL-Yat_2000_2018.csv")
