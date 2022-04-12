## SOIL WATER CONTENT AGGREGATION FUNCTION

swvl_aggregation <- function(depth,swvl1,swvl2,swvl3,swvl4){
  # check if there is swvl data
  if(all(is.na(swvl1))){stop("swvl1 data not provided")}
  if(all(is.na(swvl2))){stop("swvl2 data not provided")}
  if(all(is.na(swvl3))){stop("swvl3 data not provided")}
  if(all(is.na(swvl4))){stop("swvl4 data not provided")}
  
  #use method of aggregation depending on the presence of soil depth
  if(all(is.na(depth))){
    res <- swvl1*7/289 + swvl2*(28-7)/289 + swvl3*(100-28)/289 + swvl4*(289-100)/289
  }else{ #aggregation deppending on the depth
    depth <- depth %>% unique()
    if(depth>=289){res <- swvl1*7/289 + swvl2*(28-7)/289 + swvl3*(100-28)/289 + swvl4*(289-100)/289}
    if(depth<289& depth>100){res <- swvl1*7/depth + swvl2*(28-7)/depth + swvl3*(100-28)/depth + swvl4*(depth-100)/depth}
    if(depth<=100 & depth>27){res <- swvl1*7/depth + swvl2*(28-7)/depth + swvl3*(depth-28)/depth}
    if(depth<=27& depth>7){res <- swvl1*7/depth + swvl2*(depth-7)/depth}
    if(depth<=7){res <- swvl1}
    }
  
  return(res)

}

