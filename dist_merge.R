dist <- function(x1, y1, x2, y2) {
  ((x1-x2)^2 + (y1-y2)^2)^0.5
}

dist_merge <- function(x, y, xeast, xnorth, yeast, ynorth){
  tmp <- t(apply(x[,c(xeast, xnorth)], 1, function(x, y){
    dists <- apply(y, 1, function(x, y) dist(x[2], x[1], y[2], y[1]), x)
    cbind(1:nrow(y), dists)[dists == min(dists),,drop=F][1,]
    }
  , y[,c(yeast, ynorth)]))
  tmp <- cbind(x, min.dist=tmp[,2], y[tmp[,1],-match(c(yeast, ynorth), names(y))])
  row.names(tmp) <- NULL
  tmp
}

# modification of dist_merge to retain coordinates of y
dist_merge_mod <- function(x, y, xeast, xnorth, yeast, ynorth){
  tmp <- t(apply(x[,c(xeast, xnorth)], 1, function(x, y){
    dists <- apply(y, 1, function(x, y) dist(x[2], x[1], y[2], y[1]), x)
    cbind(1:nrow(y), dists)[dists == min(dists),,drop=F][1,]
  }
  , y[,c(yeast, ynorth)]))
  tmp <- cbind(x, 
               min.dist=tmp[,2], 
               y[tmp[,1],-match(c(yeast, ynorth), names(y))],
               y_long = (c(y[tmp[,1], match(c(yeast), names(y))]))[[1]],
               y_lat = (c(y[tmp[,1],match(c(ynorth), names(y))]))[[1]])
  row.names(tmp) <- NULL
  tmp
}