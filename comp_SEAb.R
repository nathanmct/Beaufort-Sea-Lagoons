# Function to compare all the columns in the matrix mat
# the names of the columns are in s
#
comp_SEAb <- function(mat,s) {
  
  #s <- unique(s)
  s <- as.character(levels(s))
  comp1 <- do.call(rbind, apply(combn(seq_along(s), 2), 2,
                                function(x) data.frame(g1=s[x[1]], g2=s[x[2]],
                                                       prop=sum(mat[,x[1]] > mat[,x[2]])/nrow(mat))))
  
  require(pander)
  knitr::kable(comp1)
}