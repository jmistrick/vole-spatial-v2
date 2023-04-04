#Functions from supplementary materials of Wanelik & Farine 2022 "A new method for characterising shared space use..."

#-------------------

## Function for calculating probability of detecting an individual at a distance, x, from its centroid
## Assuming a negative sigmoidal curve (defined by parameters a & b) describes the change in probability 
## of detection with increasing distance from the centroid (hereafter referred to as the 'home range')

prob_all <- function(x, a, b) {
  return(1/(1 + exp(-(a + b*x))))
}


#-------------------

## Function for generating network based on overlapping home ranges, from known individual centroids (xs and ys) and home range parameters (as & bs)
## Inputs: 	xs = vector of x location
##		ys = vector of y location 
## 		as = vector of a parameter values 
##		bs = vector of b parameter values
##    n = buffer area beyond home ranges
##    res = number of points along both x and y axes (higher = slower)
## Output: 	1. Matrix of overlapping home ranges 

get_network_2D <- function(xs, ys, as, bs, n=2, res=10) {
  
  edges <- matrix(NA,length(xs),length(xs))
  for (i in 1:(nrow(edges)-1)) {
    grid <- expand.grid(x=seq(min(xs)-n,max(xs)+n,length.out=res), ys=seq(min(ys)-n,max(ys)+n,length.out=res))
    d1 <- as.matrix(dist(cbind(c(xs[i],grid$x),c(ys[i],grid$y))))[1,-1]
    probs.1 <- prob_all(log(d1+1), as[i], bs[i])
    
    for (j in (i+1):nrow(edges)) {
      
      d2 <- as.matrix(dist(cbind(c(xs[j],grid$x),c(ys[j],grid$y))))[1,-1]
      probs.2 <- prob_all(log(d2+1), as[j], bs[j])
      
      top <- apply(cbind(probs.1,probs.2),1,min)
      bottom <- apply(cbind(probs.1,probs.2),1,max)
      
      edges[i,j] <- sum(top)/sum(bottom)
      edges[j,i] <- edges[i,j]
    }
  }
  return(edges)
}

#-------------------