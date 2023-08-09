# bx is the naive estimate, se_x is the corresponding standard error, k is the parameter in Projack
# k is recommended between 5-10. num_sim is the number of the simulations
Projack <-function(bx,se_x,k,num_sim){
  # z stands for the z-stats
  z=bx/se_x
  # obtain the number of  SNPs
  M=length(z)
  # d is the matrix for keeping the ordered z-values from the smallest to largest
  # for num_sim simulations
  d=matrix(NA,num_sim,M)
  
  # z1 is the new vector for z (serving for training part). 
  # c is to serve the role of held-out part.
  # as noted in paper, the re-ordering vectors in d matrix is based on the order of z1.
 
  for (i in 1:num_sim){
    gamma=rnorm(M,0,sd=sqrt(1/(k-1)))
    z1=z+gamma
    c=k*z-(k-1)*z1
    d[i,]=c[order(z1)]
  }
  
  # calculate the average over simulations
  c.z <- colMeans(d)
  
  #get the rank of the original z-stats
  rank_z <- rank(z)
  
  #multiply the standard error to get the estimate and return the corrected bx with the original order
  c.bx<-se_x*c.z[rank_z]
  
  #Ensure that the corrected estimate is less than the original estimate in magnitude
  ind1 <- which(abs(c.bx) > abs(bx))
  if(length(ind1) >0){c.bx[ind1] <- bx[ind1]}
  
  #Return the corrected estimate
  return(c.bx)
}  
