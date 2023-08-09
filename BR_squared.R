##bx is the naive estimate and se_x is the corresponding standard error, 
##num_boot is the number of bootstraps you can specify (1000 or greater is recommended)
##nobs is the number of observations for the G-X association
BRss <- function (bx,se_x,num_boot,nobs){
  z=bx/se_x
  M <- length(bx)
  boot_diff <- rep(0,M)
  for(i in 1:num_boot){
    #dd is the matrix with 2 columns. 
    #The first column is for saving the within-sample bootstrap estimate
    #The second column is for saving the corrected out-of-sample bootstrap estimate
    dd = matrix(NA,M,2)
    #p_de is a vector about the correlation between within-sample bootstrap estimate
    #and out-of-sample bootstrap estimate
    p_de <- -(1+2*bx^2)*sqrt(exp(1)^(-1))/(1-bx^2)
    dd[,1]=rnorm(M,mean=bx,sd=sqrt((1-bx^2)/nobs))
    dd[,2]=rnorm(M,mean=bx,sd=sqrt((1-bx^2)/(exp(1)^(-1)*nobs)*(1-p_de^2)))
    chis <- rchisq(M,df=nobs-1)
    #Simulate the variance of the within-sample bootstrap estimate from a chi-squared distribution
    vd_star <- (1-bx^2)/nobs/(nobs-1) * chis
    #Obtain the z-statistics based on the within-sample estimates
    z_star <- dd[,1]/sqrt(vd_star)
    #Get the order of z-stats based on the within-sample estimates
    ind <- order(z_star)
    #Reorder the bootstrap estimates and also the original estimate based on z-stats from the with-sample estimate
    boot.beta.e=dd[,2][ind]
    beta.d=dd[,1][ind]
    beta.N=bx[ind]
    #Correct the out-of-sample bootstrap estimate by accounting for the negative correlation
    boot_diff <- boot_diff + (beta.d-boot.beta.e)
  }
  #Get the rank of naive (original) z stats
  rank_z <- rank(z)
  #Averaging over the bootstrap difference
  boot_diff <- boot_diff/ num_boot
  #correct the naive estimate by matching the rank of re-ordered z-stats
  c.bx = bx - boot_diff[rank_z]
  #Ensure that the corrected estimate is less than the original estimate in magnitude
  ind1 <- which(abs(c.bx) > abs(bx))
  if(length(ind1) >0){c.bx[ind1] <- bx[ind1]}
  #Return the corrected estimates for the G-X association
  c.bx
}
