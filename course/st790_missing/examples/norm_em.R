#  Use the "norm" package downloaded from CRAN to
#  carry out the EM algorithm on the bivariate normal
#  example -- to install the package use
#
#  install.packages("norm") and choose a mirror site
#
#  Visit the CRAN website for the documentation

#  Load the libary for the norm package

library(norm)

#  Read in the data set -- the missing values are denoted
#  by the usual R convention NA
#  The package wants the data in a matrix, so we read it in
#  directly to a matrix

datamat <- matrix(scan("bvnormal.R.dat"),ncol=2,byrow=TRUE)

#  Run prelim.norm() to set up the data for the algorithm

prelim.em <- prelim.norm(datamat)

#  Run the EM algorithm initially with the default settings to 
#  get the form of the starting values

theta.init <- em.norm(prelim.em,showits=FALSE)

#  Extract the parameters so you can see the format

theta.init <- getparam.norm(prelim.em,theta.init,corr=TRUE)

#  Now you can replace the values by starting values you want
#  for a final run of the EM -- look at the elements of theta.init

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix

#  Use the sample complete case means as start value for mu

theta.init$mu <- c(5.042471,7.952805)

#  Use the sample complete case standard deviations

theta.init$sdv <- sqrt(c(1.026614,1.033791))

#  Use the sample complete case correlation

theta.init$r[1,2] <- 0.494137/sqrt(1.026614*1.033791)

#  Convert these to the form used by the function em.norm

theta.init <- makeparam.norm(prelim.em,theta.init)

#  Call EM algorithm again with these starting values and 
#  convergence criteria, etc, that you specify

theta.final <-
    em.norm(prelim.em,start=theta.init,showits=TRUE,maxits=100,criterion=1e-5)        
theta.final <- getparam.norm(prelim.em,theta.final,corr=TRUE)

#  Look at the final results 

mu.final <- theta.final$mu
sd.final <- theta.final$sdv
var.final <- sd.final^2
covar.final <- theta.final$r[1,2]*sd.final[1]*sd.final[2]

#  Results in same order as our homegrown EM program

theta <- round(c(mu.final,var.final[1],covar.final,var.final[2]),6)

#  Final values from running the above code, compare to our program
#  and SAS PROC MI results
                                        
#   Iterations of EM: 
#   1...2...3...4...5...6...7...8...
#   > theta
#  [1] 4.991665 7.926339 1.016241 0.514838 1.098030
