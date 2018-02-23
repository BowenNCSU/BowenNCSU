############################################################
#
#  Chapter 9, Argartroban pharmacokinetic data             
#                                                          
#  Fit nonlinear mixed effects model with within-individual         
#  power-of-mean variance function using the "pooled" 
#  GLS algorithm to obtain individual-specific parameter
#  estimates with pooled estimation of within-individual
#  variance parameters.  Then form the "response," "X,"
#  and "Z" matrices for the "linear mixed effects model"
#  and output to a file to be read into SAS PROC MIXED
#  or R lme() to estimate the population parameters beta and D
#
############################################################

#  function to print out matrices pretty

write.matrix <- function(x,file="",sep=" "){
  x <- as.matrix(x)
  p <- ncol(x)
  cat(dimnames(x)[[2]],format(t(x)),file=file,sep=c(rep(sep,p-1),
    "\n"),append=T)
}
	
#  put the individual mean function you want here with gradient matrix

meanfunc <- function(x,b1,b2,dose){

	tinf <- 240
	cl <- exp(b1)
	v <- exp(b2)
        t1 <- x<=tinf
        t2 <- tinf*(1-t1)+t1*x
        f1 <- (dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(x-tinf)/v)

#  compute analytical dervivatives -- gradient matrix X

        t3 <- (1-t1)*(x-tinf)
	temp1 <- (dose/cl)*exp(-cl*t3/v)
        temp2 <- (dose/(v^2))*exp(-cl*t3/v)
        
   meangrad <- array(0,c(length(x),2),list(NULL,c("b1","b2")))
   meangrad[,"b1"] <- temp1*(-(1-exp(-cl*t2/v))*(1/cl + t3/v)+(t2/v)*exp(-cl*t2/v))*cl 
   meangrad[,"b2"] <- temp2*((1-exp(-cl*t2/v))*t3-exp(-cl*t2/v)*t2)*v
 
   attr(f1,"gradient") <- meangrad
   f1

}

#  same mean function with no gradient

meanfunc2 <- function(x,b1,b2,dose){

	tinf <- 240
	cl <- exp(b1)
	v <- exp(b2)
        t1 <- x<=tinf
        t2 <- tinf*(1-t1)+t1*x
        f1 <- (dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(x-tinf)/v)
        f1
}

#  mean function with weights 

weightfunc <- function(x,b1,b2,dose,mut){
	
   f1 <- meanfunc2(x,b1,b2,dose)
   weightf <- f1/mut
   tinf <- 240
   cl <- exp(b1)
   v <- exp(b2)
   t1 <- x<=tinf
   t2 <- tinf*(1-t1)+t1*x
   t3 <- (1-t1)*(x-tinf)
   temp1 <- (dose/cl)*exp(-cl*t3/v)
   temp2 <- (dose/(v^2))*exp(-cl*t3/v)

#  compute analytical dervivatives -- create the gradient matrix X
        
   weightgrad <- array(0,c(length(x),2),list(NULL,c("b1","b2")))
   weightgrad[,"b1"] <- temp1*(-(1-exp(-cl*t2/v))*(1/cl + t3/v)+(t2/v)*exp(-cl*t2/v))*cl/mut
   weightgrad[,"b2"] <- temp2*((1-exp(-cl*t2/v))*t3-exp(-cl*t2/v)*t2)*v/mut
   attr(weightf,"gradient") <- weightgrad
   weightf
}
	

#  set start values, etc

# max number of iterations for fitting algorithm

cmax <- 20

#  start values

bstart <- list(b1=-6.0,b2=-2.0)

delstart <- list(delta=0.5)
p <- 2

#  name of output file

outfile <- "arg_pooledgls.Rout"
cat("ARGATROBAN DATA -- GLS ALGORITHM",file=outfile,"\n","\n","\n",append=F)

#  read in data

thedat <- read.table("argconc.dat")
colnames(thedat) <- c("obsno","id","dose","time","conc")

indiv <- thedat$id   # individual indicator
ds <- thedat$dose    # dose
xs <- thedat$time     # time
ys <- thedat$conc     # concentration

n <- length(xs)
uindiv <- unique(indiv)
m <- length(uindiv)
		
# form individual 2nd stage covariate design matrices Ai
# these are all identity matrices in this case

aimat <- NULL
r <- p

a1 <- diag(p)
	
for(i in 1:m){

      aimat <- rbind(aimat,a1)
	
}

#  "mean" function for the for use in "trick" to solve pooled
#  quadratic estimating equation - see below

pltrkfunc <- function(resid,mudot,mu,delta){
   trk <- resid*((mudot/mu)**delta)
   trkgrad <- array(0,c(length(mu),1),list(NULL,c("delta")))
   trkgrad[,"delta"] <- trk*log(mudot/mu)
   attr(trk,"gradient") <- trkgrad   #  analytic derivative
   trk
}

#  Stage 1 -- Pooled GLS estimation

#  Step 1 -- initial fit by OLS for each indiv

cat("INITIAL OLS ESTIMATION",file=outfile,"\n","\n",append=T)

bolsmat <- NULL     #  matrix to contain OLS estimates
muvec <- NULL       #  vector to contain all predicted values
residvec <- NULL    #  vector to contain all residuals
		
for (i in 1:m){

id <- uindiv[i]
y <- ys[indiv==id]
x <- xs[indiv==id]
d <- ds[indiv==id]	

olsdat <- data.frame(x,y,d)

ols.fit <- nls(y ~ meanfunc(x,b1,b2,d),olsdat,bstart)
bols <- coef(ols.fit)
	
mu <- meanfunc2(x,bols[1],bols[2],d)
resid <- y-mu

bols <- matrix(bols,ncol=p,byrow=T)	
bolsmat <- rbind(bolsmat,bols)

muvec <- c(muvec,mu)
residvec <- c(residvec,resid)	

}

bols <- cbind(uindiv,bolsmat)

muvec <- matrix(muvec,ncol=1)
residvec <- matrix(residvec,ncol=1)
		
cat("OLS estimates",file=outfile,"\n","\n",append=T)

write.matrix(round(bols,6),file=outfile)
	
#  pooled OLS estimate of sigma

sigma <- sqrt(sum(residvec*residvec)/(n-m*p))

cat(file=outfile,"\n","\n",append=T)	
cat("OLS pooled estimate of sigma ",sigma,"\n","\n",file=outfile,append=T)
			
# begin iteration loop for GLS-PL
				
cat("GLS POOLED ESTIMATION",file=outfile,"\n","\n",append=T)
	
for (k in 1:cmax){

#  compute the geometric mean and predicted values
	
#mudot <- prod(muvec)**(1/n)
mudot <- exp((1/n)*sum(log(muvec)))
mudot <- rep(mudot,length(muvec))
dummy <- rep(0,length(muvec))
pldat <- data.frame(residvec,muvec,mudot)

#  Step 2 -- estimate delta using the computational "trick" to solve 
#  the pooled quadratic estimating equation described in Sections
#  2.4.2 and 5.2.2. of Davidian and Giltinan (1995) by profiling out 
#  sigma^2 from the objective function and then by algebra casting
#  estimation of delta as a "nonlinear regression" problem with "dummy"
#  responses all equal to zero.  See this reference for details.  
	
pl.fit <- nls(dummy ~ pltrkfunc(residvec,mudot,muvec,delta),pldat,
	delstart,control=list(maxiter=500))
delta <- coef(pl.fit)
	
cat("Iteration ",k,"\n",file=outfile,append=T)
cat("Estimate of delta ",delta,"\n",file=outfile,append=T)	

#  Step 3 -- update estimates of beta for each individual

bglsmat <- NULL
new.muvec <- NULL
new.residvec <- NULL

for (i in 1:m){

id <- uindiv[i]
y <- ys[indiv==id]
x <- xs[indiv==id]
d <- ds[indiv==id]	

mu <- muvec[indiv==id]
	
mut <- mu^delta
ymut <- y/mut
	
glsdat <- data.frame(x,ymut,d)

glspl.fit <- nls(ymut ~ weightfunc(x,b1,b2,d,mut),glsdat,
	bstart,control=list(maxiter=100))

bgls <- coef(glspl.fit)

mu <- meanfunc2(x,bgls[1],bgls[2],d)
resid <- y-mu

bgls <- matrix(bgls,ncol=p,byrow=T)	
bglsmat <- rbind(bglsmat,bgls)

new.muvec <- c(new.muvec,mu)
new.residvec <- c(new.residvec,resid)	

}

bgls <- cbind(uindiv,bglsmat)
muvec <- matrix(new.muvec,ncol=1)
residvec <- matrix(new.residvec,ncol=1)
		
}

#  Compute final estimate of sigma^2

g <- muvec**delta

sigma2 <- sum((residvec/g)**2)/(n-m*p)
sigma <- sqrt(sigma2)

cat(file=outfile,"\n","\n","\n",append=T)
cat("Final GLS estimates after ",cmax," iterations",
	file=outfile,"\n","\n",append=T)

write.matrix(round(bgls,6),file=outfile)

cat(file=outfile,"\n","\n",append=T)

cat("Final pooled estimate of sigma ",sigma,"\n","\n",file=outfile,append=T)
	
#  now transform problem to allow use of mix model software to fit 2nd stage

mixdat <- NULL
	
for (i in 1:m){

id <- uindiv[i]
y <- ys[indiv==id]
x <- xs[indiv==id]
d <- ds[indiv==id]	

thisid <- rep(i,p)

#  get covariance matrix for bglsi

bglsi <- bgls[i,2:(p+1)]
mui <- meanfunc2(x,bglsi[1],bglsi[2],d)
muit <- mui^theta

xmati <- attr(weightfunc(x,bglsi[1],bglsi[2],d,muit),"gradient")

cmati <- round(sigma^2*solve(t(xmati) %*% xmati),6)	

#  cholesky decomp of inverse 

cinvhalf <- chol(round(solve(cmati),6))

#  Create the "data" for stage 2 -- response, "Xi" and "Zi"

ai <- aimat[((i-1)*p+1):(i*p),]

respi <- cinvhalf%*%bglsi

thisxi <- cinvhalf%*%ai
	
thiszi <- cinvhalf
	
thisdat <- cbind(thisid,respi,as.matrix(thisxi),as.matrix(thiszi))

mixdat <- rbind(mixdat,thisdat)

}

#  write the data to a file to read into PROC MIXED or lme()

write.table(mixdat,file="argmixsig.dat",col.names=FALSE,row.names=FALSE,append=FALSE)



