/**************************************************************************
*                                                                         *
*  CHAPTER 7 - Example of using the GLS iterative algorithm to fit        *
*  a nonlinear model with variance function                               *
*                                                                         *
*  Data on the pharmacokinetics of indomethacin on a single subject       *
*                                                                         *
*  Define a SAS macro to perform the GLS iterative algorithm with         *
*  variance function a power 2 delta of the mean                          *
*  Here, we fit the biexponential model for the mean E(Y_j | x_j) =       *
*                                                                         *
*     f(x,b) = exp(b1) exp(-exp(b2) x) + exp(b3) exp(-exp(b4) x),         *
*                                                                         *
*  which is a common PK model for a two-compartment representation of     *
*  the body, parameterized to ensure positivity of all the                *
*  estimates.  We fit the power of the mean variance function             *
*                                                                         *
*      var(Y_j | x_j) = sigma^2 f(x,b)^{2 delta}                          *
*                                                                         *
*  where delta is estimated by solving the quadratic estimating           *
*  equation in (7.21) of the notes.                                       *                       
*                                                                         *
*  SAS evidently does not have a procedure that allows the variance to    *
*  function to include an unknown parameter; thus, we implement the       *
*  algorithm manually. We start from the OLS estimate and form fixed      *
*  weights at each iteration, holding the current estimate of beta        *
*  and that of delta fixed to form the weights (variation on step 2)      *
*                                                                         *
**************************************************************************/

options ps=55 ls=80 nodate;

/**************************************************************************
*                                                                         *
*  Define a sas macro called "glsalg" to control the GLS iteration.       *
*  The macro takes as input the name of the SAS working data set, the     *
*  names of the x and y variables (here, x is a scalar, but the program   *
*  may be modified easily to pass more than one explanatory variable),    *
*  and starting values for each WLS calculation.                          *
*                                                                         *
**************************************************************************/

%macro glsalg(dset,xvar,yvar,b1,b2,b3,b4);

/**************************************************************************
*                                                                         *
*  Set up the data set for the first pass through the algorithm, which    *
*  will compute the OLS estimate.  The variable "pred1" containing the    *
*  predicted values from the previous iteration is set equal to 1 here,   *
*  so that the first pass with weights based on pred1^delta will actually *
*  be OLS (weights all = 1).                                              *
*                                                                         *
**************************************************************************/

data &dset; set &dset; 
  pred1=1;
  if &xvar=. or &yvar=. then delete;
  b1_new=&b1;
  b2_new=&b2;
  b3_new=&b3;
  b4_new=&b4;

/**************************************************************************
*                                                                         *
*  Set the initial value for delta (OLS)                                  *
*                                                                         *
**************************************************************************/

  delta1=0.0;  

/**************************************************************************
*                                                                         *
*  Call the macro "step3," defined below, that actually implements the    *
*  call to PROC NLIN to do WLS with the current set of fixed weights.     *
*  The argument sets options for supressing the printing of the results   *
*  of each call to PROC NLIN. We do 10 iterations and print out them      *
*  results for the last one                                               *
*                                                                         *
**************************************************************************/

%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(1)
%step3(0)

/**************************************************************************
*                                                                         *
*  After the final iteration, compute the final estimate of sigma         *
*  using the newest values of the weights                                 *
*                                                                         *
**************************************************************************/

data sigma; set outnlin(keep = resid pred);
data sigma2; set new(keep = delta1);
data sigma; merge sigma sigma2; 

data sigma; set sigma; wresid=resid/(pred**delta1);

proc means data=sigma noprint; 
  var wresid;
  output out=sigout uss=rss;

data sigout; set sigout; sigma=sqrt(rss/7);  * lazy 7 = n-p;

data sigout; set sigout(keep = sigma); if _n_>1 then delete;

proc print data=sigout; title2 "Final estimate of sigma";

%mend glsalg;

/**************************************************************************
*                                                                         *
*  Define the macro "step3."  The argument "foo" controls options in      *
*  calling PROC NLIN (see below).                                         *
*                                                                         *
**************************************************************************/

%macro step3(foo); 

/**************************************************************************
*                                                                         *
*  PROC NLIN prints out a summary by default.  We print out this full     *
*  output only for the final GLS fit; the "noprint" option is invoked     *
*  for all other fits.                                                    *
*                                                                         *
**************************************************************************/

%if &foo>0 %then %do;
  %let opt=noprint;
%end;
%else %do; 
  %let opt= ; 
%end;

/**************************************************************************
*                                                                         *
*  The call to PROC NLIN to implement the current WLS fit with delta      *
*  equal to delta1, the previous estimate                                 *
*                                                                         *
**************************************************************************/

proc nlin data=&dset method=gauss &opt;

  parms b1=&b1 b2=&b2 b3=&b3 b4=&b4;  
  if _iter_=-1 then do;
    b1=b1_new; b2=b2_new; b3=b3_new; b4=b4_new;
  end;

/**************************************************************************
*                                                                         *
*  Programming statements to define the mean function and its derivatives *
*                                                                         *
**************************************************************************/

  eb1=exp(b1); eb2=exp(b2); eb3=exp(b3); eb4=exp(b4);

  f=eb1*exp(-eb2*&xvar) + eb3*exp(-eb4*&xvar);

  db1=eb1*exp(-eb2*&xvar);
  db2=-eb1*eb2*&xvar*exp(-eb2*&xvar);
  db3=eb3*exp(-eb4*&xvar);
  db4=-eb3*eb4*&xvar*exp(-eb4*&xvar);

/**************************************************************************
*                                                                         *
*  The model statement and derivative statements to tell PROC NLIN the    *
*  form of the model and derivatives                                      *
*                                                                         *
**************************************************************************/

  model &yvar = f;
  der.b1=db1;
  der.b2=db2;
  der.b3=db3;
  der.b4=db4;
 
/**************************************************************************
*                                                                         *
*  The "_weight_" statement defines the values of the weights to use for  *
*  WLS.  Here, we use weights computed using the fixed value of delta1    *
*  and the values of the predicted values from the previous iteration.    *
*  Thus, the weights are FIXED (do not depend on current predicted        *
*  values).  PROC NLIN thus knows to do WLS with FIXED weights rather     *
*  than IRWLS.                                                            *
*                                                                         *
**************************************************************************/

  _weight_ = (1/abs(pred1))**(2*delta1);

/**************************************************************************
*                                                                         *
*  Form an output data set containing everything in the input data set    *
*  plus the predicted values (to set up weights for next time), the       *
*  residuals (so we can calculate sigma on the final iteration), and the  *
*  parameter estimates from this iteration (so we can print them).        *
*                                                                         *
**************************************************************************/

  output out=outnlin p=pred r=resid sse=sigma parms=b1 b2 b3 b4;
run;

/**************************************************************************
*                                                                         *
*  Now estimate delta by solving the quadratic estimating equation        *
*  (7.21).  We do this by using a computational "trick" described in      *
*  Section 2.4.2 of Davidian and Giltinan (1995) by profiling out         *
*  sigma^2 from the objective function and then by algebra casting        *
*  estimation of delta as a "nonlinear regression" problem with "dummy"   *
*  responses all equal to zero.  See this reference for details           *
*  Here, fdot = geometric mean of all the current estimated means,        *
*  and dummy = dummy responses all = 0                                    *
*                                                                         *
**************************************************************************/

data outnlin; set outnlin;
  if pred<=0 then delete;
  sigma=sqrt(sigma/7);    * lazy 7 = n-p;
  lpred=log(pred);

proc means data=outnlin noprint;
  var lpred;
  output out=outmeans mean=mlog;

data outnlin; merge outnlin outmeans(drop=_type_ _freq_);
  retain xx 0;
  if _n_=1 then xx=mlog; 
  if mlog=. then mlog=xx;

data outnlin; set outnlin;
  dummy=0;
  fdot=exp(mlog);

/*  relax the convergence criterion for variance parameters */

proc nlin data=outnlin method=gauss converge=0.00001 &opt;
  parms delta = 0 to 1.5 by 0.05;
  trk = resid*((fdot/pred)**delta);
  model dummy = trk;
  der.delta = resid*((fdot/pred)**delta)*log(fdot/pred);
  output out=new parms=delta;
run;

/**************************************************************************
*                                                                         *
*  Set up the data set containing the new predicted values, etc, for use  *
*  on the next iteration.  Also, create the data set "parms" containing   *
*  only the parameter estimates from the iteration just completed for     *
*  printing                                                               *
*                                                                         *
**************************************************************************/

data new; set new(keep = &xvar &yvar b1 b2 b3 b4 pred delta sigma
   rename=(b1=b1_new b2=b2_new b3=b3_new b4=b4_new pred=pred1
           delta=delta1 sigma=sigma1));

data parms; set new(drop = &xvar &yvar); 
  if _n_>1 then delete;

proc print data=parms;

data &dset; set new;

%mend step3;

/**************************************************************************
*                                                                         *
*  End of the macro definitions.  Now the program to read the data and    *
*  call the macros begins.                                                *
*                                                                         *
**************************************************************************/

data indo; input time conc; 
  cards;
0.25 2.05
0.50 1.04
0.75 0.81
1.00 0.39
1.25 0.30
2.00 0.23
3.00 0.13
4.00 0.11
5.00 0.08
6.00 0.10
8.00 0.06
;

/**************************************************************************
*                                                                         *
*  Call the macro "glsalg" with the name of the indomethacin data set,    *
*  the names of x (time) and y (conc), and the starting values.           *
*                                                                         *
**************************************************************************/

title1 "GLS ALGORITHM APPLIED TO THE INDOMETHACIN DATA";

%glsalg(indo,time,conc,0.69,0.69,-1.6,-1.6);

