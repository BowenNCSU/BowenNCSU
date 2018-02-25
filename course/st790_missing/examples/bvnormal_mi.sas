/*********************************************************************

    Use proc mi to carry out proper multiple imputation on bivariate 
    normal data given in EXAMPLE 2 of Chapter 3 of the notes

    We read in a simulated data set created in the R program
    bvnormal_em.R with the "NAs" replaced by "." (the SAS missing value
    indicator)

***********************************************************************/

options ls=80 nodate; run;

/*********************************************************************

  Read in the data set

***********************************************************************/

data bivariate;  infile "bvnormal.dat";
    input y1 y2;
run;

* create an id indiator for use later;

data bivariate; set bivariate;
    id=_n_;
run;

/**********************************************************************

   Call proc mi to carry out proper imputation via the multivariate
   normal using MCMC and M=10 imputed data sets.  To demonstrate the
   use of the accompanying proc mianalyze, we call proc mixed to fit
   a bivariate normal to each imputed data set and then combine the 
   results using proc mianalyze (Rubin's third "task") and get standard
   errors using Rubin's variance formula.  We use many of the default
   specifications but show these explicitly for clarity.

   proc mi statement:  out=bvn.out specifies the output data set
   containing the results of the imputation, nimpute=10 requests M=10
   imputed data sets, seed= specifies the seed for the random number
   generator for the MCMC, and the simple option displays descriptive
   statistics from the available data.

   There is also a round= option; e.g., round=0.001 in our case
   requestse that the imputed results be rounded to 3 decimal places
   so that they resemble the actual data.   However, studies by Horton
   et al. have shown that such rounding can bias results of multiple
   imputation, so this should be used with caution.

   proc mi uses MCMC by default, but if we want to change the default
   settings it uses, we have to invoke the mcmc statement; we do this to
   illlustrate.  initial=em specifies that the starting values for mean
   and covariance matrix should be found from running the EM algorithm;
   this is the default (one could specify maxiter= to set the number of
   EM iterations used -- see the em statement documentation). Other 
   options we invoke in the mcmc statement:

   chain=single -- use a single chain for all imputations (the default;
   the alternative is chain=multiple, using a separate chain for each)

   impute=full -- full is the default and is required for nonomontone
   missingness; if the missingness if monotone, impute=monotone can
   be used.
       
   nbiter=200 -- specifies the number of burn=in interations of MCMC
   before the first imputation in the chain (200 is the default)
       
   niter-100 -- specifies the number of iterations between imputations
   in a single chain

   There is also a prior= option to specify the prior information for 
   means and covariances; a noninformative prior is used by default.
   See the documentation for how to specify other priors.

   the var statement specifies the variables to be imputed.    

   And see the documentation for more options in general.
    
***********************************************************************/

proc mi data=bivariate out=bvnout seed=1518971 simple nimpute=10;
   mcmc initial=em chain=single impute=full nbiter=200 niter=100;
   var y1 y2;
run;

/**********************************************************************

    Print out the imputed data sets (first 20 obsevations of first one)

***********************************************************************/

proc print data=bvnout (obs=20);
    var _imputation_ y1 y2;
run;

/**********************************************************************

    Now create reconfigured data sets with 1 record per observation
    to illustrate how to call another proc to analyze each imputed
    data set.  Sort the data sets for use by proc mixed next.

***********************************************************************/

data bvnout_alt; set bvnout;
    array ya(2) y1 y2;
    do ind = 1 to 2;
        y = ya(ind);
        output;
        end;
run;

proc print data=bvnout_alt (obs=20);
    var _imputation_ id ind y;
run;

proc sort data=bvnout_alt;
    by _imputation_ id ind;
run;

/**********************************************************************

    Now run proc mixed on each imputed data set to obtain mean
    and covariance matrix estimates.  Use ODS to output the needed
    results (see the proc mixed documentation for the ODS variable
    names it produces).  For each imputed data set, we output data
    sets containing the mean parameter estimates (meanparms), the
    asymptotic covariance matrices of the mean parameter estimates
    (meancov), the estimates of the distinct parameters in the
    covariance matrix and their standard errors, which are computed
    because of the covtest option in the proc mixed statement
    (covparms), and the asymptotic covariance matrices of these
    covariance parameter estimates (covcov -- we actually won't use
    this; see below).

***********************************************************************/

proc mixed data=bvnout_alt asycov covtest method=ml noclprint noitprint;
    class ind id;
    by _imputation_;
    model y = ind / noint solution covb;
    repeated ind / subject=id type=un;
    ods output SolutionF=meanparms CovB=meancov CovParms=covparms AsyCov=covcov;
run;

*  Print the output data sets to have a look at them;

proc print data=meanparms; run;

proc print data=meancov; run;

proc print data=covparms; run;

proc print data=covcov; run;

/**********************************************************************

   Call proc mianalyze to combine the results for the mean.  It is set up 
   to automatically recognize certain types of data sets.  The data set
   "meanparms" contains the estimates and associated standard errors
   for the 2 mean parameters from each of the M=10 imputed data sets.
   With the parms=meanparms statement, it recognizes this.  In
   meanparms, the levels of the class variable "ind" identify each
   mean component. We specify that we want MI inference on each of these  
   mean components by identifying this class variable in the
   model statement.  
    
***********************************************************************/

proc mianalyze parms=meanparms;
    class ind;
    modeleffects ind;
run;

/**********************************************************************

   Alternatively, we can ask for full multivariate inference on the
   mean parameters.  This is done by telling proc mianalyze that we are
   are also providing the asymptotic covariance matrices for the estimates
   from each data set using the covb statement.  Instead of just
   providing MI estimates and standard errors, it will also provide
   the two full pieces of Rubin's formula, the "within" and "between"
   covariance matrices, with the wcov and bcov options.  The full
   covariance matrix, wcov + (M+1)*bcov/M, would be needed if one
   wanted, for example, to make inference on contrasts of parameters
   (e.g., the difference of the 2 means in this simple case).
    
***********************************************************************/

proc mianalyze parms=meanparms covb(effectvar=rowcol)=meancov wcov
    bcov;
    class ind;
    modeleffects ind;
run;

/**********************************************************************

   proc mianalyze is not set up to automatically do the above for
   things like covariance parameters.  To use it, we have provide
   a data set containing the individual estimates of the 3 distinct
   covariance parameters and their associated standard errors, each of
   which has its own variable name.  In the code below, we reconfigure
   the covparms data set to have 6 variables representing the 3
   distinct parameters and their standard errors.
    
***********************************************************************/

*  Reconfigure the data set;

data covparms; set covparms;
    Effect="covp";
    if CovParm="UN(1,1)" then covp=1;
    if CovParm="UN(2,1)" then covp=2;
    if CovParm="UN(2,2)" then covp=3;
    drop CovParm Subject ZValue ProbZ Effect;
run;

data covparms2(keep=c1-c3 sc1-sc3 _imputation_);
    array cc{3} c1-c3;
    array scc{3} sc1-sc3;
  do covp=1 to 3;
  set covparms;
  by _imputation_;
  cc{covp}=Estimate;
  scc{covp}=StdErr;
  if last._imputation_ then return;
end;
run;

proc print data=covparms2; run;

/**********************************************************************

   proc mianalyze will treat the variables c1-c3 as the parameters
   and the variables in the stderr statement as their associated
   standard errors.
    
***********************************************************************/
        
proc mianalyze data=covparms2;
    modeleffects c1 c2 c3;
    stderr sc1 sc2 sc3;
run;

/**********************************************************************

   Note that the above does a "univariate" analysis of the covariance
   parameters.  It is possible with more fancy coding to get the full
   asymptotic covariande matrices in the data set covcov into proc
   mianalyze do do full multivariate inference.  Alternatively, given we
   have all of these covariance matrices in a data set, the analyst could
   write his/her own code to do this instead of using proc mianalyze.
   Inference on covariance parameters is highly predicated on the
   normality assumption and often of secondary interest, so we don't
   bother with this here.
    
***********************************************************************/
