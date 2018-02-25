/*********************************************************************

    Use proc mi to carry out multivariate imputation via chained
    equations on data from a hypertension study.  The variables are:

    h1, h2, h4    hypertensive status (0 = no, 1 = yes) at 1, 2, and 4
                  months after starting treatment
    age           age range (1 = < 40, 2 = 40 - 65, 3 = > 65)
    chol          cholesterol (mg/DL) at baseline

***********************************************************************/

options ls=80 nodate; run;

/*********************************************************************

  Read in the data set

***********************************************************************/

data hyper;  infile "hyper.dat";
    input h1 h2 h4 age chol;
run;

* create an id indiator for use later;

data hyper; set hyper;
    id=_n_;
run;

/**********************************************************************

   Call proc mi to carry out imputation using MICE via the FCS (Fully 
   Conditional Specification) statement.  Here, the variables to be
   imputed are binary, categorical, and continuous.  We use M=10
   imputated data sets.  We then fit a population average model to the 
   longitudinal binary data using proc genmod and use mianalyze to
   summarize the results.

   proc mi statement:  out=hyper.out specifies the output data set
   containing the results of the imputation, and nimpute=10 requests M=10
   imputed data sets.

    Each FCS statement specifies that the given variable should be 
    regressed on all others.  This is the default -- it is also possible
    to choose to regress on a subset of variables.  E.g.,

    fcs reg(chol);

    regresses chol on all other variables;

    fcs reg(chol = h1 h2 h4);

    regresses it only on these three variables.

    We demonstrate use of the reg method for continuous variables and
    the logistic method for binary (logit link) or categorical (cumulative
    logit link)

    The var statement specifies the variables to be imputed.    

    And see the documentation for more options in general.
    
***********************************************************************/

proc mi data=hyper out=hyperout nimpute=10 seed=1591842;
    class h1 h2 h4 age;
    var h1 h2 h4 age chol;
    fcs logistic(h1);
    fcs logistic(h2);
    fcs logistic(h4);
    fcs logistic(age);
    fcs reg(chol);
run;

/**********************************************************************

    Print out the imputed data sets (first 20 obsevations of first one)

***********************************************************************/

proc print data=hyperout (obs=20);
    var _imputation_ h1 h2 h4 age chol;
run;

/**********************************************************************

    Now create reconfigured data sets with 1 record per observation
    to illustrate how to call another proc to analyze each imputed
    data set.  Sort the data sets for use by proc genmod next.

***********************************************************************/

data hyperout_alt; set hyperout;
    array y(3) h1 h2 h4;
    do time = 1 to 3;
        hyper = y(time);
        output;
        end;
    drop h1 h2 h4;
run;

data hyperout_alt; set hyperout_alt;
    month=time;
    if time=3 then month=4;

*  Create dummies for age manually to avoid weirdness with proc;
*  mianalyze;   
    
    age2=0;  if age=2 then age2=1;
    age3=0;  if age=3 then age3=1;
run;
        
proc print data=hyperout_alt (obs=20);
    var _imputation_ id time month hyper age age2 age3 chol;
run;

proc sort data=hyperout_alt;
    by _imputation_ id month;
run;

/**********************************************************************

    Now run proc genmod on each imputed data set.  We are interested
    in whether or not there is a time trend in hypertension status, 
    adjusting for age and cholesterol at baseline.

    Use ODS to output the needed results (see the proc genmod documentation 
    for the output variable names it produces).  Here, because we are
    fitting an unstructured correlation matrix, using the robust sandwich
    covariance matrix seems like overkill, so we output the model-based
    covariance matrices for use by proc mianalyze.  To output robust
    covariance matrices, use covb in place of mcovb in the repeated
    statement and replace GEENCov=gmcovb by CovB=gmcovb in the ods statement.
    
***********************************************************************/

proc genmod data=hyperout_alt descending;
    class id time;
    by _imputation_;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un withinsubject=time modelse mcovb;
    ods output GEEModPEst=gmparms ParmInfo=gmpinfo GEENCov=gmcovb;
run;

*  Print the output data sets to have a look at them;

proc print data=gmparms;
run;

proc print data=gmcovb;
run;
 
proc print data=gmpinfo;
run;

/**********************************************************************

   Call proc mianalyze to combine the results. It is set up 
   to automatically recognize certain types of data sets.  The data set
   "gmparms" contains the estimates and associated standard errors
   for the mean parameters from each of the M=10 imputed data sets.
   With the parms=gmparms statement, it recognizes this.  "gmcovb"
   contains the asymptotic covariance matrics, and "gmpinfo" contains
   parameter info.  We specify that we want MI inference on each of
   the parameters in the model in the modeleffects statement.
    
***********************************************************************/

proc mianalyze parms=gmparms covb=gmcovb parminfo=gmpinfo wcov bcov
    tcov;
    modeleffects Intercept age2 age3 month chol;
run;
