/*********************************************************************

    Use proc gee to carry out WGEE analyses.  The data are

    h1, h2, h4    hypertensive status (0 = no, 1 = yes) at 1, 2, and 4
                  months after starting treatment
    age           age range (1 = < 40, 2 = 40 - 65, 3 = > 65)
    chol          cholesterol (mg/DL) at baseline

    The covariates age and chol are always observed, while hypertensive
    status can be missing due to dropout.

    The data here exhibit monotone missingess (dropout) -- all
    individuals with missing data have monotone dropout patterns.
    
***********************************************************************/

options ls=80 nodate; run;

/*********************************************************************

  Read in the data set

***********************************************************************/

data hyper;  infile "hyper_mono.dat";
    input h1 h2 h4 age chol;
run;

* create an id indiator for use later;

data hyper; set hyper;
    id=_n_;
run;

/**********************************************************************

   First reconfigure the data set to have 1 record per observation.
    
***********************************************************************/

data hyper_alt; set hyper;
    array y(3) h1 h2 h4;
    do time = 1 to 3;
        hyper = y(time);
        output;
        end;
    drop h1 h2 h4;
run;

data hyper_alt; set hyper_alt;
    month=time;
    if time=3 then month=4;

*  Create dummies for age manually;
    
    age2=0;  if age=2 then age2=1;
    age3=0;  if age=3 then age3=1;

/*  Get the lag-1 hyper values for use in the dropout hazard model.
    Set the lag variable = 1 for the first time point; its value
    doesn't get used */
    
    prevhyper=lag(hyper);
    if month=1 then prevhyper=1;

/*  To get a separate dropout model for each time point, create
    indicator variables for each time point with dropout and
    associated covariates.  This will allow us to fit a global
    hazard model that embeds a separate hazard for each time point
    with time-specific parameters in the missmodel statement in
    proc gee.  We create dummies for months 2 and 4 so that month
    1 is the reference */

    m2=0; m4=0;
    if month=2 then m2=1;
    if month=4 then m4=1;

    age2_2=age2*m2; age2_4=age2*m4;
    age3_2=age3*m2; age3_4=age3*m4;
    chol_2=chol*m2; chol_4=chol*m4;
    prev_2=prevhyper*m2;
    prev_4=prevhyper*m4;

run;
        
proc print data=hyper_alt (obs=20);
    var id time month hyper prevhyper age age2 age3 chol;
run;

/**********************************************************************

    Now run several analyses.  We are interested in whether or not
    there is a time trend in hypertension status, adjusting for age
    and cholesterol at baseline.  We assume an unstructured working
    correlation matrix.

    First run the available case analysis in both proc genmod and 
    proc gee.

    The syntax for proc gee is similar to that for proc genmod.  The
    descending option in the proc gee statement works the same way as
    in proc genmod, and the model statement syntax and choices of
    distribution and link are the same.  The repeated statement is the
    same -- in both procs, "within" or "withinsubject" specify the ordering of
    observations when constructing a working unstructured or other
    correlation structure where the ordering is importantl either one
    can be used.  The corrw option prints out the estimated working
    correlation matrix.  
    
***********************************************************************/

title "AVAILABLE CASE ANALYSIS, PROC GENMOD"; 
proc genmod data=hyper_alt descending;
    class id time;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un withinsubject=time corrw;
run;

title "AVAILABLE CASE ANALYSIS, PROC GEE"; 
proc gee data=hyper_alt descending;
    class id time;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un within=time corrw;
run;

/**********************************************************************

    Now run the weighted analyses using proc gee.     

    The missmodel statement specifies the effects that go into a
    the linear predictor of a logistic regression model for the
    cause-specific hazards of not dropping out at each time point.
    Here, we model those as depending on the baseline covariates and
    the value of hyper at the previous time point, with separate
    parameters for the model at each time point.

    It is assumed that there are no missing data at the first time
    point (which may or may not be baseline).  If any individual has
    all outcomes missing at all time points, this is considered NOT
    to be monotone missingness, and proc gee will generate an error.
    
***********************************************************************/

*  Using subject level weighting;

title "WEIGHTED ANALYSIS, SUBJECT LEVEL WEIGHTING"; 
proc gee data=hyper_alt descending;
    class id time;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un within=time corrw;
    missmodel m4 age2_2 age2_4 age3_2 age3_4
        chol_2 chol_4 prev_2 prev_4 / type=sublevel;
run;

*  Using occasion level weighting;

title "WEIGHTED ANALYSIS, OCCASION LEVEL WEIGHTING"; 
proc gee data=hyper_alt descending;
    class id time;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un within=time corrw;
    missmodel m4 age2_2 age2_4 age3_2 age3_4
        chol_2 chol_4 prev_2 prev_4 / type=obslevel;
*    missmodel age2 age3 chol prevhyper / type=obslevel;
run;

/*********************************************************************

  Because these are simulated data, we also have access to the full
  data before dropout was induced.  Thus, we do the full data analysis
  for comparision to the weighted analyses.  Read in the data and
  reconfigure.

***********************************************************************/

title " ";
data hyperfull;  infile "hyper_mono_full.dat";
    input h1 h2 h4 age chol;
run;

data hyperfull; set hyperfull;
    id=_n_;
run;

data hyperfull_alt; set hyperfull;
    array y(3) h1 h2 h4;
    do time = 1 to 3;
        hyper = y(time);
        output;
        end;
    drop h1 h2 h4;
run;

data hyperfull_alt; set hyperfull_alt;
    month=time;
    if time=3 then month=4;

*  Create dummies for age manually;
    
    age2=0;  if age=2 then age2=1;
    age3=0;  if age=3 then age3=1;

run;

proc print data=hyperfull_alt (obs=20);
    var id time month hyper age age2 age3 chol;
run;

/*********************************************************************

   Do the analysis using both proc genmod and proc gee.
    
***********************************************************************/

title "FULL DATA ANALYSES"; 
proc genmod data=hyperfull_alt descending;
    class id time;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un withinsubject=time corrw;
run;

proc gee data=hyperfull_alt descending;
    class id time;
    model hyper = age2 age3 month chol / dist=bin link=logit;
    repeated subject=id / type=un within=time corrw;
run;
