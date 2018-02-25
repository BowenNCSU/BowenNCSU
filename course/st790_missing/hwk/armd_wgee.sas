/*********************************************************************

  Age-related macular degeneration data - analyses using weighted
  GEEs.

  There are 8 individuals with nonmonotone patterns; all the rest
  exhibit monotone dropout.  The baseline measure is available on
  everyone.

  Several analyses are performed here on the visual acuity measures
  (visual0 - visual52 in the data set).

  (a) Fit a multivariate model to the available data at 0, 
      4, 12, 24, 52 weeks for each treatment and common compound symmetric
      covariance matrix for both treatments using proc gee 

  (b) Fit the same model to only the data on individuals who exhibit
      monotone dropout (there are 8 individuals whose patterns are
      not monotone, so delete them) using proc gee
    
  (c) Fit the same model to the individuals with monotone dropout
      using WGEE methods with subject level and occasion level weighting

  (d) Using multiple imputation with the monotone option in the mcmc
      statement in proc mi to "fill in" the nonmonotone missingness
      patterns for the 8 individuals so that the imputed data exhibit solely
      monotone dropout, and repeat the analyses in (c)  

***********************************************************************/

options ls=80; run;

/*********************************************************************

  Read in the data set

***********************************************************************/

*  Read in the data -- we are interested in the visual acuity measures;

data armd; infile "armd.hwk4.dat";
    input pid line0 lost4 lost12 lost24 lost52 visual0 visual4 visual12 visual24 
    visual52 lesion trt;
run;

/* Recode treatment and drop the variables not of interest and keep the lesion
   variable for use as a covariate in the dropout hazards later.  Also
   identify individuals with nonmonotone patterns for later  */

data armd; set armd;
   if trt = 1 then treat = 0;  *  placebo;
   if trt = 4 then treat = 1;  *  active treatment;

*  identify those with nonmonotone missing pattern;
   
   nonmono=0;
   if visual4=. and visual12 ne . then nonmono=1;
   if visual12=. and visual24 ne . then nonmono=1;   
   if visual24=. and visual52 ne . then nonmono=1;

*  dichotomize lesion;

   les=0;
   if lesion>2 then les=1;
   
   drop trt line0 lost4 lost12 lost24 lost52;
run;
   
/*
proc sort data=armd; by treat; run;
proc means data=armd; by treat; var visual52; run;
proc means data=armd; by treat; where nonmono=0; var visual52; run; 
*/
/*
proc print data=armd;
    var pid nonmono visual0 visual4 visual12 visual24 visual52 les treat;
run;
*/

/**********************************************************************

   Reconfigure the data
      
***********************************************************************/

data armd_alt; set armd;
    array vis(5) visual0 visual4 visual12 visual24 visual52;
    do time = 1 to 5;
        visual = vis(time);
        output;
        end;
    drop visual0 visual4 visual12 visual24 visual52;
run;

data armd_alt; set armd_alt;

    if time = 1 then week = 0;
    if time = 2 then week = 4;
    if time = 3 then week = 12;
    if time = 4 then week = 24;
    if time = 5 then week = 52;

    drop time;
run;

/**********************************************************************

   Carry out (a) the available case analysis and (b) the same analysis
   after deleting the 8 individuals with nonmonotone dropout
      
***********************************************************************/

title "(a) AVAILABLE CASE";
proc gee data=armd_alt;
   class week pid;
   model visual = week week*treat / noint dist=normal link=identity;
   repeated subject=pid / type=cs within=week corrw;
run;

*  delete the individuals with nonmonotone dropout patterns and;
*  create additional variables for use in dropout hazard model;
*  in analyses (c) and (d);  

data armd_mono; set armd_alt;

    if nonmono=1 then delete;

*  lag-1 visual acuity for use in dropout hazard;    
    
    prevvis=lag(visual);
    if week=0 then prevvis=1;

*  Dummy variables for week to create separate dropout hazard for;
*  each week;

    w4=0; w12=0; w24=0; w52=0;
    if week=4 then w4=1;
    if week=12 then w12=1;
    if week=24 then w24=1;
    if week=52 then w52=1;
    
    prev4=prevvis*w4; prev12=prevvis*w12; prev24=prevvis*w24; prev52=prevvis*w52;

run;

/*proc print data=armd_mono; run;*/

title "(b) MONOTONE DROPOUT ONLY";
proc gee data=armd_mono;
   class week pid;
   model visual = week week*treat / noint dist=normal link=identity;
   repeated subject=pid / type=cs within=week corrw;
run;

/**********************************************************************

   Carry out (c) weighted analyses on the monotone only dataset

***********************************************************************/

title "(c) SUBJECT LEVEL WEIGHTING";
proc gee data=armd_mono;
   class week pid;
   model visual = week week*treat / noint dist=normal link=identity;
   repeated subject=pid / type=cs within=week corrw;
   missmodel w12 w24 w52 prev4 prev12 prev24 prev52 les / type=sublevel;
run;

title "(c) OCCASION LEVEL WEIGHTING";
proc gee data=armd_mono;
   class week pid;
   model visual = week week*treat / noint dist=normal link=identity;
   repeated subject=pid / type=cs within=week corrw;
   missmodel w12 w24 w52 prev4 prev12 prev24 prev52 les / type=obslevel;
run;

/**********************************************************************

   (d)(i) Call proc mi to "fill in" the nonmonotone patterns via propoer
   imputation using MCMC and create M=10 imputed data sets.  This is
   accomplished using impute=monotone in the mcmc statement.
      
***********************************************************************/

proc mi data=armd out=armdmonoout seed=1518971 nimpute=10;
   mcmc impute=monotone;  
   var visual0 visual4 visual12 visual24 visual52;
run;


/**********************************************************************

    Print out the imputed data sets (first 20 observations of first one)

***********************************************************************/

title " "; 
proc print data=armdmonoout (obs=20);
    var _imputation_ visual0 visual4 visual12 visual24 visual52;
run;

/**********************************************************************

    Now create reconfigured data sets with the variables for the
    dropout hazards as above

***********************************************************************/

data armdmono_alt; set armdmonoout;
    array vis(5) visual0 visual4 visual12 visual24 visual52;
    do time = 1 to 5;
        visual = vis(time);
        output;
        end;
    drop visual0 visual4 visual12 visual24 visual52;
run;

data armdmono_alt; set armdmono_alt;
    if time = 1 then week = 0;
    if time = 2 then week = 4;
    if time = 3 then week = 12;
    if time = 4 then week = 24;
    if time = 5 then week = 52;

*  lag-1 visual acuity for use in dropout hazard;    
    
    prevvis=lag(visual);
    if week=0 then prevvis=1;

*  Dummy variables for week to create separate dropout hazard for;
*  each week;

    w4=0; w12=0; w24=0; w52=0;
    if week=4 then w4=1;
    if week=12 then w12=1;
    if week=24 then w24=1;
    if week=52 then w52=1;
    
    prev4=prevvis*w4; prev12=prevvis*w12; prev24=prevvis*w24; prev52=prevvis*w52;
 
    drop time;
run;

proc print data=armdmono_alt (obs=20);
    var _imputation_ pid nonmono week visual treat les;
run;

proc sort data=armdmono_alt;
    by _imputation_ pid week;
run;


/**********************************************************************

    (d)(ii) Run the weighted GEE analyses on each imputed data set and
    obtain mean estimates.  Use ODS to output the need results
    (see the proc gee documentation for the ODS variable
    names it produces).  For each imputed data set, we output data
    sets containing the mean parameter estimates (meanparms) snf the
    asymptotic covariance matrices of the mean parameter estimates
    (meancov).

***********************************************************************/

*  Using subject level weighting;

title "(d)(ii) SUBJECT LEVEL IMPUTATION";
proc gee data=armdmono_alt;
   by _imputation_; 
   class week pid;
   model visual = week week*treat / noint dist=normal link=identity;
   repeated subject=pid / type=cs within=week corrw ecovb;
   missmodel w12 w24 w52 prev4 prev12 prev24 prev52 les / type=sublevel;
  ods output GEEEmpPEst=subparms ParmInfo=subinfo GEERCov=subcovb;
run;

*  Using occasion level weighting;

title "(d)(ii) OCCASION LEVEL IMPUTATION"; 
proc gee data=armdmono_alt;
   by _imputation_; 
   class week pid;
   model visual = week week*treat / noint dist=normal link=identity;
   repeated subject=pid / type=cs within=week corrw ecovb;
   missmodel w12 w24 w52 prev4 prev12 prev24 prev52 les / type=obslevel;
  ods output GEEEmpPEst=obsparms ParmInfo=obsinfo GEERCov=obscovb;
run;

*  Print the output data sets to have a look at them;


title " ";
/*
proc print data=subparms; run;
proc print data=obsparms; run;
proc print data=subcovb; run;
proc print data=obscovb; run;
proc print data=subinfo; run;
proc print data=obsinfo; run;
    */

/**********************************************************************

   (d)(iii) Call proc mianalyze to combine the results for the mean for
   each analysis  
    
***********************************************************************/

title "(d)(iii) MI SUBJECT LEVEL ANALYSIS"; 
proc mianalyze parms(classvar=level)=subparms parminfo=subinfo covb=subcovb wcov
    bcov;
    class week;
    modeleffects week week*treat;
run;

title "(d)(iii) MI OCCASION LEVEL ANALYSIS"; 
proc mianalyze parms(classvar=level)=obsparms parminfo=obsinfo covb=obscovb wcov
    bcov;
    class week;
    modeleffects week week*treat;
run;








