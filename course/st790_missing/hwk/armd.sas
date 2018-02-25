/*****************************************************************

  Age-related macular degeneration data - ML analyses under MAR

  There are 8 individuals with nonmonotone patterns; all the rest
  exhibit monotone dropout.  The baseline measure is available on
  everyone.

  Several analyses are performed here on the visual acuity measures
  (visual0 - visual52 in the data set).

  (a)  Summarize missing data patterns (see output to proc mi)

  (b-d) Fit a multivariate normal model to the available data, LOCF
      data, and complete cases with distinct means at each of 0,
      4, 12, 24, 52 weeks for each treatment and common unstructured
      covariance matrix for both treatments by ML using proc mixed 

  (f) Fit a multivariate normal model to the available data separately
      by treatment group by ML using proc mixed 
    
  (g) Fit the same models to all available data using EM algorithm
      in proc mi

*****************************************************************/

options ls=80; run;

*  Read in the data -- we are interested in the visual acuity measures;

data armd; infile "armd.dat";
    input pid line0 lost4 lost12 lost24 lost52 visual0 visual4 visual12 visual24 
    visual52 lesion trt;
run;

*  Recode treatment and drop the variables not of interest;

data armd; set armd;
   if trt = 1 then treat = 0;  *  placebo;
   if trt = 4 then treat = 1;  *  active treatment;
   drop trt line0 lost4 lost12 lost24 lost52 lesion; 
run;

/*
  Create a LOCF data set where missing values are filled in by the
  last available value
*/

data armd_locf; set armd;
    if visual4 = . then visual4 = visual0;
    if visual12 = . then visual12 = visual4;
    if visual24 = . then visual24 = visual12;
    if visual52 = . then visual52 = visual24;
run;

/* Create a complete case data set where individuals with any variable
    missing are deleted
*/

data armd_cc; set armd;
    if visual4 = . then delete;
    if visual12 = . then delete;
    if visual24 = . then delete;
    if visual52 = . then delete;
run;

/*
  We are interested in the variables visual0 - visual52.  To fit a
  normal model using proc mixed, we need to transform each data set from
  1 record per individual to 1 record per observation
*/
    
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

*  Do the same for the LOCF data set;

data armd_locf_alt; set armd_locf;
    array vis(5) visual0 visual4 visual12 visual24 visual52;
    do time = 1 to 5;
        visual = vis(time);
        output;
        end;
    drop visual0 visual4 visual12 visual24 visual52;
run;

data armd_locf_alt; set armd_locf_alt;
    if time = 1 then week = 0;
    if time = 2 then week = 4;
    if time = 3 then week = 12;
    if time = 4 then week = 24;
    if time = 5 then week = 52;
    drop time;
run;   

*  Do the same for the CC data set;

data armd_cc_alt; set armd_cc;
    array vis(5) visual0 visual4 visual12 visual24 visual52;
    do time = 1 to 5;
        visual = vis(time);
        output;
        end;
    drop visual0 visual4 visual12 visual24 visual52;
run;

data armd_cc_alt; set armd_cc_alt;
    if time = 1 then week = 0;
    if time = 2 then week = 4;
    if time = 3 then week = 12;
    if time = 4 then week = 24;
    if time = 5 then week = 52;
    drop time;
run;   

*  (b) Fit multivariate normal model to available data from both groups using ML;

title "Available case analysis, both treatments";
proc mixed data=armd_alt method=ml;
    class  week pid;
    model visual = week week*treat / noint solution;
    repeated week / subject=pid type=un r=2 rcorr=2;
run;

*  (c) Fit same multivariate normal model to complete cases using ML;

title "Complete case analysis, both treatments";
proc mixed data=armd_cc_alt method=ml;
    class  week pid;
    model visual = week week*treat / noint solution;
    repeated week / subject=pid type=un r=2 rcorr=2;
run;

*  (d) Fit same multivariate normal model to LOCF data using ML;

title "LOCF analysis, both treatments";
proc mixed data=armd_locf_alt method=ml;
    class  week pid;
    model visual = week week*treat / noint solution;
    repeated week / subject=pid type=un r=2 rcorr=2;
run;

*  (f) Now fit multivariate normal model to avaiable data by treatment;

proc sort data=armd_alt; by treat; run;

title "Separate fits by treatment, available data, ML via proc mixed";
proc mixed data = armd_alt method=ml; by treat;
    class week pid;
    model visual = week /noint solution;
    repeated week / subject=pid type=un r=2 rcorr=2;
run;

*  (g) Fit the same models to each treatment using the EM algorithm;

proc sort data=armd; by treat; run;

    title "Separate fits by treatment, available data, ML via EM algorithm";
proc mi data=armd seed=370252 simple nimpute=0; by treat;
    em initial=ac maxiter=200 converge=1e-4 itprint;
    var visual0 visual4 visual12 visual24 visual52;
run;




  


      


