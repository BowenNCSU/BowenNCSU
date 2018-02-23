/*******************************************************************

  CHAPTER 5, HIP REPLACEMENT STUDY

  Population-averaged model

  -  the repeated measurement factor weeks
  -  there is one "group" factor, gender (0=female, 1 = male)
  -  an additional among-individual covariate, age, is also available
  -  the response is haematocrit 

  These data are unbalanced both in the sense that some patients
  were not observed at all times.  We demonstrate below how to make
  sure PROCE MIXED does the bookkeeping to keep track of this.

  As discussed in the notes, the missing values are almost entirely at
  week 2, suggesting that something that had nothing to do with
  the evolution of haematocrit or the gender or age of the patients, so
  that it may be reasonable to assume MCAR missingness.  However, to be
  safe, we use maximum likelihood as discussed below.
    
*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  Read in the data set 

*******************************************************************/

data hips; infile 'hips.dat';
  input patient gender age week h;
  week2=week*week;
  time=week;    *   this variable will be used for bookkeeping below;

/*******************************************************************

  Use PROC MIXED to fit a PA quadratic model with different
  assumptions about the overall covariance matrix.  We first fit the 
  basic quadratic model (5.22) to assess covariance structure.  We
  restrict attention to covariance structures that are the same for 
  each gender for illustration.    
    
  With all but the compound symmetric structure, we have to be
  careful to communicate to PROC MIXED the fact that the data
  are imbalanced in the sense that the times are all the same
  for all patients, but some patients are not observed at some
  of the times.  In our mean model, we want WEEK, the time factor,
  to be continuous; however, PROC MIXED needs also for the time
  factor to be a classification factor so that it can properly figure out 
  the missingness pattern. We give it this information by defining
  TIME = WEEK and letting TIME be a classification factor in the
  REPEATED statement.

  Because of the missing values, we use ML here instead of REML (the
  default).  The model-based standard errors are used by default.
  Based on the discussion in Section 5.6, we should use these SEs rather
  than the robust sandwich standard errors.  In general, the robust SEs
  can be obtained by specifying the EMPIRICAL option in the PROC MIXED
  statement.
    
*******************************************************************/

*  unstructured;

title "FIT WITH UNSTRUCTURED COMMON COVARIANCE";
proc mixed data=hips method=ML;
  class patient time gender;
  model h = gender gender*week gender*week2 / noint solution chisq;
  repeated time / type = un subject=patient r= 1,10,15 rcorr=1,10,15;
run;

*  compound symmetry;

title "FIT WITH COMMON COMPOUND SYMMETRY";
proc mixed data=hips method=ML;
  class patient time gender;
  model h = gender gender*week gender*week2 / noint solution chisq;
  repeated time / type = cs subject=patient r= 1,10,15 rcorr=1,10,15;
run;

*  ar(1);

title "FIT WITH COMMON AR(1) STRUCTURE";
proc mixed data=hips method=ML;
  class patient time gender;
  model h = gender gender*week gender*week2  / noint solution chisq;
  repeated time / type = ar(1) subject=patient r= 1,10,15 rcorr=1,10,15;
run;

*  one-dependent;

title "FIT WITH COMMON ONE-DEPENDENT STRUCTURE";
proc mixed data=hips method=ML;
  class patient time gender;
  model h = gender gender*week gender*week2 / noint solution chisq covb;
  repeated time / type = toep(2) subject=patient r= 1,10,15 rcorr=1,10,15;
run;

/******************************************************************

  AIC and BIC for these four fits:

                     AIC (Smaller is Better)         574.4
                     AIC (Smaller is Better)         573.4
                     AIC (Smaller is Better)         573.4
                     AIC (Smaller is Better)         573.0
 
                     BIC (Smaller is Better)         596.8
                     BIC (Smaller is Better)         584.6
                     BIC (Smaller is Better)         584.6
                     BIC (Smaller is Better)         584.2

  According to AIC and BIC, the one-dependent structure is preferred; 
  however, there is very little difference in the fits with any of
  the strutures.

  Now we fit the model in (5.23) of the notes, which allows baseline mean 
  haematocrit to depend on patient's age in a way that is different for 
  each gender.  We could of course do something similar for the linear
  and quadratic terms.  From the results, the evidence does not suggest
  dependence of baseline mean haematocrit on age, but we go ahead 
  and illustrate use of contrast and estimate statements involving 
  age.
    
*******************************************************************/

title "FIT WITH COMMON ONE-DEPENDENT STRUCTURE";
proc mixed data=hips method=ML;
  class patient time gender;
  model h = gender gender*age gender*week gender*week2 / noint solution chisq covb;
  repeated time / type = toep(2) subject=patient r= 1,10,15 rcorr=1,10,15;

*  for illustration, we look at diff in mean response at 3 weeks between;  
*  males and females who are 65 years old and get the Wald test;

contrast 'f vs m 65 yo, wk 3' gender 1 -1 gender*age 65 -65 
                          gender*week 3 -3 gender*week2 9 -9 /chisq;

*  and estimate the mean for 65 year old females at week 2.5;

estimate 'f, 65 yo, wk 2.5' gender 1 0 gender*age 65 0 gender*week 2.5 0 gender*week2 6.25 0;

run;
