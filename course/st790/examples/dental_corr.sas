/*******************************************************************

  CHAPTER 2, EXAMPLE 1, Dental Study

  Explore overall and within-child patterns of correlation in the
  dental study data using SAS

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data are not in the correct from for use with the SAS procedures
  CORR and DISCRIM we use below.  These procedures require that the
  data be in the form of one record per individual.  The data in 
  the file dental.dat are in the form of one record per observation.
 
  In particular, the data set looks like

  1 1 8 21 0
  2 1 10 20 0
  3 1 12 21.5 0
  4 1 14 23 0
  5 2 8 21 0
      .
      .
      .

  column 1    observation number
  column 2    child id number
  column 3    age
  column 4    response (distance)
  column 5    gender indicator (0=girl, 1=boy)

  We thus create a new data set such that each record in the data
  set represents all 4 observations on each child plus gender
  identifier using PROC TRANSPOSE and use PROC PRINT to print out
  the first 5 records (so data for the first 5 children, all girls).

*******************************************************************/

data dent1; infile 'dental.dat';
  input obsno child age distance gender;
run;

proc transpose data=dent1 out=dent2 prefix=age;
   by gender child notsorted;
   var distance; 
run;

title "TRANSFORMED DATA -- 1 RECORD/INDIVIDUAL";
proc print data=dent2(obs=5); run;

/*  Alternatively, a more cumbersome way to do this is as follows 

data dent1; set dent1;
  if age=8 then age=1;
  if age=10 then age=2;
  if age=12 then age=3;
  if age=14 then age=4;
  drop obsno;
run;

proc sort data=dent1;
  by gender child;
run;

data dent2(keep=age1-age4 gender child);
  array aa{4} age1-age4;
  do age=1 to 4;
  set dent1;
  by gender child;
  aa{age}=distance;
  if last.child then return;
end;
run;

*/

/*******************************************************************

  Use  PROC CORR to obtain the sample means at each
  age (the means of the variables AGE1,...,AGE4 in DENT2) and to
  calculate the sample covariance matrix and corresponding sample
  correlation matrix separately for each group (girls and boys).
  The COV option in the PROC CORR statement asks for the sample
  covariance to be printed.  

*******************************************************************/

proc sort data=dent2; by gender; run;

title "SAMPLE COVARIANCE AND CORRELATION MATRICES BY GENDER";
proc corr data=dent2 cov;
 by gender; var age1 age2 age3 age4; 
run;

/*******************************************************************

  We now obtain the "centered" and "scaled" values that can be used
  plotting scatterplot matrices using PROC STDSIZE.  We export
  the centered/scaled values to a file so that we could import
  them into R to use the function pairs() to plot scatterplots.  

*******************************************************************/  

proc sort data=dent2; by gender child; run;

proc stdize data=dent2 method=std pstat out=dentstats(rename=(age1=csage1
      age2=csage2 age3=csage3 age4=csage4));
  var age1-age4;
  by gender;
run;

/*  Alternatively, we can do this manually as follows

proc means data=dent2 mean std noprint; by gender;
  var age1 age2 age3 age4;
  output out=dentstats mean=mage1 mage2 mage3 mage4
         std=sdage1 sdage2 sdage3 sdage4;
run;

data dentstats; merge dentstats dent2; by gender;
  csage1=(age1-mage1)/sdage1;
  csage2=(age2-mage2)/sdage2;
  csage3=(age3-mage3)/sdage3;
  csage4=(age4-mage4)/sdage4;
run;  */

title "CENTERED AND SCALED DATA BY GENDER";
proc print data=dentstats(obs=3); run;

*  create the data file -- use only the variables we need;

filename this "dentcenter.dat";
data _null_; set dentstats;
  file this;
  put gender csage1 csage2 csage3 csage4;
run;

/*******************************************************************

  One straightforward way to have SAS calculate the pooled sample
  covariance matrix and the corresponding estimated correlation matrix
  is using PROC DISCRIM.  This procedure is focused on so-called
  discriminant analysis, which is discussed in a standard text on
  general multivariate analysis.  The data are considered as
  in the form of vectors; here, the elements of a data vector are
  denoted as AGE1,...,AGE4.  We disregard other portions of the output.

*******************************************************************/

title "OBTAINING THE POOLED SAMPLE COVARIANCE MATRIX AND CORRELATION";
title2 "USING PROC DISCRIM";
proc discrim pcov pcorr data=dent2;
  class gender;
  var age1 age2 age3 age4;
run;

/*******************************************************************
  
  Another way to obtain the pooled sample covariance matrix 
  and the associated correlation matrix is to use PROC MIXED, which
  is discussed in Chapters 5-6.  PROC MIXED requires the data to
  be in the form of one observation per record, which is the original
  form in the data set DENT1.  The R and RCORR options in the REPEATED
  statement cause PROC MIXED to print out the estimates.  The MODEL
  statement specifies a different mean for each gender-age combination,
  so that the gender-specific means are correctly used in the 
  calculation.  A call to PROC MIXED of the form of that below will 
  provide the pooled sample covariance matrix as long as the data are
  balanced.  Again, we ignore the rest of the output.  

*******************************************************************/

title "OBTAINING THE POOLED SAMPLE COVARIANCE MATRIX AND CORRELATION";
title2 "USING PROC MIXED";
proc mixed data=dent1;
  class gender age child;
  model distance = gender age gender*age;
  repeated / type=un subject=child r rcorr;
run;

/*******************************************************************

  Manually create data sets of lag 1, lag 2, and lag 3 data and 
  then calculate the autocorrelation functions by gender using PROC
  of the autocorrelation function; however, for longitudinal data sets
  where the number of time points is small, the "manual" approach
  we have demonstrated here is easy to implement and understand.  

*******************************************************************/

data lag1; set dentstats;
  by child;
  pair1=csage1; pair2=csage2; output;
  pair1=csage2; pair2=csage3; output;
  pair1=csage3; pair2=csage4; output;
  if last.child then return;
  drop csage1-csage4;
run;

title "AUTOCORRELATION FUNCTION AT LAG 1";
proc print data=lag1(obs=6); run;
proc sort data=lag1; by gender;

proc corr data=lag1; by gender; 
  var pair1 pair2; 
run;

data lag2; set dentstats;
  by child;
  pair1=csage1; pair2=csage3; output;
  pair1=csage2; pair2=csage4; output;
  if last.child then return;
  drop csage1-csage4;
run;

title "AUTOCORRELATION FUNCTION AT LAG 2";
proc print data=lag2(obs=6); run;
proc sort data=lag2; by gender;

proc corr data=lag2; by gender; 
  var pair1 pair2; 
run;

data lag3; set dentstats;
  by child;
  pair1=csage1; pair2=csage4; output;
  if last.child then return;
  drop csage1-csage4;
run;

title "AUTOCORRELATION FUNCTION AT LAG 3";
proc print data=lag3(obs=6); run;
proc sort data=lag3; by gender;

proc corr data=lag3; by gender; 
  var pair1 pair2; 
run;




