/******************************************************************

  CHAPTER 7, Fitting a logistic regression model using PROC 
  GENMOD.

  We consider data with binary response (woman has suffered a 
  myocardial infarction (0=no, 1=yes), with covariates given 
  below.
    
******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data look like (first 10 records)

     1 1 33 1 0
     2 0 32 0 0
     3 1 37 0 1
     4 0 36 0 0
     5 1 50 1 1
     6 1 40 0 0
     7 0 35 0 0
     8 1 33 0 0
     9 1 33 0 0
     10 0 31 0 0

  column 1      subject id
  column 2      oral contraceptive indicator (0=no,1=yes)
  column 3      age (years)
  column 4      smoking indicator (0=no,1=yes)
  column 5      binary response -- whether MI has been suffered
                (0=no,1=yes)

******************************************************************/

data mi; infile 'infarc.dat';
  input id oral age smoke mi;
run;

/*****************************************************************

  Fit the logistic regression model using PROC GENMOD.
  We do not use a CLASS statement here, as the covariates are
  either continuous (AGE) or already in "dummy" form (ORAL, SMOKE).
  The model statement with the LINK=LOGIT option results in the 
  logistic regression model.  The DIST=BINOMIAL
  specifies the Bernoulli distribution, which is the simplest case
  of a binomial distribution.

  PROC GENMOD models by default the probability that Y=0 rather
  than Y=1.  To make PROC GENMOD model probability that 
  Y=1, as is standard, one must include the DESCENDING option in
  the PROC GENMOD statement.  An explicit statement about 
  what is being modeled will appear in the log file and output.

  The ESTIMATE statement with the EXP option will output 
  the exponential of the specified linear combination of the 
  parameters.  We use it here to compute an odds ratio.  The logistic
  model says that the ODDS of having a MI, given a woman's covariates,
  which is the ratio of the probability of having a MI to not having
  one, are

    exp(beta_0 + beta_1 oral + beta_2 age + beta_3 smoke)
    
  The ODDS RATIO compares 2 odds - for example, if we are interested
  in comparing the odds of a MI if a randomly chosen woman smokes
  (smoke=1) relative to those if she does not (smoke=0), with all her
  other covariates fixed, from above, this ratio is exp(beta_3).
  Thus, exp(beta_3) is a multiplicative factor that measures by how much
  her odds of MI change.  If beta_3>0, this factor is >1, meaning that
  the odds increase; if beta_3<0, this factor is <1, meaning that the odds
  decrease.  beta_3 is referred to as the LOG ODDS RATIO.  
    
******************************************************************/

proc genmod data=mi descending;
  model mi = oral age smoke / dist = binomial link = logit;
  estimate "smk log odds ratio" int 0 oral 0 age 0 smoke 1 / exp;
run;
