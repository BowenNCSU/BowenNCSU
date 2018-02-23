##############################################################################
#
#   Using gls() to fit the dental data with artificial missingness
#   and unstructured covariance matrix 
#
##############################################################################

library(nlme)

#  read in the data

dent1 <- matrix(scan("dental_dropout_R.dat"),ncol=5,byrow=TRUE)
child <- dent1[,2]
age <- dent1[,3]
distance <- dent1[,4]
gender <- factor(dent1[,5])

#  create time variable for each individual as required for the CorSymm
#  correlation structure

time <- rep(seq(1,4,1),max(child))

dental.un <- gls(distance ~ -1 + gender + age:gender,correlation=corSymm(form = ~ time | factor(child)),
                weights = varIdent(form = ~ 1 |age),method="ML",na.action=na.omit)

#  This is the SAS code using proc mixed to produce the same analysis

#  data dent2; infile 'dental_dropout_sas.dat';
#    input obsno child age distance gender;
#  run;

#  proc mixed method=ml data=dent2;
#    class  child;
#    model distance = gender gender*age / noint solution chisq ;
#    repeated / type=un subject=child r=2 rcorr=2;
#  run;






