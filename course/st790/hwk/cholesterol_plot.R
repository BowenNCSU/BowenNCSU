#########################################################################
#
#   Cholesterol data - high-dose chenodiol vs. placebo
#
#   1 = high-dose chenodiol, 2 = placebo   
#   
#   Make spaghetti plots
#
#########################################################################

#  Read in the data -- they are in the form of one record per individual
#  aka the "wide" format

thedat.wide <- read.table("cholesterol.dat",row.names=NULL)
colnames(thedat.wide) <- c("trt","id","month0","month6","month12","month20","month24")

#  Total number of individuals

m <- nrow(thedat.wide)

#  Reconfigure as one observation per record ("long" format) and sort
#  in id order

thedat.long <- reshape(thedat.wide,varying=c("month0","month6","month12","month20","month24"),
            v.names="chol",idvar="id",times=c(0,6,12,20,24),timevar="month",direction="long")

thedat <- thedat.long[order(thedat.long$id),]

#########################################################################
#
#   Make separate plots by treatment group
#
#########################################################################

library(ggplot2)

#   Make the treatment and id variables factors for use with ggplot 

thedat <- within(thedat, {  id <- factor(id)
    trt <- factor(trt,levels=1:2,labels=c("High-Dose Chenodiol","Placebo"))
})

#   First create the basic plot object

pp <- ggplot(thedat,aes(x=month,y=chol,group=id)) +
    theme(panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray94"),
    panel.grid.minor = element_line(colour = "gray94"),
    panel.spacing=unit(1, "lines")) 

#  Now customize

pdf("cholesterol_bytrt.pdf",width=8)
qq <- pp +  geom_line(color='dark gray') + geom_point(shape=18,color='dark gray') +
#    geom_smooth(aes(group = 1),method="loess",color="black",se=FALSE)
     stat_summary(aes(group = 1),geom = "line", fun.y = mean, size = 1.5) +                       
     facet_wrap(~ trt, ncol=2) + xlab("Month") + ylab("Serum Cholesterol (mg/dL")
    print(qq)
dev.off()
