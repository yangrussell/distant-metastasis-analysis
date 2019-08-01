# Performs univariate cox regression on a specific gene
#
# Parameters:
# x - a dataset used for univariate cox regression
# time - the DMFS_Time column of the info dataset
# status - the DMFS_Status column of the info dataset
# nn - a sequence from 1 to the length of time
#
# Example calls:
# unisurvcns(train,time,status,nn)

unisurvcns = function(x,time,status,nn) {
library(survival)

    time=time[nn]
    status=status[nn]
    x=x[,nn,drop=FALSE]
    y=Surv(time,status)
    xd=dim(x)
    maxNA=xd[2]/2 # Max number of NAs per gene
    unicox=matrix(NA,nrow=xd[1],ncol=7)
    for (i in 1:xd[1]) {
       if (i%%1000==0) {
           print(i)
       }
       z=x[i,]
       if (sum(is.na(z))<maxNA) {
           coxreg=try(coxph(y~x[i,]))
           if (class(coxreg) == "try-error") next;
           screg=summary(coxreg)
           unicox[i,1:5]=as.numeric(screg$coef)
           unicox[i,6]=screg$conf.int[3]
           unicox[i,7]=screg$conf.int[4]
       }
    }
    unicox=as.data.frame(unicox)
    names(unicox)=c("coef","HR","secoef","zscore","pval","HR_low","HR_high")

    return(unicox)
}

# Same function as above, but with age as covariate
unisurvcns_agecovar = function(x,time,status,age,nn) {
    time=time[nn] # Select specified subset
    status=status[nn]
    x=x[,nn,drop=FALSE]
    age=age[nn]
    y=Surv(time,status)
    xd=dim(x)
    maxNA=xd[2]/2 # Max number of NAs per gene
    unicox=matrix(NA,nrow=xd[1],ncol=15)
    for (i in 1:xd[1]) {
       if (i%%1000==0) {
           print(i)
       }
       z=x[i,]
       if (sum(is.na(z))<maxNA) {
           coxreg=try(coxph(y~x[i,]+age))
           if (class(coxreg) == "try-error") next;
           screg=summary(coxreg)
           unicox[i,1:5]=as.numeric(screg$coef[1,])
           unicox[i,6]=screg$conf.int[1,3]
           unicox[i,7]=screg$conf.int[1,4]
           unicox[i,8:12]=as.numeric(screg$coef[2,])
           unicox[i,13]=screg$conf.int[2,3]
           unicox[i,14]=screg$conf.int[2,4]
           unicox[i,15] = screg$sctest[3]
       }
    }
    unicox=as.data.frame(unicox)
    names(unicox)=c("coef","HR","secoef","zscore","pval","HR_low","HR_high","coef_AGE","HR_AGE","secoef_AGE","zscore_AGE","pval_AGE","HR_AGE","HR_high_AGE","overall_p_sctest")

    return(unicox)
}


# Same function as above, but with specified covariate as covariate
unisurvcns_withcovar = function(x,time,status,covar,nn) {
    time=time[nn] # Select specified subset
    status=status[nn]
    x=as.matrix(x[,nn])
    covar=covar[nn]
    y=Surv(time,status)
    xd=dim(x)
    maxNA=xd[2]/2 # Max number of NAs per gene
    unicox=matrix(NA,nrow=xd[1],ncol=15)
    for (i in 1:xd[1]) {
       if (i%%1000==0) {
           print(i)
       }
       z=x[i,]
       if (sum(is.na(z))<maxNA) {
           coxreg=try(coxph(y~x[i,]+covar))
           if (class(coxreg) == "try-error") next;
           screg=summary(coxreg)
           unicox[i,1:5]=as.numeric(screg$coef[1,])
           unicox[i,6]=screg$conf.int[1,3]
           unicox[i,7]=screg$conf.int[1,4]
           unicox[i,8:12]=as.numeric(screg$coef[2,])
           unicox[i,13]=screg$conf.int[2,3]
           unicox[i,14]=screg$conf.int[2,4]
           unicox[i,15] = screg$sctest[3]
       }
    }
    unicox=as.data.frame(unicox)
    names(unicox)=c("coef","HR","secoef","zscore","pval","HR_low","HR_high","coef_COVAR","HR_COVAR","secoef_COVAR","zscore_COVAR","pval_COVAR","HR_COVAR","HR_high_COVAR","overall_p_sctest")

    return(unicox)
}
