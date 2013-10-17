###################################################
### code chunk number 9: Cs3_code
###################################################
# Read in the data
# The data are logged already
sealData = t(harborSealnomiss)
head(sealData) #look at the first 6 lines
rownames(sealData) #look at the column names (site names)

#rename columns
years = sealData[1,] #first row is years
sealData=sealData[2:dim(sealData)[1],] #rows 2:10 are the data
# make plots of the time series (year is in first column)
par(mfrow=c(3,3))
for(i in 1:9) {
    plot(years, sealData[i,], xlab="Year", ylab="", main=names(sealData)[i])
}

# The following is an example to show you how to build and run a model
# for one of the hypotheses
# Hypothesis 5 build a simple panmictic model - 1 subpop, 1 R, 1 Q, 1 U
###############################################################
   Q.model="unconstrained" #could be left out since Q scalar
   R.model="diagonal and equal"
   U.model="unconstrained" #could be left out since U scalar
   Z.model=matrix(1,9,1)
   Z.model=factor(rep(1,9)) #or you can use factor

kem = MARSS(sealData, model=list(Z=Z.model, Q=Q.model, R=R.model, U=U.model))

# the parameter estimates; just the estimated elements
coef(kem, type="vector")
#num parameters, loglike, AIC (or use kem$AICc if you wish)
c(kem$num.params, kem$logLik, kem$AIC)

# plot the data versus predicted population states
Xpred = t(kem$states)
Xobs = sealData
A.model = coef(kem,type="matrix")$A
par(mfrow=c(3,3))
for(i in 1:9) {
    plot(years, Xobs[i,], ylab="", main=rownames(sealData)[i] )
    lines(years, Xpred[,Z.model[i]]+A.model[i], ylab="", lwd=2, col=2 )
}                                                            
mtext("Predicted (line) and Observed (points)", side=2, outer=T, line=-2)

# plot the residuals
Xpred = t(kem$states)
Xobs = sealData
par(mfrow=c(3,3))
for(i in 1:9) {
    plot(years, Xpred[,Z.model[i]] - Xobs[i,], ylab="Predicted-Observed Data", main=rownames(sealData)[i] )
}


