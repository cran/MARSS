###################################################
### chunk number 99: Cs3_code
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
   Q.constraint="unconstrained" #could be left out since Q scalar
   R.constraint="diagonal and equal"
   U.constraint="unconstrained" #could be left out since U scalar
   Z.constraint=factor(rep(1,9)) #repeat 1 nine times

kem = MARSS(sealData, constraint=list(Z=Z.constraint, Q=Q.constraint, R=R.constraint, U=U.constraint))

# the parameter estimates
kem$par$U
diag(kem$par$Q) #since 0s are on the off-diagonals, we only need to see the diagonal
diag(kem$par$R)
#num parameters, loglike, AIC (or use kem$AICc if you wish)
c(kem$num.params, kem$logLik, kem$AIC)

# plot the data versus predicted population states
Xpred = t(kem$states)
Xobs = sealData
par(mfrow=c(3,3))
for(i in 1:9) {
    plot(years, Xobs[i,], ylab="", main=rownames(sealData)[i] )
    lines(years, Xpred[,Z.constraint[i]]+kem$par$A[i], ylab="", lwd=2, col=2 )
}                                                            
mtext("Predicted (line) and Observed (points)", side=2, outer=T, line=-2)

# plot the residuals
Xpred = t(kem$states)
Xobs = sealData
par(mfrow=c(3,3))
for(i in 1:9) {
    plot(years, Xpred[,Z.constraint[i]] - Xobs[i,], ylab="Predicted-Observed Data", main=rownames(sealData)[i] )
}


