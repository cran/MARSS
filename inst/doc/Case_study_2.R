###################################################
### code chunk number 12: Cs2_Exercise1
###################################################
#Read in data
dat=t(harborSealWA) #Transpose since MARSS needs time ACROSS columns
years = dat[1,]      
n = nrow(dat)-1
dat = dat[2:nrow(dat),]
legendnames = (unlist(dimnames(dat)[1]))

#estimate parameters
Z.model = factor(c(1,1,1,1,1))
R.model = "diagonal and unequal" 
kem1 = MARSS(dat, model=
  list(Z=Z.model, R=R.model))

#make figure
matplot(years, t(dat),xlab="",ylab="index of log abundance",
  pch=c("1","2","3","4","5"),ylim=c(5,9),bty="L")
lines(years,kem1$states-1.96*kem1$states.se,type="l",
  lwd=1,lty=2,col="red")
lines(years,kem1$states+1.96*kem1$states.se,type="l",
  lwd=1,lty=2,col="red")
lines(years,kem1$states,type="l",lwd=2)
title("Observations and total population estimate",cex.main=.9)

#show params
kem1$par
kem1$logLik
kem1$AIC


###################################################
### code chunk number 18: Cs2_Exercise2
###################################################
#fit model
Z.model = factor(c(1,1,1,1,1))
R.model = "diagonal and equal" 
kem2 = MARSS(dat, model=
  list(Z=Z.model, R=R.model))

#show parameters
kem2$par$U       #population growth rate
kem2$par$Q       #process variance
kem2$par$R[1,1]  #observation variance
kem2$logLik #log likelihood
c(kem1$AIC,kem2$AIC)

#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(kem2$par$A,
  nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
xs = matrix(kem2$states,
  nrow=dim(plotdat)[1],ncol=dim(plotdat)[2],byrow=F)
resids = plotdat-matrix.of.biases-xs
par(mfrow=c(2,3))
for(i in 1:n){
  plot(resids[!is.na(resids[,i]),i],ylab="residuals")
  title(legendnames[i])
}
par(mfrow=c(1,1))


###################################################
### code chunk number 23: Cs2_Exercise3
###################################################
#fit model
Z.model = factor(c(1,1,2,2,2))
U.model = "equal" 
Q.model = "diagonal and equal"
R.model = "diagonal and equal" 
kem3 = MARSS(dat, model=list(Z=Z.model, 
  R=R.model, U=U.model, Q=Q.model))
#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(kem3$par$A,
  nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
par(mfrow=c(2,3))
for(i in 1:n){
  j=c(1,1,2,2,2)
  xs = kem3$states[j[i],]
  resids = plotdat[,i]-matrix.of.biases[,i]-xs
  plot(resids[!is.na(resids)],ylab="residuals")
  title(legendnames[i])
}
par(mfrow=c(1,1))


###################################################
### code chunk number 24: Cs2_Exercise4
###################################################
Z.model=factor(c(1,2,3,4,5))
U.model="equal"
Q.model="diagonal and equal"
R.model="diagonal and unequal"
kem=MARSS(dat, model=list(Z=Z.model, 
  U=U.model, Q=Q.model, R=R.model) )


###################################################
### code chunk number 30: Cs2_Exercises5_7
###################################################
#Exercise 5
Z.model=factor(c(1,1,2,2,2))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and unequal"
kem = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model) )

#Exercise 6
Z.model=factor(c(1,1,1,1,2))
U.model="unequal"
Q.model="equalvarcov"
R.model="diagonal and unequal"
kem = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model) )

#Exercise 7
Z.model=factor(c(1,1,2,2,3))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and unequal"
kem = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model) )


