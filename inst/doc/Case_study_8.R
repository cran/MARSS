### R code from vignette source 'Case_Study_8.xRnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)
options(width=60)
options(prompt=" ", continue=" ")
op <- par(no.readonly = TRUE)   


###################################################
### code chunk number 2: readinredddata
###################################################
head(okanaganRedds)
logRedds = log(t(okanaganRedds)[2:3,])


###################################################
### code chunk number 3: fig1
###################################################
# Code for plotting raw Okanagan redd counts
plot(okanaganRedds[,1], okanaganRedds[,2],
    xlab = "Year", ylab="Redd counts",main="", col="red")
points(okanaganRedds[,1], okanaganRedds[,3], col="blue")
legend('topleft', inset=0.1, legend=c("Aerial survey","Ground survey"), col=c("red","blue"),
pch=21)


###################################################
### code chunk number 4: reddmodel1
###################################################
Q.model="diagonal and equal"
R.model="diagonal and equal"
U.model="equal"
Z.model=factor(c(1,1)) #1 observation time series
# Fit the single state model, where the time series are assumed to be from the 
# same population. 
kem1 = MARSS(logRedds, model=list(Z = Z.model, Q = Q.model, R = R.model,
   U = U.model))


###################################################
### code chunk number 5: reddmodel2
###################################################
R.model="diagonal and unequal"
kem2 = MARSS(logRedds, model=list(Z = Z.model, Q = Q.model, R = R.model,
   U = U.model))


###################################################
### code chunk number 6: reddmodel3
###################################################
Q.constraint="diagonal and equal"
R.constraint="diagonal and equal"
U.constraint="equal"
Z.constraint=factor(c(1,2))
model3=list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem3 = MARSS(logRedds, model=model3)


###################################################
### code chunk number 7: fig2
###################################################
# Code for plotting the fit from the best model
plot(okanaganRedds[,1], logRedds[1,],
xlab = "Year", ylab="Redd counts",main="", col="red", ylim=c(0,8))
points(okanaganRedds[,1], logRedds[2,], col="blue")
lines(okanaganRedds[,1], c(kem1$states), lty=1, lwd=2)
lines(okanaganRedds[,1], c(kem1$states + 2*kem1$states.se), lty=1, lwd=1, col="grey40")
lines(okanaganRedds[,1], c(kem1$states - 2*kem1$states.se), lty=1, lwd=1, col="grey40")


###################################################
### code chunk number 8: legend
###################################################
d <- rockfish
legendnames = (unlist(dimnames(d)[2]))[2:ncol(d)]
for(i in 1:length(legendnames)) cat(paste(i,legendnames[i],"\n",sep=" "))


###################################################
### code chunk number 9: fig3
###################################################
d <- rockfish
dat = d[,2:ncol(d)] #first col is years
x = d[,1] #first col is years
n = nrow(dat) #num time series

#set up the graphical parameters to give each data a unique line, color and width
options(warn=-99)
ltys=matrix(1,nrow=n)
cols=matrix(1:4,nrow=n)
lwds=matrix(1:2,nrow=n)
pchs=matrix(as.character(c(1:n)),nrow=n)
options(warn=0)

meany = matrix(apply(dat,2,mean,na.rm=T),nrow=nrow(dat),ncol=ncol(dat),byrow=T) #take off the mean; for better plotting and to mask pattern
adj = matrix(c(1.5,0,-2,1,2,-2,-1,0,2),nrow=nrow(dat),ncol=ncol(dat),byrow=T) #just to mask pattern in data
matplot(x,dat-meany+adj,xlab="",ylab="index of log(cpue)",type="b",pch=pchs,lty=ltys,col=cols,lwd=lwds)
title("Puget Sound Total Rockfish Indices")


###################################################
### code chunk number 10: fishmodel1
###################################################
fishdat = t(rockfish[,2:dim(rockfish)[2]])
Q.model="diagonal and equal"
R.model="diagonal and unequal"
U.model="equal"
Z.model=factor(rep(1,dim(fishdat)[1]))
model1 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem1 = MARSS(fishdat, model=model1)


###################################################
### code chunk number 11: fishmodel2
###################################################
R.model="diagonal and equal"
model2 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem2 = MARSS(fishdat, model=model2)


###################################################
### code chunk number 12: fishmodel3
###################################################
R.model=matrix(list(0),9,9)
diag(R.model)=list("1","2","3","4","4","4","4","5","6")
model3 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem3 = MARSS(fishdat, model=model3)


###################################################
### code chunk number 13: getdat4
###################################################
years = rockfish[,1]
dat = rockfish[,-1]
dat4 = dat[years>=1980,]
years4=years[years>=1980]
isdata = apply(is.na(dat4),2,sum)!=length(years4)
dat4 = dat4[,isdata]


###################################################
### code chunk number 14: setargs4
###################################################
n4 = ncol(dat4)
Z4 = as.factor(rep(1,n4))
R.model4=matrix(list(0),n4,n4)
diag(R.model4)=list("1","2","2","2","2","3","4")


###################################################
### code chunk number 15: fitanalysis3
###################################################
model4 =  list(Z = Z4, Q = Q.model, R = R.model4, U = U.model)
kem4 = MARSS(t(dat4), model=model4)


###################################################
### code chunk number 16: birddat
###################################################
birddat = t(kestrel[,2:4]) 
head(kestrel)


###################################################
### code chunk number 17: fig5
###################################################
# Make a plot of the three time series
plot(kestrel[,1], kestrel[,2], xlab = "Year", ylab="Index of kestrel abundance",main="", col="red",ylim=c(0,2), pch=21)
points(kestrel[,1], kestrel[,3], col="blue", pch=22)
points(kestrel[,1], kestrel[,4], col="purple", pch=25)
legend('topright',inset=0.1, legend=c("British Columbia","Alberta","Saskatchewan"), col=c("red","blue","purple"), pch=c(21,22,25))


###################################################
### code chunk number 18: bird1
###################################################
Q.model="diagonal and equal"
R.model="diagonal and equal"
U.model="equal"
Z.model=factor(c(1,1,1))
model1=list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem1 = MARSS(birddat, model=model1, control=list(minit=100) )


###################################################
### code chunk number 19: bird1.aic
###################################################
kem1$AICc


###################################################
### code chunk number 20: bird2
###################################################
Z.model=factor(c(1,2,3)) 
model2=list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem2 = MARSS(birddat, model=model2)


###################################################
### code chunk number 21: bird2.aic
###################################################
kem2$AICc


###################################################
### code chunk number 22: bird3
###################################################
Q.model="diagonal and unequal"
model3=list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem3 = MARSS(birddat, model=model3)


###################################################
### code chunk number 23: bird3.aic
###################################################
kem3$AICc


###################################################
### code chunk number 24: bird4
###################################################
Z.model=factor(c(1,2,2)) #1 observation time series for each x time series
model4=list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem4 = MARSS(birddat, model=model4)


###################################################
### code chunk number 25: bird4.aic
###################################################
kem4$AICc


###################################################
### code chunk number 26: fig6
###################################################
# Make a plot of the predicted trajectory, confidence intervals, and the raw data in log-space
plot(kestrel[,1], kestrel[,2], xlab = "Year", ylab="Index of kestrel abundance",main="", col="red", ylim=c(0,2))
points(kestrel[,1], kestrel[,3], col="blue")
points(kestrel[,1], kestrel[,4], col="purple")
lines(kestrel[,1], c(kem4$states[1,]), lty=3, lwd=2, col="red")
lines(kestrel[,1], c(kem4$states[2,]), lty=3, lwd=2, col="blue")
lines(kestrel[,1], c(kem4$states[2,]+kem4$par$A[3,1]), lty=3, lwd=2, col="purple")
legend('topright',inset=0.1, legend=c("British Columbia","Alberta","Saskatchewan"), col=c("red","blue","purple"), pch=21)


