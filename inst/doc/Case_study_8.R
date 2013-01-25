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
    xlab = "Year", ylab="Redd counts",main="", col="red", pch=1)
points(okanaganRedds[,1], okanaganRedds[,3], col="blue", pch=2)
legend('topleft', inset=0.1, legend=c("Aerial survey","Ground survey"),
 col=c("red","blue"), pch=c(1,2))


###################################################
### code chunk number 4: reddmodel1
###################################################
model1=list()
model1$R="diagonal and equal"
model1$Z=matrix(1,2,1) #matrix of 2 rows and 1 column
model1$A="scaling" #the default
# Fit the single state model, where the time series are assumed 
# to be from thesame population. 
kem1 = MARSS(logRedds, model=model1)


###################################################
### code chunk number 5: reddmodel2
###################################################
model2=model1 #model2 is based on model1
model2$R="diagonal and unequal"
kem2 = MARSS(logRedds, model=model2)


###################################################
### code chunk number 6: reddmodel3
###################################################
model3=list()
model3$Q="diagonal and equal"
model3$R="diagonal and equal"
model3$U="equal"
model3$Z="identity"
model3$A="zero"
kem3 = MARSS(logRedds, model=model3)


###################################################
### code chunk number 7: fig2
###################################################
# Code for plotting the fit from the best model
plot(okanaganRedds[,1], logRedds[1,],
xlab = "Year", ylab="Redd counts",main="", col="red", ylim=c(0,8))
points(okanaganRedds[,1], logRedds[2,], col="blue", pch=2)
lines(okanaganRedds[,1], c(kem1$states), lty=1, lwd=2)
lines(okanaganRedds[,1], c(kem1$states + 2*kem1$states.se), lty=1, lwd=1, col="grey40")
lines(okanaganRedds[,1], c(kem1$states - 2*kem1$states.se), lty=1, lwd=1, col="grey40")


###################################################
### code chunk number 8: set-up-rockfish-data
###################################################
rec.names=paste("Rec..targeting.bottomfish.",1:4,sep="")
rec.years=apply(!is.na(rockfish[,rec.names]),1,any)
recdat = rockfish[rec.years,rec.names]
flatrecdat=apply(recdat,1,sum,na.rm=TRUE)


###################################################
### code chunk number 9: rockfish1
###################################################
matplot(rockfish[rec.years,1],recdat,ylab="",xlab="Rec CPUE")
title("Puget Sound Total Rockfish Recreational CPUE data")


###################################################
### code chunk number 10: fit-rockfish1
###################################################
kem.rock1=MARSS(flatrecdat,model=list(Q=matrix(0.01)))


###################################################
### code chunk number 11: set-up-rock-A
###################################################
A.model=array(list(0),dim=c(1,1,dim(recdat)[1]))
manage.A=!is.na(recdat[,"Rec..targeting.bottomfish.2"])
A.model[1,1,manage.A]="manag.A"
manage.B=!is.na(recdat[,"Rec..targeting.bottomfish.3"])
A.model[1,1,manage.B]="manag.B"
manage.C=!is.na(recdat[,"Rec..targeting.bottomfish.4"])
A.model[1,1,manage.C]="manag.C"


###################################################
### code chunk number 12: fit-rockfish2
###################################################
kem.rock2=MARSS(flatrecdat,model=list(A=A.model,Q=matrix(0.01)))


###################################################
### code chunk number 13: rockfish2
###################################################
matplot(rockfish[rec.years,1],recdat,ylab="",xlab="Rec CPUE")
lines(rockfish[rec.years,1],as.vector(kem.rock1$states),col="red",lty=2,lwd=2)
lines(rockfish[rec.years,1],as.vector(kem.rock2$states),col="blue",lty=1,lwd=2)
title("Puget Sound Total Rockfish Recreational CPUE data")


###################################################
### code chunk number 14: test-rock (eval = FALSE)
###################################################
##  kem.rock3a=MARSS(flatrecdat,model=list(A=A.model,tinitx=0))
##  kem.rock3b=MARSS(flatrecdat,model=list(A=A.model,tinitx=1))


###################################################
### code chunk number 15: legend
###################################################
d <- rockfish
legendnames = (unlist(dimnames(d)[2]))[2:ncol(d)]
for(i in 1:length(legendnames)) cat(paste(i,legendnames[i],"\n",sep=" "))


###################################################
### code chunk number 16: fig3
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
### code chunk number 17: fishmodel1
###################################################
fishdat = t(rockfish[,2:dim(rockfish)[2]])
Q.model="diagonal and equal"
R.model="diagonal and unequal"
U.model="equal"
Z.model=factor(rep(1,dim(fishdat)[1]))
model1 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem1 = MARSS(fishdat, model=model1)


###################################################
### code chunk number 18: fishmodel2
###################################################
R.model="diagonal and equal"
model2 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem2 = MARSS(fishdat, model=model2)


###################################################
### code chunk number 19: fishmodel3
###################################################
R.model=matrix(list(0),9,9)
diag(R.model)=list("1","2","4","4","4","4","4","5","6")
model3 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem3 = MARSS(fishdat, model=model3)


###################################################
### code chunk number 20: testing
###################################################
U.model=array(list(),dim=c(1,1,length(rockfish[,"Year"])))
U.model[1,1,rockfish[,"Year"]<1980]="pre-1980"
U.model[1,1,rockfish[,"Year"]>=1980]="post-1980"
model4 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem4 = MARSS(fishdat, model=model4)

Z.model=factor(c("trawl","trawl","rec","rec","rec","rec","rec","scuba","wdfw.trawl"))
Q.model="diagonal and unequal"
U.model="equal"
model5 = list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem5 = MARSS(fishdat, model=model5)


###################################################
### code chunk number 21: getdat4
###################################################
years = rockfish[,1]
dat = rockfish[,-1]
dat4 = dat[years>=1980,]
years4=years[years>=1980]
isdata = apply(is.na(dat4),2,sum)!=length(years4)
dat4 = dat4[,isdata]


###################################################
### code chunk number 22: setargs4
###################################################
n4 = ncol(dat4)
Z4 = as.factor(rep(1,n4))
R.model4=matrix(list(0),n4,n4)
diag(R.model4)=list("1","2","2","2","2","3","4")


###################################################
### code chunk number 23: fitanalysis3
###################################################
model4 =  list(Z = Z4, Q = Q.model, R = R.model4, U = U.model)
kem4 = MARSS(t(dat4), model=model4)


###################################################
### code chunk number 24: birddat
###################################################
birddat = t(kestrel[,2:4]) 
head(kestrel)


###################################################
### code chunk number 25: fig5
###################################################
# Make a plot of the three time series
plot(kestrel[,1], kestrel[,2], xlab = "Year", ylab="Index of kestrel abundance",main="", col="red",ylim=c(0,2), pch=21)
points(kestrel[,1], kestrel[,3], col="blue", pch=22)
points(kestrel[,1], kestrel[,4], col="purple", pch=25)
legend('topright',inset=0.1, legend=c("British Columbia","Alberta","Saskatchewan"), col=c("red","blue","purple"), pch=c(21,22,25))


###################################################
### code chunk number 26: bird1
###################################################
model.b1=list()
model.b1$R="diagonal and equal"
model.b1$Z=matrix(1,3,1)
kem.b1 = MARSS(birddat, model=model.b1, control=list(minit=100) )


###################################################
### code chunk number 27: bird1.aic
###################################################
kem.b1$AICc


###################################################
### code chunk number 28: bird2
###################################################
model.b2=list()
model.b2$Q="diagonal and equal"
model.b2$R="diagonal and equal"
model.b2$Z="identity"
model.b2$A="zero"
model.b2$U="equal"
kem.b2 = MARSS(birddat, model=model.b2)


###################################################
### code chunk number 29: bird2.aic
###################################################
kem.b2$AICc


###################################################
### code chunk number 30: bird3
###################################################
model.b3=model.b2 #is is based on model.b2
#all we change is the structure of Q
model.b3$Q="diagonal and unequal"
model.b3$U="unequal"
kem.b3 = MARSS(birddat, model=model.b3)


###################################################
### code chunk number 31: bird3.aic
###################################################
kem.b3$AICc


###################################################
### code chunk number 32: bird4
###################################################
model.b4=list()
model.b4$Q="diagonal and unequal"
model.b4$R="diagonal and equal"
model.b4$Z=factor(c("BC","AB-SK","AB-SK"))
model.b4$A="scaling"
model.b4$U="unequal"
kem.b4 = MARSS(birddat, model=model.b4)


###################################################
### code chunk number 33: bird4.aic
###################################################
kem.b4$AICc


###################################################
### code chunk number 34: fig6
###################################################
# Make a plot of the predicted trajectory, confidence intervals, and the raw data in log-space
plot(kestrel[,1], kestrel[,2], xlab = "Year", ylab="Index of kestrel abundance",main="", col="red", ylim=c(0,2), pch=21)
points(kestrel[,1], kestrel[,3], col="blue", pch=22)
points(kestrel[,1], kestrel[,4], col="purple", pch=25)
lines(kestrel[,1], c(kem.b4$states[1,]), lty=3, lwd=2, col="red")
lines(kestrel[,1], c(kem.b4$states[2,]), lty=3, lwd=2, col="blue")
lines(kestrel[,1], c(kem.b4$states[2,]+coef(kem.b4,type="matrix")$A[3,1]), lty=3, lwd=2, col="purple")
legend('topright',inset=0.1, legend=c("British Columbia","Alberta","Saskatchewan"), col=c("red","blue","purple"), pch=c(21,22,25))


