###################################################
### code chunk number 3: Cs_01_ar2.sim
###################################################
TT=50
true.2=c(r=0,b1=-1.5,b2=-0.75,q=1)
temp=arima.sim(n=TT,list(ar=true.2[2:3]),sd=sqrt(true.2[4]))
sim.ar2=matrix(temp,nrow=1)


###################################################
### code chunk number 4: Cs_02_ar2.model
###################################################
Z=matrix(c(1,0),1,2)
B=matrix(list("b1",1,"b2",0),2,2)
U=matrix(0,2,1)
Q=matrix(list("q",0,0,0),2,2)
A=matrix(0,1,1)
R=matrix(0,1,1)
pi=matrix(sim.ar2[2:1],2,1)
V=matrix(0,2,2)
model.list.2=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)


###################################################
### code chunk number 5: Cs_03_ar2.fit
###################################################
init.list=list(Q=matrix(1,1,1),B=matrix(1,2,1))
ar2=MARSS(sim.ar2[2:TT],model=model.list.2,inits=init.list)
print(cbind(true=true.2[2:4],estimated=coef(ar2,type="vector")))


###################################################
### code chunk number 6: Cs_04_ar2.gappy
###################################################
TT=50
gappy.data=sim.ar2[2:TT]
gappy.data[floor(runif(TT/2,1,TT))]=NA
ar2.gappy=MARSS(gappy.data,model=model.list.2,inits=init.list)


###################################################
### code chunk number 7: Cs_05_arima
###################################################
arima(gappy.data,order=c(2,0,0),include.mean=FALSE)


###################################################
### code chunk number 11: Cs_061_mar2.sim
###################################################
temp2=arima.sim(n=TT,list(ar=true.2[2:3]),sd=sqrt(true.2[4]))
sim.mar2=rbind(temp1,temp2)


###################################################
### code chunk number 9: Cs_06_mar2.sim
###################################################
TT=50
true.2=c(r=0,b1=-1.5,b2=-0.75,q=1)
temp1=arima.sim(n=TT,list(ar=true.2[2:3]),sd=sqrt(true.2[4]))


###################################################
### code chunk number 12: Cs_07_mar2.model
###################################################
Z=matrix(c(1,0,0,1,0,0,0,0),2,4)
B1=matrix(list(0),2,2); diag(B1)="b1"
B2=matrix(list(0),2,2); diag(B2)="b2"
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
U=matrix(0,4,1)
Q=matrix(list(0),4,4)
Q[1,1]="q"; Q[2,2]="q"
A=matrix(0,2,1)
R=matrix(0,2,2)
pi=matrix(c(sim.mar2[,2],sim.mar2[,1]),4,1)
V=matrix(0,4,4)
model.list.2m=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)


###################################################
### code chunk number 13: Cs_08_mar2.fit
###################################################
init.list=list(Q=matrix(1,1,1),B=matrix(1,2,1))
mar2=MARSS(sim.mar2[,2:TT],model=model.list.2m,inits=init.list)


###################################################
### code chunk number 14: Cs_09_mar2.compare
###################################################
model.list.2$x0=matrix(sim.mar2[1,2:1],2,1)
mar2a=MARSS(sim.mar2[1,2:TT],model=model.list.2,inits=init.list)
model.list.2$x0=matrix(sim.mar2[2,2:1],2,1)
mar2b=MARSS(sim.mar2[2,2:TT],model=model.list.2,inits=init.list)


###################################################
### code chunk number 15: Cs_10_compare.mars
###################################################
print(cbind(true=true.2[2:4],est.mar2=coef(mar2,type="vector"),est.mar2a=coef(mar2a,type="vector"),est.mar2b=coef(mar2b,type="vector")))


###################################################
### code chunk number 17: Cs_11_sim-ar3-data
###################################################
TT=100
temp=arima.sim(n=TT,list(ar=c(-1.5,-.75, .05)),sd=1)
sim.ar3=matrix(temp,nrow=1)


###################################################
### code chunk number 18: Cs_12_set-up-ar3-model
###################################################
Z=matrix(c(1,0,0),1,3)
B=matrix(list("b.1",1,0,"b.2",0,1,"b.3",0,0),3,3)
U=matrix(0,3,1)
Q=matrix(list(0),3,3); Q[1,1]="q.1"
A=matrix(0,1,1)
R=matrix(0,1,1)
pi=matrix(sim.ar3[3:1],3,1)
V=matrix(0,3,3)
model.list=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)


###################################################
### code chunk number 19: Cs_13_fit-ar3
###################################################
init.list=list(Q=matrix(1,1,1),B=matrix(1,3,1))
ar3=MARSS(sim.ar3[3:TT],model=model.list,inits=init.list)


###################################################
### code chunk number 21: Cs_14_fig-arss-model
###################################################
TT=1000 #set long
true.2ss=c(r=.5,b1=-1.5,b2=-0.75,q=.1)
temp=arima.sim(n=TT,list(ar=true.2ss[2:3]),sd=sqrt(true.2ss[4]))
sim.ar=matrix(temp,nrow=1)
noise=rnorm(TT-1,0,sqrt(true.2ss[1]))
noisy.data=sim.ar[2:TT]+noise

#set up a lag-2 model as MARSS lag-1
Z=matrix(c(1,0),1,2)
B=matrix(list("b1",1,"b2",0),2,2)
U=matrix(0,2,1)
Q=matrix(list("q",0,0,0),2,2)
A=matrix(0,1,1)
R=matrix("r")
V=matrix(0,2,2)
pi=matrix(mean(noisy.data),2,1)
model.list.2ss=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=0)

#give it some reasonable inits to run faster
init.list=list(Q=matrix(.01,1,1),B=matrix(1,2,1))
ar2ss=MARSS(noisy.data,model=model.list.2ss,inits=init.list)


###################################################
### code chunk number 22: Cs_15_fit-arss2-model
###################################################
model.list.2ss$R=matrix(0)
ar2ss2=MARSS(noisy.data,model=model.list.2ss,inits=init.list)
print(cbind(true=true.2ss,
  est.no.noise=c(NA,coef(ar2ss2,type="vector")),
  est.noisy=coef(ar2ss,type="vector")))


###################################################
### code chunk number 23: Cs_16_fit-arss2-with-arima
###################################################
arima(noisy.data,order=c(2,0,2),include.mean=FALSE)


###################################################
### code chunk number 25: Cs_17_code-to-compare-arss-estimation
###################################################
# This is the code used to make the figure comparing different ways to estimate 
# AR parameters from AR with noise data
# saved.res is just a flag.  Set to FALSE to run the if statement
if(!saved.res) {
file=paste("AR2SS",TT,".RData",sep="")
params=matrix(0,8,nsim)
#sim 2   true.2ss=c(r=.5,b1=0.8,b2=-0.2,q=.1)
#sim 1   
true.2ss=c(r=.5,b1=-1.5,b2=-0.75,q=.1)

Z=matrix(c(1,0),1,2)
B=matrix(list("b1",1,"b2",0),2,2)
U=matrix(0,2,1)
Q=matrix(list("q",0,0,0),2,2)
A=matrix(0,1,1)
V=matrix(0,2,2)
R=matrix("r")
pi=matrix(0,2,1) #since demeaned
model.list.2ss=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)

for(i in 1:nsim){
temp=arima.sim(n=TT,list(ar=true.2ss[2:3]),sd=sqrt(true.2ss[4]))
sim.ar=matrix(temp,nrow=1)
noise=rnorm(TT,0,sqrt(true.2ss[1]))
noisy.data=sim.ar+noise
noisy.data=as.vector(noisy.data-mean(noisy.data)) #demean
test.it=try(arima(noisy.data[2:TT],order=c(2,0,2),include.mean=FALSE))
while(class(test.it)=="try-error"){
temp=arima.sim(n=TT,list(ar=true.2ss[2:3]),sd=sqrt(true.2ss[4]))
sim.ar=matrix(temp,nrow=1)
noise=rnorm(TT,0,sqrt(true.2ss[1]))
noisy.data=sim.ar+noise
noisy.data=as.vector(noisy.data-mean(noisy.data)) #demean
test.it=try(arima(noisy.data[2:TT],order=c(2,0,2),include.mean=FALSE))
}
init.list=list(Q=matrix(.01,1,1),B=matrix(1,2,1))
tmp.kem=MARSS(noisy.data[2:TT],model=model.list.2ss,inits=init.list,silent=TRUE)
params[1:2,i]=coef(tmp.kem)$B
tmp.bfgs=MARSS(noisy.data[2:TT],model=model.list.2ss,inits=init.list,silent=TRUE,method="BFGS")
#if(any(is.na(tmp.bfgs$states.se)) | any(is.na(tmp.kem$states.se))) print(i)
params[3:4,i]=coef(tmp.bfgs)$B
params[5:6,i]=arima(noisy.data[2:TT],order=c(2,0,2),include.mean=FALSE)$coef[1:2]
params[7:8,i]=arima(noisy.data[2:TT],order=c(2,0,0),include.mean=FALSE)$coef[1:2]
cat(i);cat("\n")
if((i %% 25)==0) save(true.2ss,TT,params,file=file)
}
}
file=paste("C:\\Users\\ELI~1.HOL\\Desktop\\",file,sep="")
save(true.2ss,TT,params,file=file)


