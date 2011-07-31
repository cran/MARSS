MARSSresids = function(MLEobj){
#Reference page 9 in Messy Time Series
#By definition there are no residuals for t=1
TT = dim(MLEobj$model$data)[2]
m = dim(MLEobj$par$Q)[1]
n = dim(MLEobj$par$R)[1]
rt = matrix(0,m,TT)
ut = matrix(0,n,TT)
et = st.et = matrix(0,n+m,TT)
var.et = array(0,dim=c(n+m,n+m,TT))
Nt = array(0,dim=c(m,m,TT))
Mt = array(0,dim=c(n,n,TT))
Jt = matrix(0,m,TT)
G = matrix(0,n,n+m)
if(!any(takediag(MLEobj$par$R)==0)) G[,1:n] = t(chol(MLEobj$par$R))
H = matrix(0,m,n+m)
if(!any(takediag(MLEobj$par$Q)==0)) H[,(n+1):(n+m)] = t(chol(MLEobj$par$Q))
vt = MLEobj$kf$Innov
Ft = MLEobj$kf$Sigma
Tt = MLEobj$par$B
Zt = MLEobj$par$Z
for(t in seq(TT,2,-1)){
Ftinv = solve(Ft[,,t])
Kt = Tt%*%matrix(MLEobj$kf$Kt[,,t],m,n) #R is dropping the dims so we force it to be nxm
Lt = Tt - Kt%*%Zt
Jt = H-Kt%*%G  
ut[,t] = Ftinv%*%vt[,t,drop=FALSE]-t(Kt)%*%rt[,t,drop=FALSE]
rt[,t-1] = t(Zt)%*%ut[,t,drop=FALSE]+t(Tt)%*%rt[,t,drop=FALSE]
Mt[,,t] = Ftinv + t(Kt)%*%Nt[,,t]%*%Kt
Nt[,,t-1] = t(Zt)%*%Ftinv%*%Zt + t(Lt)%*%Nt[,,t]%*%Lt
et[,t] = t(G)%*%ut[,t,drop=FALSE] + t(H)%*%rt[,t,drop=FALSE]
var.et[,,t] = t(G)%*%Ftinv%*%G + t(Jt)%*%Nt[,,t]%*%Jt
st.er.et = matrix(0,n+m,n+m)
if(!any(takediag(var.et[1:n,1:n,t])==0)) st.er.et[1:n,1:n] = solve(chol(var.et[1:n,1:n,t]))
if(!any(takediag(var.et[(n+1):(n+m),(n+1):(n+m),t])==0)) st.er.et[(n+1):(n+m),(n+1):(n+m)] = solve(chol(var.et[(n+1):(n+m),(n+1):(n+m),t]))
st.et[,t] = st.er.et%*%et[,t,drop=FALSE]
}
return(list(et=et, std.et=st.et, var.et=var.et))
}