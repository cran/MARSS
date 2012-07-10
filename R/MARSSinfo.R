MARSSinfo = function(number){
if(number==21)
cat("In a variance-covariance matrix, you cannot have 0s on the diagonal. When you pass
in a qxq variance-covariance matrix (Q,R or V0) with a 0 on the diagonal, MARSS rewrites
this into H%*%V.  V is a pxp variance-covariance matrix made from only the non-zero row/cols
in your variance-covariance matrix and H is a q x p matrix with all zero rows corresponding
to the 0 diagonals in the variance-covariance matrix that you passed in. By setting, the start 
variance to 0, you have forced a 0 on the diagonal of a variance-covariance matrix and that
will break the EM algorithm (the BFGS algorithm will probably work since it uses start in a different way).
This is probably just an accident in how you specified your start values for your variance-covariance matrices.
Check the start values and compare to the model$free values.  If you did not pass in start, then MARSS's function
for generating reasonable start values does not work for your model.  So just pass in your own start values
for the variances. Note, matrices with a second dim equal to 0 are fine in test$start (and test$par).
It just means the parameter is fixed (not estimated).\n")

if(number==5)
cat("If you got an error from is.marssMLE related to your model, then the first thing to do is look at the list that you
passed into the model argument.  Often that will reveal the problem.  If not, then look at your data and make
sure it is a nxT matrix (and not a Txn) and doesn't have any weird values in it.  If the problem is still not clear 
you need to look at the model that MARSS thinks you are trying to fit.
If you used test=MARSS(foo), then test is the MLE object.  If the function exited 
without giving you the MLE object, try test=MARSS(...,fit=FALSE) to get it.  Type summary(test$model) to see
a print out of the model.  If your model is time-varying, this will be very verbose so you'll want to divert
the output to a file.  Then try this test$par=test$start, now you have filled in the par element of the MLE object.
Try parmat(test,t=1) to see all the parameters at t=1 using the start as the par values.  This might reveal
the problem too.  Note, matrices with a second dim equal to 0 are fine in test$start (and test$par).
It just means the parameter is fixed (not estimated).
\n")

if(number=="denom not invertible")
cat("First check your data and covariates (if you have them).  Make sure you didn't make a mistake when entering
the data.  For example, a row of all NAs or two covariates are the same.  Then look at your model by
passing in fit=FALSE to the MARSS() call.  Are you trying to estimate B but you set Q to zero?  That won't work.  Note if you are estimating D, your error will say problems in A update.  If you are estimating C, your error will say problems in U update.  This is because in the MARSS algorithms, the models with D and C are rewritten into a simpler MARSS model with time-varying A and U.  If you have set R=0, you might get this error if you are trying to estimate A.  Did you set a VO (say, diagonal) that is inconsisent with MARSSkf(fit)$V0T (the covariance matrix implied by the model)?  That can cause problems with the Q update.
\n")

}
