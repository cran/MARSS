#######################################################################################################
#   MARSScheckdims function
#   Utility function to check dims in MARSS lists
#   Called by is.marssm(), is.marssMLE(), MARSScheckpar()
#######################################################################################################
MARSScheckdims <- function(el, target, n, m)
{
      if (el == "Z") { 
	flag <- !isTRUE(all.equal( dim(target$Z), c(n, m) ) ) 
      }
      if (el == "A") {
	flag <- !isTRUE(all.equal( dim(target$A), c(n, 1) ) )
      }
      if (el == "R") {
	flag <- !isTRUE(all.equal( dim(target$R), c(n, n) ) )
      }
      if (el %in% c("B", "Q", "V0")) {
	flag <- !isTRUE(all.equal( dim(target[[el]]), c(m, m) ) )
      }
      if (el %in% c("U", "x0")) {
	flag <- !isTRUE(all.equal( dim(target[[el]]), c(m, 1) ) )
      }

      flag
}
