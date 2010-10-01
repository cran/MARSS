\name{allowed}
\alias{allowed}
\alias{allowed.methods}
\alias{describe.marssm}
\alias{kem.methods}
\alias{optim.methods}
\alias{alldefaults}
\alias{model.elem}
\alias{model.elem.w.V0}

\title{ MARSS function defaults and allowed methods }
\description{
  Defaults and allowed fitting methods for the \code{\link{MARSS}} function are specified in the file \code{MARSSsettings.R}. These are hidden thus to see them preface with \code{MARSS:::}.
}
\details{
\code{allowed} is a list with the allowed constraints for each fitting method (used in \code{\link{checkPopWrap}}.  \code{allowed.methods} is a vector with the allowed \code{method} arguments for the \code{\link{MARSS}} function.  \code{kem.methods} and \code{optim.methods} are vectors of \code{method} arguments that fall in each of these two categories; used by \code{\link{MARSS}}.  \code{alldefaults} is a list that specifies the defaults for \code{\link{MARSS}} arguments if they are not passed in; used by \code{\link{checkPopWrap}}.  \code{model.elem} and \code{model.elem.w.V0} specify the parameters names used in lists such are par lists.
}


