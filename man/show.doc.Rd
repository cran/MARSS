\name{ show.doc }
\alias{show.doc}
\title{ Documentation Utility }
\description{
  Utility in the \code{\link{MARSS-package}} to open documentation for R packages.  This is used to open files and show the directory of the \code{doc} directory of a package.   This is similar to but a bit more flexible than the R utility \code{\link{RShowDoc}}.
}
\usage{
show.doc(pkg, filename, dir="doc")
}
\arguments{
  \item{pkg}{ The name of an R package. Need not be in quotes. }
  \item{filename}{ A file name or ``manual'', ``index'' or ``dir''.  Need not be in quotes. Using ``manual'', ``index'', and ``dir'' has special behavior (no quotes required).  ``manual'' opens (tries to open) the file manual.pdf.  ``index'' tries to open the file index.html.  ``dir'' lists the files in the directory (specified by the \code{dir} argument).}
  \item{dir}{ The package subdirectory in which the file is located.  Typically documentation is 
  in the \code{doc} subdirectory and this is the default.  However, other package directories can be specified.  Use \code{dir="."} to see the directories available. }
}

\author{ 
  Eli Holmes   eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ See Also \code{\link{RShowDoc}} }
\examples{
\dontrun{
#MARSS is used as the package in these examples, but you can replace
#it with any package name you have loaded into R.  show.doc is not
#MARSS specific

#list the files in the "doc" subdirectory
show.doc(MARSS, dir)
#show the file manual.pdf in the "doc" subdirectory
show.doc(MARSS, manual) 
#show the index.html file in the "doc" subdirectory
show.doc(MARSS, index)
#Open the file Case_study_2.R in the "doc" subdirectory
show.doc(MARSS, Case_study_2.R)
#show the subdirectories in the package directory
show.doc(MARSS, dir, dir=".")
}
}

