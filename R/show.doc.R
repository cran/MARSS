show.doc = function(pkg, filename, dir="doc") {
filename=as.character(substitute(filename))
pkg=as.character(substitute(pkg))
demo.path =.find.package(pkg, NULL, verbose = FALSE)
if(filename=="manual") {
  demo.path =file.path(demo.path, dir)
  demo.path =file.path(demo.path, paste(filename, ".pdf", sep=""))
  browseURL(paste("file://",demo.path, sep=""))
  return()
}
if(filename=="dir") {
  demo.path =file.path(demo.path, dir)
  return(dir(demo.path))
}
if(filename=="index") {
  demo.path =file.path(demo.path, dir)
  demo.path =file.path(demo.path, paste(filename, ".html", sep=""))
  browseURL(paste("file://",demo.path, sep=""))
  return()
}
demo.path =file.path(demo.path, dir)
demo.path =file.path(demo.path, filename)
tmp = unlist(strsplit(filename, "\\."))
if(length(tmp)==1) extension=""
else extension=tmp[length(tmp)]
if(extension!="pdf") file.show(demo.path, title=filename)
else browseURL(paste("file://", demo.path, sep=""))
}

