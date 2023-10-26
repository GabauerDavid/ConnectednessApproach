
#' @title Dynamic total connectedness plot
#' @description Visualize dynamic total connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param width The width of the graphics region in inches
#' @param height The height of the graphics region in inches
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @import graphics
#' @import grDevices
#' @export
PlotTCI = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), width=10, height=5, ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$TCI
  date = as.Date(rownames(x))
  t = length(date)
  k = ncol(x)
  lower = ylim[1]
  upper = ylim[2]
  
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/TCI.pdf"), width=width, height=height)
  par(mfrow=c(1,1), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  if (length(dim(dca$NET))>2) {
    x_ = x
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    plot(date, x_[,1], type="l", main="", las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
    grid(NA, NULL, lty=2)
    polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,1])),col=1, border=1)
    for (j in 1:dim(x_)[2]) {
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j, border=j)
    }
    legend("topleft", colnames(x_), fill=1:dim(x_)[2], bty="n")
    for (j in 1:ncol(x_)) {
      lines(date, x_[,j],col=j)
    }
    abline(h=0, lty=3)
    box()
  } else {
    if (is.null(lower)) {
      lower = 0
    }
    if (is.null(upper)) {
      upper = 100
    }
    plot(date, as.numeric(x), type="l", main="", las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
    grid(NA, NULL, lty=2)
    polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x)),col=1, border=1)
    if (!is.null(ca)) {
      for (il in 1:length(ca)) {
        lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$TCI, col=il+1)
        gTCI = ca[[il]]$gTCI
        if (!is.null(gTCI)) {
          for (ij in 1:ncol(gTCI)) {
            lines(as.Date(rownames(ca[[il]]$TCI)), gTCI[,ij], col=ij+2)
          }
        }
      }
      if (length(ca)==1) {
        if (ca[[1]]$config$approach=="Internal" || ca[[1]]$config$approach=="External") {
          legend("topleft", c("TCI", paste("TCI",ca[[1]]$config$approach), colnames(gTCI)), fill=1:(ncol(gTCI)+2), bty="n")
        } else if (ca[[1]]$config$approach=="Inclusive" || ca[[1]]$config$approach=="Exclusive") {
          legend("topleft", c("TCI", paste("TCI", ca[[1]]$config$approach)), fill=1:2, bty="n")
        }
      }
    }
    box()
  }
  if (!is.null(path)) dev.off()
}
