
#' @title Dynamic net total directional connectedness plot
#' @description Visualize dynamic net total directional connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @export
PlotNET = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$NET
  date = as.Date(rownames(x))
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  lower = ylim[1]
  upper = ylim[2]

  k_row = ceiling(sqrt(k))
  k_col = ceiling(k/k_row)
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/NET.pdf"), width=10, height=7)
  par(mfcol=c(k_row,k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  if (length(dim(dca$NET))>2) {
    for (i in 1:k) {
      x_ = x[,i,]
      if (is.null(lower)) {
        lower = min(c(min(x_), min(apply(x_,1,sum))))
      }
      if (is.null(upper)) {
        upper = max(c(max(x_), max(apply(x_,1,sum))))
      }
      plot(date, x_[,1], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
      grid(NA, NULL, lty=2)
      for (j in 1:ncol(x_)) {
        polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j, border=j)
      }
      for (j in 1:ncol(x_)) {
        lines(date, x_[,j],col=j)
      }
      abline(h=0, lty=3)
      legend("topleft", colnames(x_), fill=c(1:(ncol(x_)+1)), bty="n")
      box()
    }
  } else {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (i in 1:k) {
      plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper))#, ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[,i])),col=1, border=1)
      if (!is.null(ca)) {
        for (il in 1:length(ca)) {
          lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$NET[,i], col=il+1)
        }
      }
      abline(h=0, lty=3)
      box()
    }
  }
  if (!is.null(path)) dev.off()
}
