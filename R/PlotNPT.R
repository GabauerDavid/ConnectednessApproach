#' @title Dynamic net pairwise transmission plot
#' @description Visualize dynamic net total directional connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param width The width of the graphics region in inches
#' @param height The height of the graphics region in inches
#' @param ... Arguments to be passed to methods, such as graphidcal parameters (see par).
#' @return Return connectedness plot
#' @export
PlotNPT = function(dca, ca=NULL, path=NULL, width=10, height=7, ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$config$approach)) {
    ca = list(ca)
  }
  x = dca$NPT
  if (is.null(x)) {
    stop(paste(dca$config$approach, "has no NPT."))
  }
  date = as.Date(rownames(x))
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  k_row = ceiling(sqrt(k))
  k_col = ceiling(k/k_row)

  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/NPT.pdf"), width=width, height=height)
  par(mfcol=c(k_row,k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  if (length(dim(dca$NET))>2) {
    for (i in 1:k) {
      x_ = x[,i,]
      plot(date, apply(x_,1,sum), type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(0,k-1))#, ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(apply(x_,1,sum))),col=1, border=1)
      for (j in ncol(x_):1) {
        polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,j])),col=j+1, border=j+1)
      }
      for (j in 1:ncol(x_)) {
        lines(date, x_[,j],col=j+1)
      }
      lines(date, apply(x_,1,sum), col=1)
      abline(h=0, lty=3)
      legend("topleft", colnames(x_), fill=1:(ncol(x_)+1), bty="n")
      box()
    }
  } else {
    for (i in 1:k) {
      plot(date, x[,i], type="l", main=NAMES[i], las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(0,k-1), ...)
      grid(NA, NULL, lty=2)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[,i])),col=1, border=1)
      if (!is.null(ca)) {
        for (il in 1:length(ca)) {
          lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$NPT[,i], col=il+1)
        }
      }
      abline(h=0, lty=3)
      box()
    }
  }
  if (!is.null(path)) dev.off()
}

