#' @title Dynamic pairwise connectedness plot
#' @description Visualize dynamic pairwise connectedness
#' @param dca Connectedness object
#' @param ca Compare dca object with a single connectedness object or a list of of connectedness objects
#' @param path Path where plots should be saved
#' @param ylim A vector including the lower and upper limit of the y-axis
#' @param selection Indidcator of the illustrated series
#' @param ... Arguments to be passed to methods, such as graphidcal parameters (see par).
#' @return Return connectedness plot
#' @export
PlotPCI = function(dca, ca=NULL, path=NULL, ylim=c(NULL, NULL), selection=NULL, ...) {
  message("The pairwise connectedness index is implemented according to:\n Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management, 60, 100680.")
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (length(ca)>0 && !is.null(ca$approach)) {
    ca = list(ca)
  }
  x = dca$PCI
  if (is.null(x)) {
    stop(paste(dca$approach, "has no PCI."))
  }
  date = as.Date(dimnames(x)[[3]])
  t = length(date)
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  lower = ylim[1]
  upper = ylim[2]
  
  kk = k*(k-1)/2
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (!is.null(path)) pdf(file=paste0(path, "/PCI.pdf"), width=10, height=7)
  if (is.null(selection)) {
    k_row = ceiling(sqrt(kk))
    k_col = ceiling(kk/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  } else {
    k_row = ceiling(sqrt(k))
    k_col = ceiling(k/k_row)
    par(mfcol=c(k_row, k_col), oma=c(0,0,0,0) + 0.5, mar = c(1,1,1,1) + .5, mgp=c(1, 0.4, 0))
  }
  if (dca$approach!="Frequency") {
    if (is.null(lower)) {
      lower = min(x)
    }
    if (is.null(upper)) {
      upper = max(x)
    }
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            plot(date, x[i,j,], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x[i,j,])),col=1, border=1)
            if (!is.null(ca)) {
              for (il in 1:length(ca)) {
                lines(as.Date(rownames(ca[[il]]$TCI)), ca[[il]]$PCI[i,j,], col=il+1)
              }
            }
            abline(h=0, lty=3)
            box()
          }
        }
      }
    }
  } else {
    for (j in 1:k) {
      for (i in 1:k) {
        if (i>j) {
          if (i==selection || j==selection || is.null(selection)) {
            x_ = x[i,j,,]
            x[which(x>=99.99999,arr.ind=T)] = 0
            if (is.null(lower)) {
              lower = min(x)
            }
            if (is.null(upper)) {
              upper = max(x)
            }
            plot(date, x_[,1], type="l", main=paste(NAMES[j],"-",NAMES[i]), las=1, xlab="", ylab="", xaxs="i", yaxs="i", tck=-0.02, ylim=c(lower,upper), ...)
            grid(NA, NULL, lty=2)
            polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,1])),col=1, border=1)
            for (l in ncol(x_):1) {
              polygon(c(date,rev(date)),c(c(rep(0,t)),rev(x_[,l])),col=l, border=l)
            }
            for (l in 1:ncol(x_)) {
              lines(date, x_[,l],col=l)
            }
            abline(h=0, lty=3)
            legend("topleft", colnames(x_), fill=1:ncol(x_), bty="n")
            box()
          }
        }
      }
    }
  }
  if (!is.null(path)) dev.off()
}
