#' @title Network plot
#' @description Visualize net pairwise or pairwise connectedness measures
#' @param x NPDC or PCI matrix
#' @param path Path where plots should be saved
#' @param method Either visualizing NPDC or PCI
#' @param name_length Length of variable names in the network plot
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @return Return connectedness plot
#' @export
#' @import igraph
PlotNetwork = function(dca, method="NPDC", path=NULL, name_length=NULL, ...) {
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
  if (method=="NPDC") {
    x = dca$NPDC
  } else if (method=="PCI") {
    x = dca$PCI
  } else {
    stop("This method does not exists")
  }
  date = as.Date(dimnames(x)[[3]])
  t = length(date)
  k = ncol(x)
  NAMES = dimnames(x)[[1]]
  if (is.null(NAMES)) {
    NAMES = 1:k
  } else {
    NAMES = colnames(x)
    if (!is.null(name_length)) {
      NAMES = substr(NAMES, 1, name_length)
    }
  }
  
  oldpar = par(no.readonly=TRUE)
  on.exit(par(oldpar)) 
  if (length(dim(x))==4) {
    kk = dim(x)[4]
    k1 = ceiling(sqrt(kk))
    k2 = ceiling(kk/k1)
  } else {
    kk = k1 = k2 = 1
    x = array(x, c(k,k,t,1))
  }
  par(mfrow = c(k1,k2), oma = c(0,0,0,0), mar = c(0,0,0,0), mgp = c(0, 0, 0))
  if (!is.null(path)) pdf(file=paste0(path, "/NetworkPlot.pdf"), width=10, height=10)
    for (ijk in 1:kk) {
      x_ = t(apply(x[,,,ijk], 1:2, mean))
      x_ = ifelse(x_<0, 0, x_)
      colnames(x_) = rownames(x_) = NAMES
      diag(x_) = 0
      x_ = x_ - min(x_)
      x_ = x_ / max(x_)
      m = 10 * x_
      if (nrow(which(m!=0, arr.ind=TRUE))>(k*(k-1)/2)) {
        gr = graph.adjacency(m, mode="undirected", weighted=TRUE)
        lo = layout_in_circle(gr)
        net = graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)
      } else {
        gr = graph.adjacency(m, mode="undirected", weighted=TRUE)
        lo = layout_in_circle(gr)
        net = graph.adjacency(m, mode="directed", weighted=TRUE, diag=FALSE)
      }
      plot.igraph(net,vertex.label=V(net)$name, layout=lo, vertex.label.cex=0.99, vertex.size=20,
                  vertex.label.color="white", edge.color="black", mark.col="steelblue4",
                  edge.width=E(net)$weight, edge.arrow.size=1)
    }
  if (!is.null(path)) dev.off()
}
