#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.generator(): Data generator                                      #
#-----------------------------------------------------------------------#

#' Data generator
#'
#' Implements the data generation from multivariate normal distributions with different graph structures, including \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"band"} and \code{"scale-free"}.
#'
#' @param n The number of observations (sample size), as an integer of at least \code{2}. The default value is \code{200}.
#' @param d The positive integer number of variables (dimension). The default value is \code{50}.
#' @param graph The graph structure with 5 options: \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"band"} and \code{"scale-free"}.
#' @param v A finite positive number used for the off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with \code{u}. The default value is \code{0.3}.
#' @param u A finite positive number added to the diagonal elements of the precision matrix to control the magnitude of partial correlations. The default value is \code{0.1}.
#' @param g A finite positive integer. For \code{"cluster"} or \code{"hub"} graph, \code{g} is the number of clusters or hubs; values above \code{d} are treated as \code{d}. The default is about \code{d/20} if \code{d >= 40} and \code{min(2, d)} otherwise. For \code{"band"} graph, \code{g} is the bandwidth and defaults to \code{1}; values above \code{d - 1} add no edges. Not applicable to \code{"random"} or \code{"scale-free"} graph.
#' @param prob For \code{"random"} or \code{"cluster"} graph, a finite number in \code{[0, 1]} giving the probability that a pair of nodes has an edge. The default is \code{min(1, 3/d)} for \code{"random"}. For \code{"cluster"}, the default is \code{min(1, 6*g/d)} if \code{d/g <= 30} and \code{0.3} otherwise. Not applicable to \code{"hub"}, \code{"band"}, or \code{"scale-free"} graph.
#' @param vis A single non-missing logical. If \code{TRUE}, visualize the adjacency matrix, graph pattern, covariance matrix, and empirical correlation matrix. The default is \code{FALSE}.
#' @param verbose A single non-missing logical. If \code{FALSE}, tracing output is disabled. The default is \code{TRUE}.
#' @details
#' Given the adjacency matrix \code{theta}, the graph patterns are generated as below:\cr\cr
#' (I) \code{"random"}: Each pair of off-diagonal elements are randomly set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with probability \code{prob}, and \code{0} otherwise. It results in about \code{d*(d-1)*prob/2} edges in the graph.\cr\cr
#' (II)\code{"hub"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each group is associated with a "center" row \code{i} in that group. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} if \code{j} also belongs to the same group as \code{i} and \code{0} otherwise. It results in \code{d - g} edges in the graph.\cr\cr
#' (III)\code{"cluster"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with the probability \code{prob}if both \code{i} and \code{j} belong to the same group, and \code{0} otherwise. It results in about \code{g*(d/g)*(d/g-1)*prob/2} edges in the graph.\cr\cr
#' (IV)\code{"band"}: Let \code{b=min(g,d-1)}. The off-diagonal elements are set to be \code{theta[i,j]=1} if \code{1<=|i-j|<=b} and \code{0} otherwise. It results in \code{(2d-1-b)*b/2} edges in the graph.\cr\cr
#' (V) \code{"scale-free"}: The graph is generated using B-A algorithm. The initial graph has two connected nodes and each new node is connected to only one node in the existing graph with the probability proportional to the degree of the each node in the existing graph. It results in \code{d-1} edges for \code{d>=2}; for \code{d=1}, the graph is empty.
#'
#' The adjacency matrix \code{theta} has all diagonal elements equal to \code{0}. To obtain a positive definite precision matrix, the smallest eigenvalue of \code{theta*v} (denoted by \code{e}) is computed. Then we set the precision matrix equal to \code{theta*v+(|e|+0.1+u)I}. The covariance matrix is then computed to generate multivariate normal data. Every graph type returns an empty graph with zero sparsity when \code{d=1}.
#' @return
#' An object with S3 class "sim" is returned:
#' \item{data}{
#'   The \code{n} by \code{d} matrix for the generated data
#' }
#' \item{sigma}{
#'   The covariance matrix for the generated data
#' }
#' \item{omega}{
#'   The precision matrix for the generated data
#' }
#' \item{sigmahat}{
#'   The empirical correlation matrix for the generated data
#' }
#' \item{theta}{
#'   The adjacency matrix of true graph structure (in sparse matrix representation) for the generated data
#' }
#' \item{sparsity}{
#'   The proportion of possible off-diagonal entries present in the graph; zero when \code{d=1}
#' }
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}
#' @examples
#' ## band graph with bandwidth 3
#' L = huge.generator(graph = "band", g = 3)
#' plot(L)
#'
#' ## random sparse graph
#' L = huge.generator(vis = TRUE)
#'
#' ## random dense graph
#' L = huge.generator(prob = 0.5, vis = TRUE)
#'
#' ## hub graph with 6 hubs
#' L = huge.generator(graph = "hub", g = 6, vis = TRUE)
#'
#' ## hub graph with 8 clusters
#' L = huge.generator(graph = "cluster", g = 8, vis = TRUE)
#'
#' ## scale-free graphs
#' L = huge.generator(graph="scale-free", vis = TRUE)
#' @export
huge.generator = function(n = 200, d = 50, graph = "random", v = NULL, u = NULL, g = NULL, prob = NULL, vis = FALSE, verbose = TRUE){
  if(length(n) == 1 && is.numeric(n) && is.finite(n) &&
     n > .Machine$integer.max)
    stop("n exceeds R's supported integer range.")
  n = .huge_validate_positive_integer(n, "n")
  if(n < 2)
    stop("n must be at least 2 so empirical correlation is defined.")
  if(length(d) == 1 && is.numeric(d) && is.finite(d) &&
     d > .Machine$integer.max)
    stop("d exceeds R's supported integer range.")
  d = .huge_validate_positive_integer(d, "d")
  graph.types = c("random", "hub", "cluster", "band", "scale-free")
  if(length(graph) != 1 || !is.character(graph) || is.na(graph) ||
     !(graph %in% graph.types))
    stop(paste(
      "graph must be exactly one of",
      "\"random\", \"hub\", \"cluster\", \"band\", or \"scale-free\"."
    ))

  .validate_flag = function(value, name){
    if(length(value) != 1 || !is.logical(value) || is.na(value))
      stop(sprintf("%s must be TRUE or FALSE.", name))
    value
  }
  vis = .validate_flag(vis, "vis")
  verbose = .validate_flag(verbose, "verbose")

  .validate_positive_number = function(value, name){
    if(length(value) != 1 || !is.numeric(value) || !is.finite(value) ||
       value <= 0)
      stop(sprintf("%s must be a finite positive number.", name))
    as.numeric(value)
  }
  v = .validate_positive_number(if(is.null(v)) 0.3 else v, "v")
  u = .validate_positive_number(if(is.null(u)) 0.1 else u, "u")

  grouped.graph = graph == "hub" || graph == "cluster"
  if(grouped.graph || graph == "band"){
    if(is.null(g)){
      g = 1L
      if(grouped.graph){
        if(d > 40)  g = ceiling(d/20)
        if(d <= 40) g = min(2, d)
      }
    } else {
      if(length(g) == 1 && is.numeric(g) && is.finite(g) &&
         g > .Machine$integer.max)
        stop("g exceeds R's supported integer range.")
      g = .huge_validate_positive_integer(g, "g")
    }
    if(grouped.graph) g = min(g, d)
  } else {
    g = 1L
  }

  .validate_probability = function(value){
    if(length(value) != 1 || !is.numeric(value) || !is.finite(value) ||
       value < 0 || value > 1)
      stop("prob must be a finite number in [0, 1].")
    as.numeric(value)
  }
  .transform_prob = function(p) sqrt(p/2)*(p<0.5)+(1-sqrt(0.5-0.5*p))*(p>=0.5)

  if(graph == "random"){
    prob = if(is.null(prob)) min(1, 3/d) else .validate_probability(prob)
    prob = .transform_prob(prob)
  }

  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)  prob = 0.3
      if(d/g <= 30) prob = min(1, 6*g/d)
    } else {
      prob = .validate_probability(prob)
    }
    prob = .transform_prob(prob)
  }

  if(verbose) cat("Generating data from the multivariate normal distribution with the", graph,"graph structure....")

  # partition variables into groups only for graph types that use them
  if(grouped.graph){
    g.large = d%%g
    g.small = g - g.large
    n.small = floor(d/g)
    n.large = n.small+1
    g.list = c(rep(n.small,g.small),rep(n.large,g.large))
    g.ind = rep(seq_len(g),g.list)
  }

  # build the graph structure
  theta = matrix(0,d,d);
  if(graph == "band"){
    for(i in seq_len(min(g, d - 1L))){
      index = seq_len(d - i)
      theta[cbind(index, index + i)] = 1
      theta[cbind(index + i, index)] = 1
    }
  }
  if(graph == "cluster"){
    for(i in 1:g){
       tmp = which(g.ind==i)
       tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
       tmp2 = tmp2 + t(tmp2)
       theta[tmp,tmp][tmp2<prob] = 1
    }
  }
  if(graph == "hub"){
    for(i in 1:g){
       tmp = which(g.ind==i)
       theta[tmp[1],tmp] = 1
       theta[tmp,tmp[1]] = 1
    }
  }
  if(graph == "random"){
    tmp = matrix(runif(d^2,0,0.5),d,d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
  }

  if(graph == "scale-free" && d > 1L){
  out = .Call("_huge_SFGen", 2, d, PACKAGE= "huge")
  theta = matrix(as.numeric(out$G),d,d)
  }
  diag(theta) = 0
  omega = theta*v

  # make omega positive definite and standardized
  diag(omega) = abs(min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values)) + 0.1 + u
  # omega is symmetric positive definite: invert via Cholesky (chol2inv),
  # ~2x faster than solve() and symmetric by construction
  sigma = cov2cor(chol2inv(chol(omega)))
  omega = chol2inv(chol(sigma))

  # generate multivariate normal data
  x = mvrnorm(n,rep(0,d),sigma)
  sigmahat = .huge_fast_cor(x)

  # graph and covariance visulization
  if(vis == TRUE){
  old.par = .huge_graphics_state()
  on.exit(.huge_restore_graphics_state(old.par), add = TRUE)
  par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  image(theta, col = gray.colors(256),  main = "Adjacency Matrix")
  image(sigma, col = gray.colors(256), main = "Covariance Matrix")
  g = graph_from_adjacency_matrix(theta, mode="undirected", diag=FALSE)
  layout.grid = layout_with_fr(g)
  plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")
  image(sigmahat, col = gray.colors(256), main = "Empirical Correlation Matrix")
  }
  if(verbose) cat("done.\n")

  sparsity = if(d == 1L) 0 else sum(theta)/(d*(d-1))
  sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega, theta = Matrix(theta,sparse = TRUE), sparsity=sparsity, graph.type=graph)
  class(sim) = "sim"
  return(sim)
}

#' Print function for S3 class "sim"
#'
#' Print the information about the sample size, the dimension, the pattern and sparsity of the true graph structure.
#'
#' @param x An object with S3 class \code{"sim"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.generator}}
#' @export
print.sim = function(x, ...){
  cat("Simulated data generated by huge.generator()\n")
  cat("Sample size: n =", nrow(x$data), "\n")
  cat("Dimension: d =", ncol(x$data), "\n")
    cat("Graph type = ", x$graph.type, "\n")
    d = ncol(x$data)
    sparsity = if(d == 1L) 0 else sum(x$theta)/d/(d-1)
    cat("Sparsity level:", sparsity,"\n")
}

#' Plot function for S3 class "sim"
#'
#' Visualize the covariance matrix, the empirical correlation matrix, the adjacency matrix and the graph pattern of the true graph structure.
#'
#' @param x An object with S3 class \code{"sim"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.generator}} and \code{\link{huge}}
#' @export
plot.sim = function(x, ...){
     old.par = .huge_graphics_state()
     on.exit(.huge_restore_graphics_state(old.par), add = TRUE)
     par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
     image(as.matrix(x$theta), col = gray.colors(256),  main = "Adjacency Matrix")
  image(x$sigma, col = gray.colors(256), main = "Covariance Matrix")
  g = graph_from_adjacency_matrix(x$theta, mode="undirected", diag=FALSE)
  layout.grid = layout_with_fr(g)

  plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")
  image(x$sigmahat, col = gray.colors(256), main = "Empirical Correlation Matrix")
}
