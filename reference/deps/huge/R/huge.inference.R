#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.inference(): graph inference                                     #
#-----------------------------------------------------------------------#

#' Graph inference
#'
#' Implements the inference for high dimensional graphical models, including Gaussian and Nonparanormal graphical models
#' We consider the problems of testing the presence of a single edge and the hypothesis is that the edge is absent.
#'
#' For Nonparanormal graphical model we provide Score test method and Wald Test. However it is really slow for inferencing on Nonparanormal model, especially for large data. Gaussian inference supports one variable, while Nonparanormal inference requires at least two. Nonparanormal score-test diagonal p-values do not represent edges and may be undefined; every tested off-diagonal p-value must be finite.
#'
#'
#' @param data A finite numeric \code{n} by \code{d} data matrix with at least two observations and no constant columns.
#' @param T A finite \code{d} by \code{d} estimate of the inverse correlation matrix with a positive diagonal.
#' @param adj A finite numeric or logical \code{d} by \code{d} adjacency matrix corresponding to the graph.
#' @param alpha The significance level in \code{(0, 1]}. The default value is \code{0.05}.
#' @param type The type of input data. There are 2 options: \code{"Gaussian"} and \code{"Nonparanormal"}. The default value is \code{"Gaussian"}.
#' @param method For a Nonparanormal model, the test method: \code{"score"} or \code{"wald"}. The default is \code{"score"}. Ignored for Gaussian inference.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @return
#' An object is returned:
#' \item{data}{
#'   The \code{n} by \code{d} data matrix from the input.
#' }
#' \item{p}{
#'   The \code{d} by \code{d} p-value matrix of hypothesis.
#' }
#' \item{error}{
#'   The type I error of hypothesis at alpha significance level.
#' }
#' @examples
#' #generate data
#' L = huge.generator(n = 50, d = 12, graph = "hub", g = 4)
#'
#' #graph path estimation using glasso
#' est = huge(L$data, method = "glasso")
#'
#' #inference of Gaussian graphical model at 0.05 significance level
#' T = tail(est$icov, 1)[[1]]
#' out1 = huge.inference(L$data, T, L$theta)
#'
#' #inference of Nonparanormal graphical model using score test at 0.05 significance level
#' T = tail(est$icov, 1)[[1]]
#' out2 = huge.inference(L$data, T, L$theta, type = "Nonparanormal")
#'
#' #inference of Nonparanormal graphical model using wald test at 0.05 significance level
#' T = tail(est$icov, 1)[[1]]
#' out3 = huge.inference(L$data, T, L$theta, type = "Nonparanormal", method = "wald")
#'
#' #inference of Nonparanormal graphical model using wald test at 0.1 significance level
#' T = tail(est$icov, 1)[[1]]
#' out4 = huge.inference(L$data, T, L$theta, 0.1, type = "Nonparanormal", method = "wald")
#' @references
#' 1.Q Gu, Y Cao, Y Ning, H Liu. Local and global inference for high dimensional nonparanormal graphical models.\cr
#' 2.J Jankova, S Van De Geer. Confidence intervals for high-dimensional inverse covariance estimation. \emph{Electronic Journal of Statistics}, 2015.\cr
#' @export
huge.inference = function(data, T, adj, alpha = 0.05, type = "Gaussian", method = "score"){
  if(length(type) != 1 || !is.character(type) || is.na(type) ||
     !(type %in% c("Gaussian", "Nonparanormal")))
    stop("type must be exactly one of \"Gaussian\" or \"Nonparanormal\".")
  if(type == "Nonparanormal" &&
     (length(method) != 1 || !is.character(method) || is.na(method) ||
      !(method %in% c("score", "wald"))))
    stop("method must be exactly one of \"score\" or \"wald\".")
  if(length(alpha) != 1 || !is.numeric(alpha) || !is.finite(alpha) ||
     alpha <= 0 || alpha > 1)
    stop("alpha must be a finite number in (0, 1].")

  if(!(is.matrix(data) || inherits(data, "Matrix")))
    stop("data must be a numeric matrix.")
  data.input = data
  data = as.matrix(data)
  if(!is.numeric(data))
    stop("data must be a numeric matrix.")
  d = ncol(data)
  n = nrow(data)
  if(n < 1 || d < 1)
    stop("data must be a non-empty numeric matrix.")
  if(any(!is.finite(data)))
    stop("data must contain only finite values.")
  if(n < 2)
    stop("Inference requires at least two observations.")
  if(type == "Nonparanormal" && d < 2)
    stop("Nonparanormal inference requires at least two variables.")
  constant = vapply(seq_len(d), function(j) {
    all(data[, j] == data[1, j])
  }, logical(1))
  if(any(constant))
    stop("Inference data contains a constant column.")

  if(!(is.matrix(T) || inherits(T, "Matrix")))
    stop(sprintf("T must be a numeric %d by %d matrix.", d, d))
  T = as.matrix(T)
  if(!is.numeric(T) || !identical(dim(T), c(d, d)))
    stop(sprintf("T must be a numeric %d by %d matrix.", d, d))
  if(any(!is.finite(T)))
    stop("T must contain only finite values.")

  if(!(is.matrix(adj) || inherits(adj, "Matrix")))
    stop(sprintf("adj must be a numeric or logical %d by %d matrix.", d, d))
  adj = as.matrix(adj)
  if((!is.numeric(adj) && !is.logical(adj)) ||
     !identical(dim(adj), c(d, d)))
    stop(sprintf("adj must be a numeric or logical %d by %d matrix.", d, d))
  if(any(!is.finite(adj)))
    stop("adj must contain only finite values.")
  if(any(diag(T) <= 0))
    stop("T must have a positive diagonal.")

  if(type == "Gaussian")
  {
    U = tryCatch(
      suppressWarnings(.huge_fast_cor(data)),
      error = function(e) NULL
    )
    if(is.null(U) || any(!is.finite(U)))
      stop("Gaussian inference cannot form a finite correlation matrix.")
    # De-biased estimator W = 2T - T U T with asymptotic standard deviation
    # sigma[j,k] = sqrt(T[j,j] T[k,k] + T[j,k]^2) (Jankova & van de Geer).
    variance = outer(diag(T), diag(T)) + T^2
    if(any(!is.finite(variance)) || any(variance <= 0))
      stop(paste(
        "Gaussian inference variance must be finite and positive;",
        "check the scale of T."
      ))
    W = 2*T - T%*%U%*%T
    sigma = sqrt(variance)
    p = 2*(1 - pnorm(abs(sqrt(n)*W/sigma)))
  }
  if(type == "Nonparanormal")
  {
    diag.outer = outer(diag(T), diag(T))
    if(any(!is.finite(diag.outer)) || any(diag.outer <= 0))
      stop("Products of T diagonal entries must remain finite and positive.")
    x=data
    U<-matrix(0,d,d)
    G=list()
    Test<-matrix(0,d,d)
    for(i in 1:n)
      G[[i]]<-matrix(0,d,d)
    Temp_jk<-matrix(0,n,n)
    for(j in 1:d)
    {
      for(k in 1:d)
      {
        if(j==k)
        {
          U[j, k] = 1
          next
        }

        for(i1 in 1:n)
        {
          for(i2 in 1:n)
          {
            Temp_jk[i1, i2] = sign(x[i1, j] - x[i2, j]) *
              sign(x[i1, k] - x[i2, k])
            G[[i1]][j, k] = G[[i1]][j, k] - pi/2*Temp_jk[i1, i2]

          }
        }
        U[j, k] = sin(pi/2*sum(Temp_jk)/((n-1)*n))
        for(i in 1:n)
          G[[i]][j, k] = G[[i]][j, k]/(n-1) + asin(U[j, k])
      }
    }
    #F
    F<-apply(U,1,function(x) sqrt(1-x^2))
    #R
    R<-matrix(0,d^2,d^2)
    for (i in 1:n)
      R<-R + as.matrix(as.vector(F*G[[i]]))%*%as.vector(F*G[[i]])
    R<-R/n

    #kronecker product of T
    T_k<-kronecker(T, T)

    if(method == "score")
    {
      S<-matrix(0,d,d)
      sigma<-matrix(0,d,d)
      #ST_n
      ST_n<-matrix(0,d,d)
      for(j in 1:d)
      {
        for(k in 1:d)
        {
          #S
          ej<-matrix(0,d,1)
          ek<-matrix(0,d,1)
          ej[j] = 1
          ek[k] = 1
          T_h<-T
          T_h[j, k]=0
          S[j, k] = t(ej)%*%t(T_h)%*%U%*%T_h%*%ek/(T[j, j]*T[k, k])

          idx <- (k - 1) * d + j

          #w
          temp1<-T_k[,idx]
          temp1<-temp1[-idx]
          w<-as.matrix((-temp1)/T_k[idx, idx])

          #sigma
          temp2<-R[idx,]
          temp2<-as.matrix(temp2[-idx])
          sigma[j, k] = sqrt(R[idx, idx] - 2*t(temp2)%*%w + t(w)%*%R[-idx, -idx]%*%w)

          ST_n[j, k] = S[j, k]*sqrt(n)/(2*sigma[j, k])

        }
      }

      #p-value
      p<-2*(1 - pnorm(abs(ST_n)))
      rm(temp1,temp2,ST_n,S)
    }

    if(method == "wald")
    {
      T_W<-matrix(0,d,d)
      sigma<-matrix(0,d,d)
      #W_n
      W_n<-matrix(0,d,d)
      temp1<-T%*%U
      temp2<-U%*%T
      for(j in 1:d)
      {
        for(k in 1:d)
        {
          idx <- (k - 1) * d + j

          #w
          temp3<-T_k[,idx]
          temp3<-temp3[-idx]
          w<-as.matrix((-temp3)/T_k[idx,idx])

          #sigma
          temp4<-R[idx,]
          temp4<-as.matrix(temp4[-idx])
          sigma[j, k] = sqrt(R[idx, idx] - 2*t(temp4)%*%w + t(w)%*%R[-idx, -idx]%*%w)

          #T_W
          ej<-matrix(0,d,1)
          ek<-matrix(0,d,1)
          ej[j] = 1
          ek[k] = 1
          T_W[j ,k] = (T[j,k]*t(ej)%*%temp1%*%ek + T[j, k]*t(ej)%*%temp2%*%ek - t(ej)%*%t(T)%*%temp2%*%ek)/(t(ej)%*%temp1%*%ek + t(ej)%*%temp2%*%ek - 1)
          W_n[j, k] = T_W[j, k]*sqrt(n)/(2*sigma[j, k]*T[j, j]*T[k, k])

        }
      }


      #p-value
      p<-2*(1 - pnorm(abs(W_n)))
      rm(temp1,temp2,temp3,temp4,T_W,W_n)
    }
    rm(R,F,G,w,T_k)
  }

  offdiag = row(p) != col(p)
  finite.p = if(type == "Gaussian") {
    all(is.finite(p))
  } else {
    all(is.finite(p[offdiag]))
  }
  if(!finite.p)
    stop(paste(
      "Inference produced non-finite edge p-values;",
      "the inputs are numerically degenerate."
    ))

  error=0
  for(j in 1:d)
  {
    for(k in 1:d)
    {
      if(j==k)
        next
      if(p[j, k]<alpha && adj[j, k]==0){
        error=error+1
      }
    }
  }
  error=error/d^2


  inf = list()
  inf$data = data.input
  inf$p = p
  inf$error = error

  rm(U,p)
  return(inf)
}
