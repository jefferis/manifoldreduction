#' Compute Optimal Manifold Representation of Points
#'
#' @param coords d x N matrix of d-dimensional data points
#' @param no_iterations The number of iterations to use
#' @param Verbose Whether to print status messages
#' @param maxDim The local dimensionality below which points are fixed. The
#'   default value of \code{1.2} enforces a quasi-1D manifold.
#' @param knntouse Number of nearest neighbours to consider when calculating
#'   interactions.
#' @param A spatial parameter that determines the interactions between
#'   neighbouring points according to the tradeoff F(M,Pm) = D + lambda*I (see
#'   original paper for details).
#' @param neighbourhood_size Number of nearest neighbours to consider when
#'   calculating local dimensionality (default 20).
#' @param solvemethod An integer from -1 to 5 determining the solver used to
#'   compute the local dimensionality.
#'
#' @description This implements the algorithm described in Optimal Manifold
#'   Representation of Data: An Information Theoretic Approach by Denis Chigirev
#'   and William Bialek, which attempts to reduce higher dimensional data onto a
#'   lower dimensional manifold (by default 1D).
#'
#' @return a list with the following elements: \itemize{
#'
#'   \item gamma the d x N coordinates of the output points on manifold
#'
#'   \item P
#'
#'   \item lambda determines the tradeoff F(M,Pm) = D + lambda*I
#'
#'   \item dimension the local dimensionality of the manifold
#'
#'   }
#' @importFrom RcppEigen fastLmPure
#' @export
#' @references
#' \href{http://papers.nips.cc/paper/2399-optimal-manifold-representation-of-data-an-information-theoretic-approach.pdf}{Optimal
#' Manifold Representation of Data: An Information Theoretic Approach}
manifold_reduction<-function(coords, lambda=2.0, neighbourhood_size=20L,
                             knntouse=75L, no_iterations=45L, maxDim=1.2,
                             Verbose=TRUE, solvemethod=0L){
  # % K is the number of points of the low dimensional manifold
  # % xx is the original data (should be between 0 and 1)
  # % lamba determines the tradeoff F(M,Pm) = D + lambda*I
  # % gamma will be new manifold positions
  gamma=coords;
  n=dim(coords)[1]
  K=dim(coords)[2]
  if(n>K)
    warning("dimensionality is greater than number of points.",
            " Perhaps you need to transpose your input!")

  P=rep_len(1/K,K);
  # a mask on the points in coords
  xx=P;

  # points outside this region will be ignored when calculating dimensionality
  neighbourhood_threshold=50*lambda;
  dimension=rep_len(3,K)

  log2to20=log(2:neighbourhood_size)
  moveInd=1:K

  for (z in 1:no_iterations){
    if(Verbose) message('iteration ',z,' out of ',no_iterations)
    #   % number of nearest neighbours to consider - in general the
    #   % interaction between points falls off very rapidly due to a
    #   % negative exponential.  Therefore it makes sense only to consider
    #   % a few close neighbours and set the interaction of all other
    #   % points to 0.
    #   % kpoints=min([K max([75 ceil(1.5*K^(1/3))])]);
    kpoints=pmin(K, knntouse)

    if(z>5 && K>=neighbourhood_size) {
      nnres=nabor::knn(t(gamma), k=neighbourhood_size)
      # nb original algorithm returned squared distance
      nndist=t(nnres$nn.dists[,2:neighbourhood_size])^2
      for(i in 1:K) {
        if(nndist[(neighbourhood_size-1L),i]<=neighbourhood_threshold){
          numzeros=sum(nndist[,i]==0)
          if(numzeros>0){
            #some points are right on top of this one so let's say
            dimension[i]=0
          } else {
            X=cbind(rep_len(1,(neighbourhood_size-1L)), log(sqrt(nndist[,i])))
            if(solvemethod<0) linfit=qr.solve(X, log2to20)
            else linfit=fastLmPure(X, log2to20, method = solvemethod)$coefficients
            # set dimensionality of this point to gradient
            dimension[i]=linfit[2]
          }
        } else {
          dimension[i]=0
        }
      }
      # Vectorised calculation of moveInd
      moveInd=which(dimension>maxDim & nndist[(neighbourhood_size-1L),]<=neighbourhood_threshold)
    }

    gammaNew=matrix(0, n, K)
    Pnew=rep(0,K)
    xx=xx/sum(xx)
    if(Verbose)
      message('kpoints: ',kpoints,' moveInd: ',length(moveInd))

    # find kpoints nearest neighbours from gamma for each coord
    nnres = nabor::knn(t(gamma), query=t(coords), k=kpoints)
    nndist=t(nnres$nn.dists)^2
    nnidx=t(nnres$nn.idx)

    # Precompute since it is unchanged inside the loop
    negexpdist=exp( -nndist / (2*lambda^2) );
    # Iterate over all points in current region
    for (u in 1:K){
      nnidxsForThisPoint=nnidx[,u]

      # -ve exponential of distance/space constant
      Px=P[nnidxsForThisPoint]*negexpdist[,u]
      # normalise so weight of all points is 1
      Px=Px/sum(Px)

      # If first item in Px has gone out of range
      # then zero corresponding point in mask
      if (Px[1]<0 || Px[1]>1){
        xx[u]=0
        Px=rep(0,kpoints)
      }
      # Add to Pnew the xth fraction of Px
      Pnew[nnidxsForThisPoint]=Pnew[nnidxsForThisPoint]+xx[u]*Px
      # add to every point in gammaNew a fraction of
      # the original coords of current point * Px weight
      gammaNew[,nnidxsForThisPoint]=gammaNew[,nnidxsForThisPoint]+
        outer(coords[,u]*xx[u], Px)
    }

    Pnew=Pnew+10^(-30)
    Pnew=Pnew/sum(Pnew)

    for(i1 in 1:n)
      gammaNew[i1,]=gammaNew[i1,]/Pnew

    P[moveInd]=Pnew[moveInd]
    gamma[,moveInd]=gammaNew[,moveInd]
  }
  return(list(gamma=gamma,P=P,lambda=lambda,dimension=dimension))
}
