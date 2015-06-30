#' Compute Optimal Manifold Representation of Points
#'
#' @param xcoords d x N matrix of d-dimensional data points
#' @param no_iterations The number of iterations to use
#' @param Verbose Whether to print status messages
#' @param solvemethod An integer from -1 to 5 determining the solver used to
#'   compute the local dimensionality.
#' @param maxDim The local dimensionality below which points are fixed. The
#'   default value of \code{1.2} enforces a quasi-1D manifold.
#' @param knntouse Number of nearest neighbours to consider when calculating
#'   interactions.
#'
#' @description This implements the algorithm described in Optimal Manifold
#'   Representation of Data: An Information Theoretic Approach Denis Chigirev
#'   and William Bialek which attempts to reduce higher dimensional data onto a
#'   lower dimensional manifold (by default 1D).
#'
#' @return a list with the following elements: \itemize{
#'
#'   \item gamma the d x N coordinates of the output points on manifold
#'
#'   \item P
#'
#'   \item lamba determines the tradeoff F(M,Pm) = D + lambda*I
#'
#'   \item dimension the local dimensionality of the manifold
#'
#'   }
#' @importFrom RcppEigen fastLmPure
#' @references
#' \href{http://papers.nips.cc/paper/2399-optimal-manifold-representation-of-data-an-information-theoretic-approach.pdf}{Optimal
#' Manifold Representation of Data: An Information Theoretic Approach}
manifold_reduction<-function(xcoords, no_iterations=45L, Verbose=TRUE,
                             solvemethod=0L, maxDim=1.2, knntouse=75L){
  # % K is the number of points of the low dimensional manifold
  # % xx is the original data (should be between 0 and 1)
  # % lamba determines the tradeoff F(M,Pm) = D + lambda*I
  # % gamma will be new manifold positions
  gamma=xcoords;
  n=dim(xcoords)[1]
  K=dim(xcoords)[2]

  P=rep_len(1/K,K);
  # a mask on the points in xcoords
  xx=P;

  lambda=2;
  dimension=rep_len(3,K)

  maxDim=1.2; # points with local dimensionality < maxDim will be fixed

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

    if(z>5 && K>=20) {
      nnres=nabor::knn(t(gamma), k=20)
      # nb original algorithm returned squared distance
      nndist=t(nnres$nn.dists)^2
      log2to20=log(2:20)
      for(i in 1:K) {
        if(nndist[20,i]<=100){
          numzeros=sum(nndist[2:20,i]==0)
          if(numzeros>0){
            #some points are right on top of this one so let's say
            dimension[i]=0
          } else {
            X=cbind(rep_len(1,19), log(sqrt(nndist[2:20,i])))
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
      moveInd=which(dimension>maxDim & nndist[20,]<=100)
    }

    gammaNew=matrix(0, n, K)
    Pnew=rep(0,K)
    xx=xx/sum(xx)
    message('kpoints: ',kpoints,' moveInd: ',length(moveInd))

    # find kpoints nearest neighbours from gamma for each xcoord
    nnres = nabor::knn(t(gamma), query=t(xcoords), k=kpoints)
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
      for (i1 in 1:n)
        gammaNew[i1,nnidxsForThisPoint]=gammaNew[i1,nnidxsForThisPoint]+
        ((xcoords[i1,u]*xx[u])*Px)
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
