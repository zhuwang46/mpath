C similar to and derived from R/glmreg_fit
C input: etastart, mustart - changed
      subroutine glmreg_fit_fortran(x, y, weights, n,m,start, etastart,
     +         mustart, offset, nlambda, lambda, alpha, 
     +         gam, rescale, standardize, penaltyfactor,
     +         thresh, epsbino, maxit, 
     +         eps, theta, family, penalty, trace, beta, b0, yhat)
      implicit none
      integer n,m, i,j, penalty, nlambda, family, standardize, maxit,
     +     innermaxit, trace, rescale, good,
     +     satu, convout(nlambda), startv
      double precision x(n, m), y(n), weights(n),start(m+1),etastart(n),
     +     mustart(n), offset(n), lambda(nlambda), lam(m, nlambda), 
     +     alpha, gam, thresh, epsbino, eps, theta, penaltyfactor(m),
     +     beta(m, nlambda), b0(nlambda), sumpen, outpll(maxit,nlambda),
     +     wt(n),eta(n), dev, meanx(m), normx(m),xd(m), nulldev,
     +     penfac(m), resdev(nlambda), yhat(n), mu(n), sumwt,
     +     crossprod_beta(nlambda), meany, meanoffset
    
      if(family.EQ.1) then
            rescale = 0
      endif
C### this theta is not useful but required as input in Fortran glmlink subroutine
      if(family .NE. 4)then
              theta = 1
      endif
C      call  deveval(n, y, mustart, theta, weights, family, dev)
      startv = 1
      do 90 j=1, nlambda
      b0(j) = start(1)
        do 100 i=1, m
         beta(i, j) = start(i+1)
100     continue 
90    continue 
C    resdev <- rep(0, nlambda)
C    yhat <- matrix(0, nobs, nlambda)
C    penfac <- penalty.factor/sum(penalty.factor) * m
C    lam <- penfac %o% lambda
      sumpen = sum(penaltyfactor)
      do 110 i=1, m
      penfac(i) = penaltyfactor(i)/sumpen * m
110   continue 
C generate outer product of two vectors penfac and lambda, return lam
      call outprod(m, penfac, nlambda, lambda, lam)
      if(family.EQ.1)then 
        innermaxit = maxit
        maxit = 1
      else 
              innermaxit = 1
      endif
                    
      sumwt = sum(weights)
      do 130 i=1, n
      wt(i) = weights(i)/sumwt
130   continue
C    wtnew <- weights/sum(weights)
C    meanx <- drop(wtnew %*% x)
C     compute weighted column averages meanx = x^(transpose) * wt
      call DGEMV('T',n, m, 1.0D0, x, n, wt, 1, 0.0D0, meanx, 1)

C    if(standardize){
C        xx <- scale(x, meanx, FALSE) # centers x
C        xx <- sqrt(wtnew) * xx
C        one <- rep(1,n)  
C        normx <- sqrt(one %*% xx^2)
C        xx <- scale(x, meanx, normx)
C    }
C choose family=2 (any interger except 1 to avoid centering/scaling y) below
      if(standardize .EQ. 1)then
              call preprocess(x,y,n,m,weights, 2,standardize,normx,xd,
     +     mu)
      endif
      call outloop(x,y,weights,wt,n,m,penalty,nlambda,lam,alpha,gam,
     +            theta,rescale,mustart,eta,offset,family,standardize,
     +            nulldev,thresh,maxit,innermaxit,eps,trace,start,
     +            startv,beta,b0,resdev,yhat,
     +              convout, satu, good, epsbino,outpll) 
      if (standardize .EQ. 1)then
         do 200 j=1, nlambda
         do 250 i=1, m
                beta(i, j)=beta(i, j)/normx(i)
 250      continue
 200      continue
C     compute crossproduct, like crossprod(meanx, beta) in R
C http://www.tat.physik.uni-tuebingen.de/~kley/lehre/ftn77/tutorial/blas.html
         call DGEMV('T',m, nlambda, 1.0D0, beta, m, meanx, 1, 0.0D0, 
     +             crossprod_beta, 1)
         do 300 j=1, nlambda
                b0(j)=b0(j) - crossprod_beta(j)
 300      continue
         if (family .EQ. 1)then
                meany=sum(y)/n
            do 400 j=1, nlambda
                b0(j)=b0(j) + meany - meanoffset
 400        continue
         endif
      endif
C update mustart
      call linkinv(n, yhat, family, mustart)

      return
      end