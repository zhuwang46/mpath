C used in nclreg.R
      subroutine nclreg_fortran(x, y, weights, n,m,start, etastart,
     +         mustart, offset, nlambda, lambda, alpha, 
     +         gam, standardize, penaltyfactor, maxit, eps, family,
     +         penalty, trace, beta, b0, yhat, iter,
     +         del, rfamily, B, s, rescale, thresh, epsbino, theta)
      implicit none
      integer n,m,i,ii,k,j, penalty,nlambda,family, standardize, maxit,
     +  trace, iter, rfamily, rescale
      double precision x(n, m), y(n), weights(n),start(m+1),etastart(n),
     +     mustart(n), offset(n), lambda(nlambda),
     +     alpha, gam, eps, penaltyfactor(m), thresh, epsbino,  theta,
     +     beta(m, nlambda), b0(nlambda), beta_1(m), b0_1,
     +     yhat(n), d, del, lambda_i, 
     +     fk_old(n), s, B, h(n), fk(n), a, los(iter,nlambda), 
     +     pll(iter, nlambda)
    
      call dblepr("     del=", -1, del, 1)
      call intpr("     iter=", -1, iter, 1)
      i=1
      b0_1=0
      do 5 j=1, m
      beta_1(j)=0
5     continue
10    if(i .LE. nlambda)then
             call intpr("i=", -1, i, 1)
       k = 1
       d = 10
500        if(d .GT. del .AND. k .LE. iter)then
C             if(trace .EQ. 1)then
               call intpr("  k=", -1, k, 1)
               call dblepr("     d=", -1, d, 1)
C             endif
           call dcopy(n, yhat, 1, fk_old, 1)
           call compute_h(rfamily, n, y, fk_old, s, B, h)
           call dblepr("h=", -1, h, n)
C          check if h has NAN value
           do 30 ii=1, n
            if(h(ii) .NE. h(ii))then
                    exit
            endif
30         continue
            lambda_i = lambda(i)/B
           call glmreg_fit_fortran(x, h, weights, n, m, start, 
     +          etastart, mustart, offset, 1, lambda_i, alpha,
     +          gam, rescale, standardize, penaltyfactor, thresh,
     +          epsbino, maxit, eps, theta, family,  
     +          penalty, trace, beta_1, b0_1, yhat)
           call dcopy(n, yhat, 1, fk, 1)
           call dblepr("beta_1=", -1, beta_1, m)
           call dblepr("yhat=", -1, yhat, 5)
           call dblepr("etastart=", -1, etastart, 5)
           call dblepr("mustart=", -1, mustart, 5)
           call dblepr("start=", -1, start, m+1)
C           call dcopy(n, yhat, 1, etastart, 1)
           start(1) = b0_1
           do 100 j=1, m
           start(j+1)=beta_1(j)
100         continue
C may add call penGLM and loss values in the future
           a = 0
           do 120 ii=1, n
            a=a+(fk_old(i) - fk(i))**2
120         continue
            d = a
            k = k + 1
           goto 500
           endif
           b0(i) = b0_1
           do 200 j=1, m
            beta(j, i) = beta_1(j)
200        continue 
      i = i + 1
      goto 10
      endif

      return
      end
