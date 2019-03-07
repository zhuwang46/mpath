C used in ncvreg.R
      subroutine ncvreg_fortran(x, y, weights, n,m,start, etastart,
     +         mustart, offset, nlambda, lambda, alpha, 
     +         gam, standardize, penaltyfactor,
     +         maxit, eps, penalty, trace, beta, b0, yhat, iter,
     +         del, rfamily, s)
      implicit none
      integer n,m,i,ii,k,j, penalty,nlambda,family, standardize, maxit,
     +  trace, iter, rfamily
      double precision x(n, m), y(n), weights(n),start(m+1),etastart(n),
     +     mustart(n), offset(n), lambda(nlambda), lam(m, nlambda), 
     +     alpha, gam, eps, penaltyfactor(m),
     +     beta(m, nlambda), b0(nlambda), beta_1(m), b0_1,
     +     penfac(m), yhat(n), mu(n), d, del,
     +     fk_old(n), s, B, h(n), fk(n), a, los(iter,nlambda), 
     +     pll(iter, nlambda)
    
      i=1
      b0_1=0
      do 5 j=1, m
      beta_1(j)=0
5     continue
10    if(i .LE. nlambda)then
       k = 1
500        if(d .GT. del .AND. k .LE. iter)then
           call dcopy(n, yhat, 1, fk_old, 1)
           call compute_h(rfamily, y, fk_old, s, B, h)
C          check if h has NAN value
           do 30 ii=1, n
            if(h(ii) .NE. h(ii))then
                    exit
            endif
30         continue
           call glmreg_fit_fortran(x, h, weights, n, m, start, 
     +          etastart, mustart, offset, 1, lambda(i)/B, alpha,
     +          gam, 0, 0, penaltyfactor, 0, 0, maxit, eps, 0, 1, 
     +          penalty, trace, beta_1, b0_1, yhat)
C           call dcopy(n, yhat, 1, fk, 1)
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
