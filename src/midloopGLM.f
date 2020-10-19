CCC   Different from midloop, internal subroutine, used only for
C     family!=1;
C     Middle loop: Update the quadratic approximation likelihood
C     input:
C     n
C     m
C     x
C     y
C     mu
C     offset
C     family
C     nlambda
C     lamk
C     alpha
C     gam
C     maxit
C     del
C     thresh
C     wt              used only for family=1
C     trace
C     maxit: not used (and maxit=1) for family=1      
C     output:
C     beta
C     b0
C     yhat
C     dev
C     eta, mu

      subroutine midloopGLM(n,m,x,y,yold,weights, mu, eta, offset,
     +     family, penalty,lamk, 
     +     alpha, gam, theta, rescale, standardize, intercept, eps,
     +     maxit, thresh, nulldev, wt, beta, b0,yhat,dev,trace,convmid, 
     +     satu, ep, pll, activeset, jk)
      
      implicit none
      integer standardize, trace, penalty, maxit, i, j, jj, nmid, n, 
     +     family, m,converged,convmid, satu, intercept,
     +     rescale,activeset(m), jk, ii
      double precision x(n,m),y(n), mu(n), z(n), eta(n), wt(n), w(n), 
     +     del,olddev,weights(n),yold(n), 
     +     thresh, nulldev, dev, theta, thetaold, wtw(n),lamk(m),alpha, 
     +     gam, eps, beta(m), betaold(m), b0, b0old, yhat(n), ep, 
     +     pll(maxit), offset(n)

      do 5 jj=1, maxit
         pll(jj)=0
 5    enddo
      jj = 1
      dev = nulldev
      del = nulldev * 10
      convmid = 0
      satu = 0
 1000 if (jj .LE. maxit .AND. convmid .EQ. 0 .AND. satu .EQ. 0) then
         if(trace.EQ.1)then 
            call intpr("family=", -1, family, 1)
            call intpr("Middle loop: Update the quadratic approximation 
     +likelihood function", -1, 1, 0)
            call intpr(" middle loop iteration jj=", -1, jj, 1)
            call dblepr("convergency criterion at beginning", -1, del,1)
            call intpr("convmid", -1, convmid,1)
         endif
         call DCOPY(m, beta, 1, betaold, 1)
         
         b0old = b0
         thetaold = theta
         call glmlink(n,mu,family,theta,w, ep)
         call zeval(n, y, eta, mu, w, family, z)

         do 10 i=1, n
            wtw(i)=wt(i) * w(i)
            z(i)=z(i) - offset(i)
 10      continue
         call loop_glm(x, y, z, n, m, w, mu, penalty, thresh, eps, 
     +        standardize, intercept,family,beta, b0, lamk, alpha, gam, 
     +        weights,trace,nmid, rescale, converged, theta,
     +        pll(jj),activeset, jk)
         do 220 i = 1, n
            yhat(i) = b0
            do 230 ii = 1, jk
            j=activeset(ii)
                  yhat(i) = yhat(i) + x(i,j) * beta(j)
 230        continue
 220     continue
C        call DCOPY(n, yhat, 1, eta, 1)
         do 350 i = 1, n
            yhat(i) = yhat(i) + offset(i)
            eta(i) = yhat(i)
 350     continue
         call linkinv(n, eta, family, mu)
         olddev = dev
C     compute deviance dev
         call deveval(n, yold, mu, theta, weights, family, dev)
         if(family .EQ. 2)then
           if(dev/nulldev .LT.0.01)then
C           if(dev/nulldev .LT.0.01 .OR. nulldev .LT. 0.01)then
C          because dev in deveval depends on weights, if we scale
C          weights/n, then dev can be very small, but maybe not that
C          important      
C           if(dev/nulldev .LT.0.01 .OR. dev.LT.0.0001)then
C               call dblepr("saturated model, residual deviance = ", -1, 
C     +              dev,1)
C               call dblepr("saturated model, null deviance = ", -1, 
C     +              nulldev, 1)
C               call intpr("family", -1, family, 1)
C                call rwarn("saturated model for family='binomial', exiti
C     +ng coordinate descent iteration. Consider decreasing lambda.zero.m
C     +in.ratio if zipath is called. This may be caused by larger maxit.em
C     in zipath!") 
                satu = 1
C               call DCOPY(m, betaold, 1, beta, 1)
C               b0 = b0old
C               theta = thetaold
            endif
         endif
         del = dabs(dev - olddev)
         convmid = converged
         olddev = dev
         jj=jj + 1
         goto 1000
      endif
      if(trace.EQ.1)then
         call intpr("  Iterations used in the middle loop:", -1, jj-1,1)
         call dblepr("deviance difference at the end of middle loop "
     +        , -1, del, 1)
      endif
      
      return
      end
C     link function for glm model
C     input:
C     n is the length of mu
C     mu(n) is a vector of length n
C     family is an integer: 1 for gaussian, 2 for binomial, 
C     3 for poisson and 4 for negative binomial 
C     theta: a parameter for negative binomial family

C     output:
C     mu
C     w       
      subroutine glmlink(n, mu, family, theta, w, ep)
      
      implicit none
      integer i, n, family
      double precision mu(n), theta, ep, w(n)
C     parameter (ep=1e-5)

      do 10 i=1, n
         if(family.EQ.1)then
            w(i) = 1
         else if(family.EQ.2)then
            if(1-mu(i) .LT. ep)then
               mu(i) = 1
               w(i) = ep
            else if(mu(i) .LT. ep) then
               mu(i) = 0
               w(i) = ep
            else
               w(i)=mu(i)*(1-mu(i))
            endif
         else if(family.EQ.3)then
            w(i) = mu(i)
         else if(family.EQ.4)then
            if(mu(i) .LT. ep)then
               mu(i) = ep
C     else if(mu(i) .GT. 1000)then  ### changed 5/29/15
C     mu(i) = 1000
            endif
            w(i) = mu(i)/(1 + 1/theta * mu(i))
         endif
 10   continue
      
      return 
      end

C     link inverse function
C     input:
C     n
C     eta
C     family is an integer: 1 for gaussian, 2 for binomial, 
C     3 for poisson and 4 for negative binomial 
C     output:
C     mu      
      subroutine linkinv(n, eta, family, mu)

      implicit none
      integer n, family, i
      double precision eta(n), mu(n)
      
      do 10 i= 1, n
         if(family.EQ.1)then
            mu(i) = eta(i)
         else if(family.EQ.2)then
            mu(i) = 1.0D0/(1+exp(-eta(i)))
            if(mu(i) .LT. 1e-5)then
               mu(i) = 1e-5
            else if(mu(i) .GT. 1-1e-5)then
               mu(i) = 1-1e-5
            endif
C         else if(family.EQ.3)then
C            mu(i) = exp(eta(i))
         else if(family.EQ.3 .OR. family.EQ.4)then
            if(eta(i) .LT. log(1e-5))then
               mu(i) = 1e-5
            else
C     else if(eta(i) .LT. log(1000D0))then  ### changed 5/29/15
               mu(i) = exp(eta(i))
C     else 
C     mu(i) = 1000
            endif
         endif
 10   continue
      return
      end

C     evaluate z value
C     input:
C     n
C     y(n)
C     eta(n)
C     mu(n)
C     family is an integer: 1 for gaussian, 2 for binomial, 
C     3 for poisson and 4 for negative binomial 
C     output:
C     z(n)      
      subroutine zeval(n, y, eta, mu, w, family, z)
      
      implicit none
      integer n, i, family
      double precision y(n), eta(n), mu(n), w(n), z(n)

      do 10 i= 1, n
         if(family.EQ.1)then
            z(i) = y(i)
         else if(family.EQ.2)then
            z(i) = eta(i) + (y(i) - mu(i))/w(i)
         else if(family.EQ.3 .OR. family.EQ.4)then
            z(i) = eta(i) + (y(i)-mu(i))/mu(i)
         endif
 10   continue
      return
      end

C     deviance evaluation
C     input:
C     n
C     y(n)
C     mu(n)
C     theta
C     weights(n)
C     family is an integer: 1 for gaussian, 2 for binomial, 
C     3 for poisson and 4 for negative binomial 

C     output
C     dev

      subroutine deveval(n, y, mu, theta, weights, family, dev)
      implicit none

      integer n, family, i
      double precision dev, y(n), mu(n), theta, weights(n),tmp
      integer cisnan
      external cisnan
      
      dev = 0.0D0
      do 10 i=1, n
         if(family.EQ.1)then
            dev = dev + weights(i) * (y(i) - mu(i))**2 
         else if(family.EQ.2)then
            tmp = 0
            if(y(i) .GT. 0)then
               tmp = tmp + y(i)*dlog(y(i))
            endif
            if(mu(i) .GT. 0)then
               tmp = tmp - y(i)*dlog(mu(i))
            endif
            if(1-y(i) .GT. 0)then
               tmp = tmp + (1-y(i))*dlog(1-y(i))
            endif
            if(1-mu(i) .GT. 0)then
               tmp = tmp - (1-y(i))*dlog(1-mu(i))
            endif
            dev = dev + 2*weights(i)*tmp
         else if(family.EQ.3)then
            dev=dev+2*(weights(i)*(-y(i)+mu(i)+y(i)*dlog(max(1.0D0,y(i))
     +           /mu(i))))
         else 
            dev=dev+2*(weights(i)*(y(i)*dlog(max(1.0D0,y(i))/mu(i))-
     +           (y(i)+theta)*dlog((y(i)+theta)/(mu(i) + theta)))) 
         endif
         if(cisnan(dev).NE.0) then
             call intpr("dev is NA in Fortran src/deveval, check if some
     + columns of x have the same values", -1, 1, 1)
            call intpr("i=", -1, i, 1)
            call dblepr("y(i)=", -1, y(i), 1)
            call dblepr("mu(i)=", -1, mu(i), 1)
            call dblepr("theta=", -1, theta, 1)
            call dblepr("dev=", -1, dev, 1)
            exit
         endif
 10   continue
      return
      end

CCC   ref R/loglik.R
      subroutine loglikFor(n, y, mu, theta, w, family, ll)
      implicit none

      integer n, family, i, y0
      double precision ll, y(n), mu(n), theta, w(n)
      double precision rlgamma
      external rlgamma

      ll = 0
      do 10 i=1,n
         if(family .EQ. 4)then  !negbin
            if(y(i) .EQ. 0)then
               y0 = 1
            else 
               y0=0
            endif
            ll=ll + w(i) * (rlgamma(theta+y(i))-rlgamma(theta)-
     +           rlgamma(y(i)+1)+
     +           theta*log(theta) + y(i)*log(mu(i)+y0) - (theta + y(i))*
     +           log(theta +  mu(i)))
         else if(family .EQ. 1)then !gaussian
            ll=ll - w(i)*(y(i)-mu(i))**2 !changed 11/21/15
         else if(family .EQ. 2) then !binomial
            if(mu(i) .GT. 0 .AND. mu(i) .LT. 1)then
               ll=ll + w(i)*(y(i)*log(mu(i)/(1-mu(i)))+log(1-mu(i)))
            endif
         else if(family .EQ. 3) then !poisson
            ll=ll + w(i)*(-mu(i) + y(i)*log(mu(i))-rlgamma(y(i)+1))
         endif
 10   continue
      return
      end
      

C     linear predict value eta = a0 + x * b + offset
C input: n, m, nlambda, x, b, a0, offset, family
C output: eta, mu
      subroutine pred(n,m, nlambda,x,b,a0,offset, family, eta,mu)
      implicit none
      integer n,m,nlambda,family,k,i,j
      double precision b(m,nlambda),a0(nlambda),eta(n,nlambda)
      double precision mu(n,nlambda),x(n,m), offset(n)
      external linkinv

      do 30 k=1,nlambda
         do 20 i=1,n
            eta(i,k) = a0(k) + offset(i)
            do 10 j=1,m
               eta(i,k) = eta(i,k) + x(i,j) * b(j,k)
 10         continue
            call linkinv(1, eta(i,k), family, mu(i,k))
 20      continue
 30   continue
      return
      end

      subroutine penEval(theta, lone, ltwo, gam, penalty, res)
      implicit none
      double precision theta, lone, ltwo, gam, res
      integer penalty
      if(penalty .EQ. 2)then    !mnet
         if(abs(theta) .LE. gam * lone)then 
            res = lone * abs(theta) - theta**2/(2*gam)+0.5*ltwo*theta**2
         else 
            res = 0.5*gam*lone**2 + 0.5*ltwo*theta**2
         endif
      else if(penalty .EQ. 3)then !scad
         if(abs(theta) .LE. lone)then 
            res= lone*abs(theta) + 0.5*ltwo*theta**2
         else if(abs(theta) .LE. gam * lone)then 
            res = ((gam * lone * abs(theta) - 
     +           0.5*(theta**2+lone**2))/(gam-1) + 0.5*ltwo*theta**2)
         else 
            res = lone**2*(gam**2-1)/(2*(gam-1))+0.5*ltwo*theta**2
         endif
      else if(penalty .EQ. 1)then !enet
         res = abs(theta)*lone + 0.5*ltwo*theta**2
      endif

      return
      end

      subroutine penGLM(start, m,lambda, alpha, gam, penalty, pen)
      implicit none
      integer m, j, penalty
      double precision start(m), lambda(m), alpha, gam, res, pen

      pen = 0
      do 10 j=1, m
         call penEval(start(j), lambda(j)*alpha, 
     +        lambda(j)*(1-alpha), gam, penalty, res)
         pen = pen + res
 10   enddo
      
      return
      end
