C     loss subroutine converted from R functions ziPoisson and ziNegBin
C     in package mpath/zipath.R
C     family=3 - ziPoisson; 4 - ziNegBin
C     output: los is the average loss

      subroutine ziloss(n, y, weights, fc, fz, family, theta, los,
     +       loglik0, loglik1)
      implicit none
      integer family, i, n
      double precision x, y(n), weights(n), fc(n), fz(n), mu, 
     +     phi, theta, los, dpois, dnbinom, loglik0(n), loglik1(n)
      external dpois, dnbinom

      los = 0.d0
      do 20 i=1, n
         if(family .EQ. 3)then
            mu=dexp(fc(i))
            call linkinv(1, fz(i), 2, phi)
            if(i==915)then
                call intpr("i=", -1, i, 1)
                call dblepr("y(i)=", -1, y(i), 1)
                call dblepr("mu=", -1, mu, 1)
                call dblepr("phi=", -1, phi, 1)
                call dblepr("dlog(phi+dexp(dlog(1.0-phi)-mu))",
     +          -1, dlog(phi+dexp(dlog(1.0-phi)-mu)), 1)
                call dblepr("dpois(y(i), mu, 1)",
     +               -1,dpois(int(y(i)),mu,1),1) 
                call dblepr("dlog(1-phi)+dpois(y(i), mu, 1))", -1, 
     +           dlog(1-phi)+dpois(int(y(i)), mu, 1), 1) 
            endif
               loglik0(i)=dlog(phi+dexp(dlog(1.0-phi)-mu))
               loglik1(i)=dlog(1-phi)+dpois(int(y(i)), mu, 1)
            if(y(i) < 1)then
               los=los+weights(i)*dlog(phi+dexp(dlog(1.0-phi)-mu))
            else
               los=los+weights(i)*(dlog(1-phi)+dpois(int(y(i)), mu, 1))
            endif
         else if(family .EQ. 4)then
            mu=dexp(fc(i))
            call linkinv(1, fz(i), 2, phi)
            if(y(i) < 1)then
               los=los+weights(i)*dlog(phi+dexp(dlog(1.0-phi)+
     +                 dnbinom(0, theta, mu, 1)))
            else
               los=los+weights(i)*(dlog(1-phi)+
     +          dnbinom(int(y(i)),theta,mu,1))
            endif
         endif
 20   continue

      return
      end
