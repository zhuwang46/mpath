C     used in zipath.R
C     inputs: family: 3 (poisson), 4 (negbin)
      subroutine zipath_fortran(x, z, y, y1, weights, n, kx, kz, 
     +     start_count, start_zero, mustart_count, mustart_zero, 
     +     offsetx, offsetz, nlambda, lambda_count,
     +     lambda_zero, alpha_count, alpha_zero,  
     +     gam_count, gam_zero, standardize, penaltyfactor_count, 
     +     penaltyfactor_zero, maxit, eps, family,
     +     penalty, trace, coefc, coefz, yhat, iter,
     +     del, los, pll, rescale, thresh, epsbino, 
     +     theta, active)
      implicit none
      integer n,m,i,ii,k,j,jj,kx, kz, penalty,nlambda,family, 
     +     standardize, maxit, y1(n), trace, iter, rfamily, 
     +     rescale, jk_count, jk_zero, active, activeset_count(kx), 
     +     activeset_zero(kz), stopit,m_count_act,
     +     m_zero_act, AllocateStatus, DeAllocateStatus, 
     +     varsel_count(kx), varsel_count_old(kx),
     +     varsel_zero(kz), varsel_zero_old(kz)
      double precision x(n, kx), z(n, kz), weights(n), 
     +     start_count(kx+1), dpois,  
     +     start_zero(kz+1), etastart_count(n), etastart_zero(n),
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count(nlambda),
     +     lambda_zero(nlambda), alpha_count, alpha_zero, gam_count, 
     +     gam_zero, eps, penaltyfactor_count(kx), y(n),
     +     penaltyfactor_zero(kz), wt(n), probi(n), thresh, epsbino, 
     +     theta, coefc(kx+1, nlambda), coefz(kz+1, nlambda), b0x, b0z,
     +     yhat(n), d, del, fk(n),
     +     fk_old(n), a, los(iter,nlambda), 
     +     pll(iter, nlambda), penval
      double precision, dimension(:, :), allocatable :: x_act, z_act
      double precision, dimension(:), allocatable :: start_count_act,
     +     start_zero_act, betax, betaz,
     +     penaltyfactor_count_act, penaltyfactor_zero_act
      external :: dpois, gfunc

      b0x=0
      b0z=0
      stopit = 0
      m_count_act = kx
      m_zero_act = kz
      jk_count = kx
      jk_zero = kz

      if(family .EQ. 3)then
         do ii=1, n
            if(y1(ii) .EQ. 1)then
               probi(ii)=0
            else
              probi(ii)=mustart_zero(ii) 
              probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +                  mustart_count(ii)))
            endif
         enddo
      endif
      allocate(start_count_act(kx+1), stat=AllocateStatus)
      allocate(start_zero_act(kz+1), stat=AllocateStatus)
      allocate(penaltyfactor_count_act(kx), stat=AllocateStatus)
      allocate(penaltyfactor_zero_act(kz), stat=AllocateStatus)
      allocate(x_act(n, kx), stat=AllocateStatus)
      allocate(z_act(n, kz), stat=AllocateStatus)
      call copymatrix(n, kx, x, x_act)
      call copymatrix(n, kz, z, z_act)

      allocate(betax(m_count_act), stat=AllocateStatus)
      allocate(betaz(m_zero_act), stat=AllocateStatus)
      do 5 j=1, kx
         betax(j)=0
         activeset_count(j)=j
 5    continue
      do 105 j=1, kz
         betaz(j)=0
         activeset_zero(j)=j
 105    continue
      call dcopy(kx+1, start_count, 1, start_count_act, 1)
      call dcopy(kz+1, start_zero, 1, start_zero_act, 1)
      call dcopy(kx, penaltyfactor_count, 1, penaltyfactor_count_act, 1)
      call dcopy(kz, penaltyfactor_zero, 1, penaltyfactor_zero_act, 1)
      do 101 j=1, kx
         varsel_count_old(j)=j
         varsel_count(j)=j
 101  continue
      do 102 j=1, kz
         varsel_zero_old(j)=j
         varsel_zero(j)=j
 102  continue
      i=1
 10   if(i .LE. nlambda)then
C         if(trace .EQ. 1)then
            call intpr("i=", -1, i, 1)
C         endif
         k = 1
         d = 10
         call dblepr("del=", -1, del, 1)
 500     if(d .GT. del .AND. k .LE. iter)then
            a = 0
C      if(trace .EQ. 1)then
            call intpr("  k=", -1, k, 1)
            call intpr("  activeset", -1, activeset_count, jk_count)
            call dblepr("     d=", -1, d, 1)
C      endif
C            call dcopy(n, yhat, 1, fk_old, 1)
            do 30 ii=1, n
               wt(ii)=weights(ii)*(1-probi(ii))/n 
 30         continue
            call dblepr("wt", -1, wt, 10)
           call intpr("before glmreg_fit_fortran",-1,1, 1)
           call dblepr("start_count_act",-1,start_count_act,m_count_act)
           call dblepr("etastart_count",-1,etastart_count, 10)
           call dblepr("mustart_count",-1,mustart_count, 10)
            call glmreg_fit_fortran(x_act,y,wt,n,m_count_act,
     +           start_count_act,etastart_count,mustart_count,offsetx,
     +           1, lambda_count(i), alpha_count, gam_count, rescale, 
     +           0, penaltyfactor_count_act, thresh,
     +           epsbino, maxit, eps, theta, family,  
     +           penalty, trace, betax, b0x, yhat)
            call dcopy(n, yhat, 1, mustart_count, 1)
            call gfunc(yhat, n, family, epsbino, etastart_count)
           call intpr("after glmreg_fit_fortran",-1,1, 1)
           call dblepr("start_count_act",-1,start_count_act,m_count_act)
           call dblepr("etastart_count",-1,etastart_count, 10)
           call dblepr("mustart_count",-1,mustart_count, 10)
           call dblepr("yhat",-1, yhat, 10)
           a=a+(start_count_act(1) - b0x)**2
           call dblepr("betax",-1,betax, m_count_act)
           call dblepr("b0x",-1,b0x, 1)
            call dcopy(n, yhat, 1, fk, 1)
            start_count_act(1) = b0x
            if(jk_count .GT. 0)then
               do 100 j=1, m_count_act
               a=a+(start_count_act(j+1) - betax(j))**2
                  start_count_act(j+1)=betax(j)
 100           continue
            endif
C            call dblepr("after gfunc,etastart_zero",-1,etastart_zero, n)
C            call dblepr("start_zero_act", -1, start_zero_act,m_zero_act)
C            call dblepr("lambda_zero(i)", -1, lambda_zero(i), 1)
C            call dblepr("alpha_zero(i)", -1, alpha_zero, 1)
C            call dblepr("penaltyfactor_zero_act", -1,
C     +                  penaltyfactor_zero_act, m_zero_act)
C            call dblepr("z_act", -1, z_act(1, 1:m_zero_act), m_zero_act)
            call dblepr("mustart_zero", -1, mustart_zero, 10)
            call glmreg_fit_fortran(z_act,probi,weights,n,m_zero_act,
     +           start_zero_act, etastart_zero, mustart_zero,offsetz,
     +           1, lambda_zero(i), alpha_zero, gam_zero, rescale,
     +           standardize, penaltyfactor_zero_act, thresh,
     +           epsbino, maxit, eps, theta, 2,  
     +           penalty, trace, betaz, b0z, yhat)
            call dcopy(n, yhat, 1, mustart_zero, 1)
            call gfunc(yhat, n, 2, epsbino, etastart_zero)
           call dblepr("betaz",-1,betaz, m_zero_act)
           call dblepr("b0z",-1,b0z, 1)
           call dblepr("yhat",-1, yhat, 10)
           call dblepr("mustart_zero",-1, mustart_zero, 10)
           if(family .EQ. 3)then
             do ii=1, n
              if(y1(ii) .EQ. 0)then
               probi(ii)=yhat(ii) 
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +                  mustart_count(ii)))
              endif
             enddo
           endif
            call dblepr("probi updated", -1, probi, 10)
            
           a=a+(start_zero_act(1) - b0z)**2
           start_zero_act(1) = b0z
            if(jk_zero .GT. 0)then
               do 1100 j=1, m_zero_act
               a=a+(start_zero_act(j+1) - betaz(j))**2
                  start_zero_act(j+1)=betaz(j)
1100           continue
            endif
            if(trace .EQ. 1)then
C     call loss(n, y, fk, cost, rfamily, s, los(k, i))
               penval = 0.d0
               call penGLM(betax, m_count_act, 
     +              lambda_count(i)*penaltyfactor_count_act, 
     +              alpha_count, gam_count, penalty, penval)
               if(standardize .EQ. 1)then
                  pll(k, i)=los(k, i) + n*penval
               else 
                  pll(k, i)=los(k, i) + penval
               endif
            endif
C            do 120 ii=1, n
C               a=a+(fk_old(i) - fk(i))**2
C 120        continue
            d = a
            k = k + 1
            goto 500
         endif
C     needs to store the b0 result here
C     b0(i) = b0x
         coefc(1, i) = b0x
         if(jk_count .GT. 0)then
            do 200 ii = 1, m_count_act
               coefc(1+varsel_count(ii), i) = betax(ii)
 200        continue
         endif
         coefz(1, i) = b0z
         if(jk_zero .GT. 0)then
            do 210 ii = 1, m_zero_act
               coefz(1+varsel_zero(ii), i) = betaz(ii)
 210        continue
         endif
         if(active .EQ. 1)then
            call find_activeset(m_count_act, betax, eps, activeset_count
     +           , jk_count)
C     this activeset is relative to the current x_act, but how about
C     relative to x instead? compute varsel_count for true index in x
            if(jk_count .NE. m_count_act .AND. jk_count .GT. 0)then
               deallocate(start_count_act, stat=DeAllocateStatus)
               allocate(start_count_act(jk_count+1),stat=AllocateStatus)
               deallocate(penaltyfactor_count_act,stat=DeAllocateStatus)
               allocate(penaltyfactor_count_act(jk_count),stat=
     +                 AllocateStatus)
               start_count_act(1) = b0x
               do 35 ii=1, jk_count
                  start_count_act(ii+1)=betax(activeset_count(ii))
                  varsel_count(ii)=varsel_count_old(activeset_count(ii))
                  varsel_count_old(ii)=varsel_count(ii)
                  penaltyfactor_count_act(ii)=
     +                    penaltyfactor_count(varsel_count(ii))
 35            continue
               deallocate(betax, stat=DeAllocateStatus) 
               allocate(betax(jk_count), stat=AllocateStatus)
               deallocate(x_act, stat=DeAllocateStatus)
               allocate(x_act(n, jk_count), stat=AllocateStatus)
C     update x_act matrix
               do 55 jj=1, n
                  do 45 ii=1, jk_count
                     x_act(jj, ii) = x(jj, varsel_count(ii))
 45               continue
 55            continue
               m_count_act = jk_count
            endif
            call find_activeset(m_zero_act, betaz, eps, activeset_zero
     +           , jk_zero)
            if(jk_zero .NE. m_zero_act .AND. jk_zero .GT. 0)then
               deallocate(start_zero_act, stat=DeAllocateStatus)
               allocate(start_zero_act(jk_zero+1),stat=AllocateStatus)
               deallocate(penaltyfactor_zero_act,stat=DeAllocateStatus)
               allocate(penaltyfactor_zero_act(jk_zero),stat=
     +                 AllocateStatus)
               start_zero_act(1) = b0z
               do 1135 ii=1, jk_zero
                  start_zero_act(ii+1)=betax(activeset_zero(ii))
                  varsel_zero(ii)=varsel_zero_old(activeset_zero(ii))
                  varsel_zero_old(ii)=varsel_zero(ii)
                  penaltyfactor_zero_act(ii)=
     +                    penaltyfactor_zero(varsel_zero(ii))
1135           continue
               deallocate(betaz, stat=DeAllocateStatus) 
               allocate(betax(jk_zero), stat=AllocateStatus)
               deallocate(z_act, stat=DeAllocateStatus)
               allocate(z_act(n, jk_zero), stat=AllocateStatus)
C     update x_act matrix
               do 1155 jj=1, n
                  do 1145 ii=1, jk_zero
                     z_act(jj, ii) = z(jj, varsel_zero(ii))
1145               continue
1155            continue
               m_zero_act = jk_zero
            endif
         endif
         i = i + 1
         goto 10
      endif
      deallocate(betax, start_count_act, x_act, 
     +     penaltyfactor_count_act)
      deallocate(betaz, start_zero_act, z_act, 
     +     penaltyfactor_zero_act)

      return
      end
