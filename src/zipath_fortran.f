C     used in zipath.R
C     inputs: family: 3 (poisson), 4 (negbin)
C     theta
C     kx: number of variables of x having no intercept column
C     kz: number of variables of z having no intercept column
C     outputs: coefc, coefz, thetaout
      subroutine zipath_fortran(x, z, y, y1, weights, n, kx, kz, 
     +     start_count, start_zero, mustart_count, mustart_zero, 
     +     offsetx, offsetz, nlambda, lambda_count,
     +     lambda_zero, alpha_count, alpha_zero,  
     +     gam_count, gam_zero, standardize, penaltyfactor_count, 
     +     penaltyfactor_zero, maxit, eps, family,
     +     penalty, trace, coefc, coefz, yhat, iter,
     +     del, los, pll, rescale, thresh, epsbino, 
     +     theta_fixed, maxit_theta, theta, thetaout, active)
      implicit none
      integer n,m,i,ii,k,j,jj,kx, kz, penalty,nlambda,family, 
     +     standardize, maxit, y1(n), trace, iter, rfamily, 
     +     rescale, jk_count, jk_zero, active, activeset_count(kx), 
     +     activeset_zero(kz), stopit,m_count_act, maxit_theta,
     +     m_zero_act, AllocateStatus, DeAllocateStatus, 
     +     varsel_count(kx), varsel_count_old(kx), nact, conv,
     +     varsel_zero(kz), varsel_zero_old(kz), theta_fixed
      double precision x(n, kx), z(n, kz), weights(n), 
     +     start_count(kx+1), dpois, dnbinom, 
     +     start_zero(kz+1), etastart_count(n), etastart_zero(n),
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count(nlambda), thetastart, thetaout(nlambda),
     +     lambda_zero(nlambda), alpha_count, alpha_zero, gam_count, 
     +     gam_zero, eps, penaltyfactor_count(kx), y(n),
     +     penaltyfactor_zero(kz), wt(n), probi(n), thresh, epsbino, 
     +     theta, coefc(kx+1, nlambda), coefz(kz+1, nlambda), b0_x, b0z,
     +     yhat(n), d, del, los(iter,nlambda), theta0(nlambda),
     +     pll(iter, nlambda), penval, betaxall(kx), betazall(kz)
      double precision, dimension(:, :), allocatable :: x_act, z_act
      double precision, dimension(:), allocatable :: start_count_act,
     +     start_zero_act, betax, betaz,
     +     penaltyfactor_count_act, penaltyfactor_zero_act
      external :: dpois, dnbinom, gfunc

      b0_x=0
      b0z=0
      stopit = 0
C     find current active set
      call find_activeset(kx, start_count(2: (kx+1)), eps, 
     +     activeset_count, jk_count)
      call find_activeset(kz, start_zero(2: (kz+1)), eps, 
     +     activeset_zero, jk_zero)
C     When all coef are zero except intercept, choose a predictor
      if(jk_count==0)then
         jk_count = 1
         activeset_count(1)=1
      endif
      if(jk_zero==0)then
         jk_zero = 1
         activeset_zero(1)=1
      endif
      m_count_act = jk_count
      m_zero_act = jk_zero

      call gfunc(mustart_count, n, family, epsbino, etastart_count)
      call gfunc(mustart_zero, n, 2, epsbino, etastart_zero)
      do ii=1, n
         if(y1(ii) .EQ. 1)then
            probi(ii)=0
         else
            probi(ii)=mustart_zero(ii) 
            if(family .EQ. 3)then
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +              mustart_count(ii)))
            else if(family .EQ. 4)then
               probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dnbinom(0,
     +              theta, mustart_count(ii)))
            endif
         endif
      enddo

      allocate(start_count_act(jk_count+1), stat=AllocateStatus)
      allocate(start_zero_act(jk_zero+1), stat=AllocateStatus)
      allocate(penaltyfactor_count_act(jk_count),stat=AllocateStatus)
      allocate(penaltyfactor_zero_act(jk_zero), stat=AllocateStatus)
      allocate(x_act(n, jk_count), stat=AllocateStatus)
      allocate(z_act(n, jk_zero), stat=AllocateStatus)
C     call dblepr("start_count", -1, start_count, kx+1)
C     call intpr("activeset_count", -1, activeset_count, kx)
C     call dblepr("start_zero", -1, start_zero, kz+1)
C     call intpr("activeset_zero", -1, activeset_zero, kz)
      do jj=1, n
         do ii=1, jk_count
            x_act(jj, ii)=x(jj, activeset_count(ii))
         enddo
      enddo
      do jj=1, n
         do ii=1, jk_zero
            z_act(jj, ii)=z(jj, activeset_zero(ii))
         enddo
      enddo
      allocate(betax(jk_count), stat=AllocateStatus)
      allocate(betaz(jk_zero), stat=AllocateStatus)
      do 5 j=1, jk_count
         betax(j)=0
C     activeset_count(j)=j
 5    continue
      do 105 j=1, jk_zero
         betaz(j)=0
C     activeset_zero(j)=j
 105  continue
      start_count_act(1)=start_count(1)
      do ii=1, jk_count
         start_count_act(ii+1)=start_count(1+activeset_count(ii))
         penaltyfactor_count_act(ii)=
     +        penaltyfactor_count(activeset_count(ii))
      enddo
      start_zero_act(1)=start_zero(1)
      do ii=1, jk_zero
         start_zero_act(ii+1)=start_zero(1+activeset_zero(ii))
         penaltyfactor_zero_act(ii)=
     +        penaltyfactor_zero(activeset_zero(ii))
      enddo
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
            call intpr("Fortran lambda iteration i=", -1, i, 1)
         if(trace .EQ. 1)then
            call intpr("Fortran lambda iteration i=", -1, i, 1)
            call dblepr("x_act", -1, x_act(1:5, 1), 5)
            call dblepr("z_act", -1, z_act(1:5, 1), 5)
            call dblepr("y", -1, y, 5)
            call dblepr("y1", -1, y1, 5)
            call dblepr("weights", -1, weights, 5)
            call intpr("m_count_act", -1, m_count_act, 1)
            call intpr("m_zero_act", -1, m_zero_act, 1)
            call dblepr("start_count_act", -1, start_count_act,
     +           m_count_act+1)
            call dblepr("start_zero_act", -1, start_zero_act,
     +           m_zero_act+1)
            call dblepr("mustart_count=", -1, mustart_count, 5)
            call dblepr("mustart_zero=", -1, mustart_zero, 5)
            call dblepr("offsetx=", -1, offsetx, 5)
            call dblepr("offsetz=", -1, offsetz, 5)
            call dblepr("penaltyfactor_count_act=", -1,
     +           penaltyfactor_count_act, m_count_act)
            call dblepr("penaltyfactor_zero_act=", -1,
     +           penaltyfactor_zero_act, m_zero_act)
            call dblepr("betax", -1, betax, m_count_act)
            call dblepr("betaz", -1, betaz, m_zero_act)
            call dblepr("yhat", -1, yhat, 5)
            call intpr("iter=", -1, iter, 1)
            call dblepr("los=", -1, los, 5)
            call dblepr("pll=", -1, pll, 5)
            call dblepr("del=", -1, del, 1)
         endif
         j = 1
         nact = 10
         conv=0
         if(j < nact .AND. conv==0)then
             call intpr("active set, j=", -1, j, 1)
            call zi(x_act, z_act, y, y1, probi, weights, n, m_count_act,
     +           m_zero_act, start_count_act, start_zero_act, 
     +           mustart_count, mustart_zero, offsetx, offsetz, 
     +           lambda_count(i), lambda_zero(i), alpha_count, 
     +           alpha_zero, gam_count, gam_zero, standardize, 
     +           penaltyfactor_count_act, penaltyfactor_zero_act,
     +           maxit, eps, family, penalty, trace, yhat, iter, del,
     +           los, pll, rescale, thresh, epsbino, theta_fixed, 
     +           maxit_theta, theta, betax, b0_x, betaz, b0z)
            call zi(x, z, y, y1, probi, weights, n, kx,
     +           kz, start_count, start_zero, mustart_count, 
     +           mustart_zero, offsetx, offsetz, lambda_count(i), 
     +           lambda_zero(i), alpha_count, alpha_zero, gam_count,
     +           gam_zero, standardize, penaltyfactor_count, 
     +           penaltyfactor_zero, maxit, eps, family, penalty, 
     +           trace, yhat, 1, del, los, pll, rescale, thresh, 
     +           epsbino, theta_fixed, maxit_theta, theta, betaxall, 
     +           b0_x, betazall, b0z)
       call dblepr("after zi Fortran, betaxall=", -1, betaxall, kx)
       call dblepr("after zi Fortran, betazall=", -1, betazall, kz)
C     call dblepr("after zi Fortran, b0_x=", -1, b0_x, 1)
C     call dblepr("after zi Fortran, b0z=", -1, b0z, 1)
            call find_activeset(kx, betaxall, eps, activeset_count
     +           , jk_count)
            call find_activeset(kz, betazall, eps, activeset_zero
     +           , jk_zero)
C          check if converged here!
            call intpr("activeset_count", -1, activeset_count, kx)
            call intpr("activeset_zero", -1, activeset_zero, kz)
            ii=1
            if(ii <= jk_count .AND. conv==0)then
C                if(varsel_count_old(activese_count(ii)
C     +             varsel_count(ii))then
C                conv=0
C                ii = ii + 1
C                endif
            endif
C            if(conv==0)then
C                allocate x_act, start_count_act 
C             endif
            ii=1
            if(ii <= jk_zero .AND. conv==0)then
C                if((abs(betazall(ii)) > eps .AND. varind_zero(ii)==0)
C     +            .OR. (abs(betazall(ii)) <= eps .AND.
C     +             varind_zero(ii)==1))then
C                conv=0
C                ii = ii + 1
C                endif
            endif
C            if(conv==0)then
C                allocate z_act, start_zero_act 
C             endif

C            do ii=1, jk_count
C               varsel_count(ii)=varsel_count_old(activeset_count(ii))
C               varsel_count_old(ii)=varsel_count(ii)
C               if(activeset_count(ii).NE.activeset_count_old(ii))then
C                 conv=0
C               endif
C            enddo
         endif
         coefc(1, i) = b0_x
         if(jk_count .GT. 0)then
            do 200 ii = 1, m_count_act
               coefc(1+varsel_count(ii), i) = betax(ii)
 200        continue
         endif
         thetaout(i)=theta
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
     +              AllocateStatus)
               start_count_act(1) = b0_x
               do 35 ii=1, jk_count
                  start_count_act(ii+1)=betax(activeset_count(ii))
                  varsel_count(ii)=varsel_count_old(activeset_count(ii))
                  varsel_count_old(ii)=varsel_count(ii)
                  penaltyfactor_count_act(ii)=
     +                 penaltyfactor_count(varsel_count(ii))
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
            call find_activeset(m_zero_act, betaz, eps, activeset_zero,
     +           jk_zero)
            if(jk_zero .NE. m_zero_act .AND. jk_zero .GT. 0)then
               deallocate(start_zero_act, stat=DeAllocateStatus)
               allocate(start_zero_act(jk_zero+1),stat=AllocateStatus)
               deallocate(penaltyfactor_zero_act,stat=DeAllocateStatus)
               allocate(penaltyfactor_zero_act(jk_zero),stat=
     +              AllocateStatus)
               start_zero_act(1) = b0z
               do 1135 ii=1, jk_zero
                  start_zero_act(ii+1)=betax(activeset_zero(ii))
                  varsel_zero(ii)=varsel_zero_old(activeset_zero(ii))
                  varsel_zero_old(ii)=varsel_zero(ii)
                  penaltyfactor_zero_act(ii)=
     +                 penaltyfactor_zero(varsel_zero(ii))
 1135          continue
               deallocate(betaz, stat=DeAllocateStatus) 
               allocate(betax(jk_zero), stat=AllocateStatus)
               deallocate(z_act, stat=DeAllocateStatus)
               allocate(z_act(n, jk_zero), stat=AllocateStatus)
C     update x_act matrix
               do 1155 jj=1, n
                  do 1145 ii=1, jk_zero
                     z_act(jj, ii) = z(jj, varsel_zero(ii))
 1145             continue
 1155          continue
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
