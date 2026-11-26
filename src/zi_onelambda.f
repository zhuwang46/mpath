C     fit a zero inflated penalized regression for a single penalty
C     parameter. A building block used in zipath_fortran.f
C     inputs: family: 3 (poisson), 4 (negbin)
C     theta, start_count_act, start_zero_act
C     m_count_act: number of count model variables of x having no intercept column
C     m_zero_act: number of zero model variables of z having no intercept column
C     outputs: betax, b0_x, betaz, b0z, theta, start_count_act, start_zero_act
      subroutine zi_onelambda(x_act, z_act, y, y1, probi, weights, n, 
     +     m_count_act, m_zero_act, start_count_act, start_zero_act,
     +     mustart_count, mustart_zero, offsetx, offsetz,
     +     intercept,lambda_count,lambda_zero, alpha_count,alpha_zero,  
     +     gam_count, gam_zero, penaltyfactor_count_act, 
     +     penaltyfactor_zero_act, maxit, eps, family,
     +     penalty, trace, yhat, iter, del, rescale, thresh, 
     +     epsbino, theta_fixed, maxit_theta, theta, 
     +     betax, b0_x, betaz, b0z, satu)
      implicit none
      integer n,ii,k,j,penalty, family, maxit, y1(n), trace, iter, 
     +     rescale, stopit,m_count_act, maxit_theta, intercept,
     +     m_zero_act, theta_fixed, satu
      double precision weights(n), dpois, dnbinom, 
     +     etastart_count(n), etastart_zero(n),
     +     mustart_count(n), mustart_zero(n), offsetx(n), offsetz(n), 
     +     lambda_count, thetastart, lambda_zero, alpha_count, 
     +     alpha_zero, gam_count,  gam_zero, eps,
     +     penaltyfactor_count_act(m_count_act), y(n),
     +     penaltyfactor_zero_act(m_zero_act), wt(n), probi(n), thresh, 
     +     epsbino, theta, b0_x, b0z,yhat(n), d, del, theta0, penval,
     +     x_act(n, m_count_act), z_act(n, m_zero_act), pll_old, pll,
     +     start_count_act(m_count_act+1), start_zero_act(m_zero_act+1),
     +     betax(m_count_act), betaz(m_zero_act), los_old,los,
     +     zero_offset(n),penalty_count(m_count_act),
     +     penalty_zero(m_zero_act),sum_penalty
      external :: dpois, dnbinom, gfunc, ziloss

      stopit = 0
      do ii=1,n
         zero_offset(ii)=0.0D0
      enddo
C     Match the penalty-factor normalization used by glmreg_fit_fortran.
      if(m_count_act.GT.0)then
         sum_penalty=sum(penaltyfactor_count_act)
         if(sum_penalty.GT.0.0D0)then
            do ii=1,m_count_act
               penalty_count(ii)=lambda_count*
     +              penaltyfactor_count_act(ii)*m_count_act/sum_penalty
            enddo
         else
            penalty_count=0.0D0
         endif
      endif
      if(m_zero_act.GT.0)then
         sum_penalty=sum(penaltyfactor_zero_act)
         if(sum_penalty.GT.0.0D0)then
            do ii=1,m_zero_act
               penalty_zero(ii)=lambda_zero*penaltyfactor_zero_act(ii)
     +              *m_zero_act/sum_penalty
            enddo
         else
            penalty_zero=0.0D0
         endif
      endif
      k = 1
      d = 10
      call gfunc(mustart_count, n, family,epsbino,etastart_count)
      call gfunc(mustart_zero, n, 2, 0.0D0, etastart_zero)
C     The predictors obtained from mustart already include offsets.
      call ziloss(n, y, zero_offset, zero_offset, weights,
     +     etastart_count, etastart_zero, family, theta, los_old)
      penval = 0.0D0
      if(m_count_act.GT.0) call penGLM(start_count_act(2), m_count_act,
     +     penalty_count,
     +     alpha_count, gam_count, penalty, penval)
      pll_old=los_old - penval
      penval = 0.0D0
      if(m_zero_act.GT.0) call penGLM(start_zero_act(2), m_zero_act,
     +     penalty_zero,
     +     alpha_zero, gam_zero, penalty, penval)
      pll_old=pll_old - penval
 500  if(d .GT. del .AND. k .LE. iter)then
         if(trace .EQ. 1)then
            call intpr("  EM algorithm iteration k=", -1, k, 1)
            call dblepr("     d=", -1, d, 1)
         endif
         do 30 ii=1, n
            wt(ii)=weights(ii)*(1-probi(ii)) 
 30      continue
         if(family .NE. 4 .OR. theta_fixed .EQ. 1)then
            call glmreg_fit_fortran(x_act,y,wt,n,m_count_act,
     +           start_count_act,etastart_count,mustart_count,
     +           offsetx,1, lambda_count, alpha_count, gam_count,
     +           rescale,0, intercept, penaltyfactor_count_act, thresh,
     +           epsbino, maxit, eps, theta, family,  
     +           penalty, 0, betax, b0_x, yhat, satu)
         else
            thetastart = theta 
            call glmregnb_fortran(x_act,y,wt,n,m_count_act,offsetx,
     +           1, lambda_count, penalty,alpha_count, gam_count, 
     +           rescale, 0, intercept, penaltyfactor_count_act, thresh,
     +           maxit_theta, maxit, eps, epsbino, start_count_act, 
     +           etastart_count, mustart_count, thetastart, 0, 
     +           theta0, trace, betax, b0_x, theta, yhat)
         endif
C     yhat: the fitted mean values, obtained by transforming the
C     linear predictors by the inverse of the link function.
C         call dcopy(n, yhat, 1, mustart_count, 1)
C     etastart_count: linear predictors, obtained by transforming the
C     mean values by the link function.
C         call gfunc(mustart_count, n, family,epsbino,etastart_count)
         start_count_act(1) = b0_x
         if(m_count_act .GT. 0)then
            do 100 j=1, m_count_act
               start_count_act(j+1)=betax(j)
 100        continue
         endif
         do ii=1, n
            yhat(ii)=0
         enddo
         call glmreg_fit_fortran(z_act,probi,weights,n,m_zero_act,
     +        start_zero_act, etastart_zero, mustart_zero,offsetz,
     +        1, lambda_zero, alpha_zero, gam_zero, rescale,
     +        0, intercept, penaltyfactor_zero_act, thresh,
     +        epsbino, maxit, eps, theta, 2,  
     +        penalty, 0, betaz, b0z, yhat, satu)
C         call dcopy(n, yhat, 1, mustart_zero, 1)
C         call gfunc(mustart_zero, n, 2, 0.0D0, etastart_zero)
         do ii=1, n
            if(y1(ii) .EQ. 0)then
               probi(ii)=mustart_zero(ii) 
               if(family .EQ. 3)then
                  probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*dpois(0,
     +                 mustart_count(ii), 0))
               else if(family .EQ. 4)then
                  probi(ii)=probi(ii)/(probi(ii)+(1-probi(ii))*
     +                 dnbinom(0, theta, mustart_count(ii), 0))
               endif
            endif
         enddo
         
         start_zero_act(1) = b0z
         if(m_zero_act .GT. 0)then
            do 1100 j=1, m_zero_act
               start_zero_act(j+1)=betaz(j)
 1100       continue
         endif
         call ziloss(n, y, zero_offset, zero_offset, weights,
     +        etastart_count, etastart_zero, family, theta, los)
         call penGLM(betax, m_count_act, 
     +        penalty_count,
     +        alpha_count, gam_count, penalty, penval)
         pll=los - penval
         call penGLM(betaz, m_zero_act, 
     +        penalty_zero,
     +        alpha_zero, gam_zero, penalty, penval)
         pll=pll - penval
         d=abs((pll-pll_old)/pll_old)
         pll_old = pll
         k = k + 1
         goto 500
      endif

      return
      end
