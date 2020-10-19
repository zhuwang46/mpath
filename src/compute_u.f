C used in ccglmreg_onelambda, step 3 in the COCO Algorithm
C input dfun, if dfun=1-3, for regression; otherwise, if dfun=4-7 for classification
C output u
      subroutine compute_u(dfun, n, y, fk_old, u)
              implicit none
              integer dfun, i, n
              double precision y(n), fk_old(n), u(n)

      do 10 i=1, n
      if(dfun .LE. 3)then
              u(i)=y(i)-fk_old(i)
        else 
              u(i)=y(i)*fk_old(i)
      end if
 10           continue
             
      return
      end 
