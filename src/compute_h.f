C converted from R mpath/R/compute.h function
C input rfamily, y, fk_old, s, B -- no changes
C output h
      subroutine compute_h(rfamily, n, y, fk_old, s, B, h)
              implicit none
              integer rfamily, i, n
              double precision y(n), fk_old(n), s, B, u(n), g(n), h(n)

      if(rfamily .EQ. 1)then
              do 10 i=1, n
              u(i)=y(i)-fk_old(i)
              call gradient(rfamily, n, u(i), s, g(i))
              h(i) = g(i)/B + fk_old(i)
 10   continue
      else
         if(rfamily .EQ. 2 .OR. rfamily .EQ. 3 .OR. rfamily .EQ. 4)then
              do 30 i=1, n
              u(i)=y(i)*fk_old(i)
              call gradient(rfamily, n, u(i), s, g(i))
              h(i) = -y(i)*g(i)/B + fk_old(i)
  30  continue
         endif
      endif
             
            return
            end 
