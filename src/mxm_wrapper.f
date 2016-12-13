c
c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.
c
      subroutine mxm(a,n1,b,n2,c,n3)

      real a(n1,n2),b(n2,n3),c(n1,n3)

      call mxmf2(a,n1,b,n2,c,n3)

      return
      end
