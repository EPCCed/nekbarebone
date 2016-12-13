c-----------------------------------------------------------------------
      subroutine dssum(f)
      include 'GLOBAL'
      
      real f(1)

      call gs_op(gs_h,f,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w

      return
      end
c-----------------------------------------------------------------------
      subroutine proxy_setupds(nx)
      include 'TOTAL'

      integer dof
      integer*8 glo_num(lx1*ly1*lz1*lelt)
      real*8 t0, t1

      call comm_time(t0)
      call set_vert_box(glo_num,nx) ! Set global-to-local map
    
      ntot      = nx*nx*nx*nelt ! assumes nx=ny=nz
      call gs_setup(gs_h,glo_num,ntot,comms_h,mp) ! Initialize gather-scatter
      dof = ntot*mp
      call comm_time(t1)
      t1 = t1 - t0
      if (comms_id.eq.0) then
         write(6,1) t1,gs_h,nx,dof
    1    format('   setupds time',1pe11.4,' seconds ',2i3,i12)
      endif

      return
      end
c-----------------------------------------------------------------------
c     Set up global numbering for elements in a box
      subroutine set_vert_box(glo_num,nx)
      include 'TOTAL'

      integer*8 glo_num(1),ii,kg,jg,ig ! The latter 3 for proper promotion
      integer e,ex,ey,ez,eg

      nn = nx-1  ! nn := polynomial order

      do e=1,nelt
        eg = lglel(e)                              
        call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)  
        do k=0,nn
        do j=0,nn
        do i=0,nn
           kg = nn*(ez-1) + k                     
           jg = nn*(ey-1) + j                     
           ig = nn*(ex-1) + i
           ii = 1 + ig + jg*(nn*nelx+1) + kg*(nn*nelx+1)*(nn*nely+1) 
           ll = 1 + i + nx*j + nx*nx*k + nx*nx*nx*(e-1)
           glo_num(ll) = ii
        enddo
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
      integer ex,ey,ez,eg

      nelxy = nelx*nely
 
      ez = 1 +  (eg-1)/nelxy
      ey = mod1 (eg,nelxy)
      ey = 1 +  (ey-1)/nelx
      ex = mod1 (eg,nelx)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_glo_num(glo_num,nx)
      include 'TOTAL'

      integer*8 glo_num(lx1*ly1*lz1,lelt)

      integer e

      do e=1,nelt
         call outmat_e_i8(glo_num(1,e),e,nx)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_e_i8(gn,e,nx)
      include 'TOTAL'

      integer*8 gn(lx1,ly1,lz1)

      integer e

      write(6,*)
      write(6,2) e
      write(6,*)

      do k0=3,1,-2

         k1=k0+1
            write(6,*) k0,k1
         do j=nx,1,-1
            write(6,1) ((gn(i,j,k),i=1,nx),k=k0,k1)
         enddo
         write(6,*)

      enddo
    1 format('gn:',4i5,3x,4i5)
    2 format('gn: element: ',i4)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_glo_num_general(glo_num,nx)
      include 'TOTAL'

      integer*8 glo_num(1)

      integer e

      io = 0
      do e=1,nelt
         write(6,*)
         write(6,2) e
         write(6,*)

         do k=nx,1,-1
              write(6,*) 'k = ',k
            do j=nx,1,-1
                write(6,1) (glo_num(igo),igo=1+io,io+nx)
                io = io+nx
            enddo
         enddo

      enddo
    1 format('gn:',6i5)
    2 format('gn: element: ',i4)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_r(x,name5)
      include 'TOTAL'
      
      character*5 name5

      real x(lx1*ly1*lz1,lelt)

      integer e

      do e=1,nelt
         call outmat_e_r(x(1,e),name5,e)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_e_r(x,name5,e)
      include 'TOTAL'
      
      character*5 name5

      real x(lx1,ly1,lz1)

      integer e

      write(6,*)
      write(6,2) e,name5
      write(6,*)

      do k0=3,1,-2

         k1=k0+1
         do j=ny1,1,-1
            write(6,1) ((x(i,j,k),i=1,4),k=k0,k1)
         enddo
         write(6,*)

      enddo
    1 format('mat: ',4f8.3,3x,4f8.3)
    2 format('mat: element: ',i4,2x,a5)
 
      return
      end
c-----------------------------------------------------------------------
