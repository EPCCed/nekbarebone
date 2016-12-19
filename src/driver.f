c-----------------------------------------------------------------------
      program nekbarebone
      include 'TOTAL'
      
      real x(lt),f(lt),r(lt),w(lt)
      real p(lt),z(lt),c(lt),g(6,lt)
    
      logical ifbrick
      integer iel0,ielN,ielD   ! element range per proc.
      integer nx0,nxN,nxD      ! poly. order range
      integer npx,npy,npz      ! processor decomp
      integer mx,my,mz         ! element decomp
      integer icount,niter,n

      call init_comms

      call read_param(ifbrick,
     &                iel0,ielN,ielD,
     &                nx0,nxN,nxD,
     &                npx,npy,npz,
     &                mx,my,mz)

      icount = 0
            
      call comm_time(tstart)
      do nx1=nx0,nxN,nxD
        call init_dim
        do nelt=iel0,ielN,ielD
          call init_mesh(ifbrick,npx,npy,npz,mx,my,mz)
          call proxy_setupds(nx1) 
           
          call set_multiplicity(c) ! inverse of counting matrix
          call proxy_setup(ah,bh,ch,dh,zh,wh,g)
                      
          niter = 100
          n = nx1*ny1*nz1*nelt

          call set_f(f,c,n)

          call cg(x,f,g,c,r,w,p,z,n,niter)

          call comm_barrier(comms_h)

          call set_timer_flop_cnt(0)
          call cg(x,f,g,c,r,w,p,z,n,niter)
          call set_timer_flop_cnt(1)

          call gs_free(gs_h)
           
          icount = icount + 1
          mfloplist(icount) = mflops*comms_np
        enddo
      enddo
      call comm_barrier(comms_h) 
      call comm_time(tstop) 

      avmflop = 0.0
      do i = 1,icount
        avmflop = avmflop+mfloplist(i)
      enddo
      if(icount.ne.0) then
        avmflop=avmflop/icount
      endif
      if(comms_id.eq.0) then
        write(6,1) avmflop
      endif
    1    format('Avg MFlops = ',1pe12.4)

      if (comms_id.eq.0) then
        ttotal = tstop-tstart
        write(6,*) ' '
        write(6,*) ' '
        write(6,'(1(A,1p1e13.5,A,/))') 
     &       'total elapsed time: ',ttotal, ' sec'
      endif 

      call comm_finalize
      call exit(0)
      end
c--------------------------------------------------------------
      subroutine set_f(f,c,n)
      real f(n),c(n)

      do i=1,n
        arg  = 1.e9*(i*i)
        arg  = 1.e9*cos(arg)
        f(i) = sin(arg)
      enddo

      call dssum(f)
      call col2 (f,c,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_dim
      include 'SIZE'
      include 'GLOBAL'
      
C     Transfer array dimensions to common

      ny1=nx1
      nz1=nx1
 
      ndim=ldim

      return
      end
c-----------------------------------------------------------------------
      subroutine init_mesh(ifbrick,npx,npy,npz,mx,my,mz)
      include 'TOTAL'
      
      logical ifbrick
      integer e,eg,offs,npx,npy,npz,mx,my,mz
 
c     Trigger reset of mask
      cmask(-1) = 1.0

      if(.not.ifbrick) then   ! A 1-D array of elements of length P*lelt
        nelx = nelt*comms_np
        nely = 1
        nelz = 1

        npx=comms_np
        npy=1
        npz=1

        mx=nelt
        my=1
        mz=1

        if(comms_id.eq.0) then 
          write(6,*)
          write(6,*) 'Processor Distribution: npx,npy,npz=',npx,
     &               npy,npz
          write(6,*) 'Element Distribution: nelx,nely,nelz=',nelx,
     &               nely,nelz
          write(6,*) 'Local Element Distribution: mx,my,mz=',mx,
     &               my,mz
        endif
   
        do e=1,nelt
          eg = e + comms_id*nelt
          lglel(e) = eg
        enddo
      else              ! A 3-D block of elements 
        !xyz distribution of total proc if user-provided isn't valid
        if(npx*npy*npz.ne.comms_np) then
          call cubic(npx,npy,npz,comms_np) 
        endif

        !xyz distribution of total NELT if user-provided isn't valid
        if(mx*my*mz.ne.nelt)  then
          call cubic(mx,my,mz,nelt) 
        endif
      
        nelx = mx*npx
        nely = my*npy 
        nelz = mz*npz

        if(comms_id.eq.0) then 
          write(6,*)
          write(6,*) 'Processor Distribution:  npx,npy,npz=',npx,
     &               npy,npz
          write(6,*) 'Element Distribution: nelx,nely,nelz=',nelx,
     &               nely,nelz
          write(6,*) 'Local Element Distribution: mx,my,mz=',mx,
     &               my,mz
        endif

        e = 1
        offs = (mod(comms_id,npx)*mx) +
     &         npx*(my*mx)*(mod(comms_id/npx,npy)) + 
     &         (npx*npy)*(mx*my*mz)*(comms_id/(npx*npy))
        do k = 0,mz-1
        do j = 0,my-1
        do i = 0,mx-1
          eg = offs+i+(j*nelx)+(k*nelx*nely)+1
          lglel(e) = eg
          e        = e+1
        enddo
        enddo
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cubic(mx,my,mz,np)

      mx = np
      my = 1
      mz = 1
      ratio = np

      iroot3 = np**(1./3.) + 0.000001
      do i = iroot3,1,-1
        iz = i
        myx = np/iz
        nrem = np-myx*iz

        if (nrem.eq.0) then
          iroot2 = myx**(1./2.) + 0.000001
          do j=iroot2,1,-1
            iy = j
            ix = myx/iy
            nrem = myx-ix*iy
            if (nrem.eq.0) goto 20
          enddo
   20     continue

          if (ix.lt.iy) then
            it = ix
            ix = iy
            iy = it
          endif

          if (ix.lt.iz) then
            it = ix
            ix = iz
            iz = it
          endif

          if (iy.lt.iz) then
            it = iy
            iy = iz
            iz = it
          endif

          if ( REAL(ix)/iz.lt.ratio) then
            ratio = REAL(ix)/iz
            mx = ix
            my = iy
            mz = iz
          endif
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_multiplicity (c)       ! Inverse of counting matrix
      include 'GLOBAL'
      
      real c(1)

      n = nx1*ny1*nz1*nelt

      call rone(c,n)
      call gs_op(gs_h,c,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w

      do i=1,n
        c(i) = 1./c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_timer_flop_cnt(iset)
      include 'GLOBAL'
      
      real*8 time0,time1
      save time0,time1

      if (iset.eq.0) then
        flop_a  = 0
        flop_cg = 0
        call comm_time(time0)
      else
        call comm_time(time1)
        time1 = time1-time0
        if (time1.gt.0) mflops = (flop_a+flop_cg)/(1.e6*time1)
        if (comms_id.eq.0) then
          write(6,*)
          write(6,1) nelt,comms_np,nx1, nelt*comms_np
          write(6,2) mflops*comms_np, mflops
          write(6,3) flop_a,flop_cg
          write(6,4) time1
        endif
    1   format('nelt = ',i7, ', np = ', i9,', nx1 = ', i7,
     &         ', elements =', i10 )
    2   format('Tot MFlops = ', 1pe12.4, ', MFlops      = ', e12.4)
    3   format('Setup Flop = ', 1pe12.4, ', Solver Flop = ', e12.4)
    4   format('Solve Time = ', e12.4)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine read_param(ifbrick,
     &                      iel0,ielN,ielD,
     &                      nx0,nxN,nxD,
     &                      npx,npy,npz,
     &                      mx,my,mz)
      include 'SIZE'
      include 'GLOBAL'
      
      logical ifbrick
      integer iel0,ielN,ielD
      integer nx0,nxN,nxD
      integer npx,npy,npz
      integer mx,my,mz,meth

      !open .rea
      ifbrick = .false.  
      npx=0
      npy=0
      npz=0
      mx =0
      my =0
      mz =0

      if(comms_id.eq.0) then
        open(unit=9,file='data.rea',status='old') 
        read(9,*,err=100) ifbrick
        read(9,*,err=100) iel0,ielN,ielD
        read(9,*,err=100) nx0,nxN,nxD
        read(9,*,err=100,iostat=ii) npx,npy,npz !optional
        read(9,*,err=100,iostat=ii) mx,my,mz !optional
        close(9)
      endif

      call comm_barrier(comms_h)
      
      call comm_bcast(comms_h,ifbrick,4,0)

      call comm_bcast(comms_h,iel0,4,0)! ELEMENT range
      call comm_bcast(comms_h,ielN,4,0)
      call comm_bcast(comms_h,ielD,4,0)

      call comm_bcast(comms_h,nx0,4,0) ! POLY Order range
      call comm_bcast(comms_h,nxN,4,0)
      call comm_bcast(comms_h,nxD,4,0)

      call comm_bcast(comms_h,npx,4,0) ! PROC decomp
      call comm_bcast(comms_h,npy,4,0)
      call comm_bcast(comms_h,npz,4,0)

      call comm_bcast(comms_h,mx,4,0)  ! NELT decomp
      call comm_bcast(comms_h,my,4,0)
      call comm_bcast(comms_h,mz,4,0)

      if(iel0.gt.ielN .or. nx0.gt.nxN) goto 200
      if(ielN.gt.lelt .or. nxN.gt.lx1) goto 210
      if(ielD.gt.ielN .or. nxD.gt.nxN) goto 220
      
      if(comms_id.eq.0) write(6,*) "ifbrick: ",ifbrick

      return

  100 continue
      write(6,*) "ERROR READING    data.rea.....ABORT"
      write(6,*) "CHECK PARAMETERS data.rea.....ABORT"
      call exitt(1)

  200 continue
      write(6,*) "ERROR data.rea :: iel0 > ielN or nx0 > nxN :: ABORT"
      call exitt(1)

  210 continue
      write(6,*) "ERROR data.rea : ielN>lelt or nxN>lx1(SIZE) :: ABORT"
      call exitt(1)
  
  220 continue
      write(6,*) "WARNING data.rea : STRIDE   nxD>nxN or ielD>ielN !!!"
  
      nx0=4
  
      return
      end
c-----------------------------------------------------------------------
