c-----------------------------------------------------------------------
      subroutine init_comms
      include 'SIZE'
      include 'GLOBAL'

      integer nval
      real eps, oneeps
      
      call comm_init
      call comm_world(comms_h)
      
      call comm_id(comms_h,comms_id)
      call comm_np(comms_h,comms_np)

      ! check upper tag size limit
      call comm_tag_ub(comms_h,nval)
      if (nval.ne.0.AND.nval.lt.(10000+max(lp,lelg))) then
         if(comms_id.eq.0) write(6,*) 'ABORT: TAG_UB too small!'
         call exitt(1)
      endif
            
      IF (comms_np.GT.lp) THEN
         WRITE(6,*) 
     &   'ERROR: Code compiled for a max of',lp,' processors.'
         WRITE(6,*) 
     &   'Recompile with LP =',comms_np,' or run with fewer processors.'
         WRITE(6,*) 
     &   'Aborting in routine init_nek_comms.'
         call exitt(1)
      endif

      ! set word size for REAL
      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) then
         wdsize=8
      else
         if (comms_id.eq.0) 
     &     write(6,*) 'ABORT: single precision mode not supported!'
         call exitt(1)
      endif

      if (wdsize.eq.8) then
        call comm_type_dp(comms_real)
      else
        call comm_type_real(comms_real)
      endif
      call comm_type_int(comms_int)
      call comm_type_int8(comms_int8)
      
      ! set word size for INTEGER
      ! HARDCODED since there is no secure way to detect an int overflow
      isize = 4

      ! set word size for LOGICAL
      lsize = 4

      ! set word size for CHARACTER
      csize = 1

      if (comms_id.eq.0) then 
         write(6,*) 'Number of processors:',comms_np
         WRITE(6,*) 'REAL    wdsize      :',wdsize
         WRITE(6,*) 'INTEGER wdsize      :',isize
      endif

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Global vector commutative operations for reals, ints and int8s
c-----------------------------------------------------------------------
      subroutine gop(x, w, op, n)
      include 'GLOBAL'

      real x(n), w(n)
      character*3 op

      if (op.eq.'+  ') then
         call comm_allreduce_add(comms_h,comms_real,x,n,w)
      elseif (op.EQ.'M  ') then
         call comm_allreduce_max(comms_h,comms_real,x,n,w)
      elseif (op.EQ.'m  ') then
         call comm_allreduce_min(comms_h,comms_real,x,n,w)
      elseif (op.EQ.'*  ') then
         call comm_allreduce_mul(comms_h,comms_real,x,n,w)
      else
         write(6,*) comms_id,' OP ',op,' not supported.  ABORT in GOP.'
         call exitt(1)
      endif

      call copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine igop(x, w, op, n)
      include 'GLOBAL'

      integer x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call comm_allreduce_add(comms_h,comms_int,x,n,w)
      elseif (op.EQ.'M  ') then
        call comm_allreduce_max(comms_h,comms_int,x,n,w)
      elseif (op.EQ.'m  ') then
        call comm_allreduce_min(comms_h,comms_int,x,n,w)
      elseif (op.EQ.'*  ') then
        call comm_allreduce_mul(comms_h,comms_int,x,n,w)
      else
        write(6,*) comms_id,' OP ',op,' not supported.  ABORT in igop.'
        call exitt(1)
      endif

      call icopy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine i8gop(x, w, op, n)
      include 'GLOBAL'

      integer*8 x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call comm_allreduce_add(comms_h,comms_int8,x,n,w)
      elseif (op.EQ.'M  ') then
        call comm_allreduce_max(comms_h,comms_int8,x,n,w)
      elseif (op.EQ.'m  ') then
        call comm_allreduce_min(comms_h,comms_int8,x,n,w)
      elseif (op.EQ.'*  ') then
        call comm_allreduce_mul(comms_h,comms_int8,x,n,w)
      else
        write(6,*) comms_id,' OP ',op,' not supported.  ABORT in igop.'
        call exitt(1)
      endif

      call i8copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt(err)
      include 'GLOBAL'

      integer err

      call comm_barrier(comms_h)
      call comm_finalize
      call exit(err)

      return
      end
c-----------------------------------------------------------------------
