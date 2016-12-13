#include <stddef.h> /* for size_t */
#include <stdlib.h> /* for exit */
#include <stdio.h> /* for exit */
#include <string.h> /* memcpy */
#include <limits.h> /* for gs identities */
#include <float.h>  /* for gs identities */
#include "name.h"
#include "fail.h"
#include "types.h"
#include "tensor.h"
#include "gs_defs.h"
#include "gs_local.h"
#include "comm.h"


void comm_init()
{
#ifdef MPI
  MPI_Init(0,0);
#endif
}

void comm_finalize()
{
#ifdef MPI
  MPI_Finalize();
#endif
}


void comm_dup(const comm_hdl ch, comm_hdl *chp)
{
#ifdef MPI
  comm_hdl chd;
  MPI_Comm_dup(ch, &chd);
  *chp = chd;
#endif
}

void comm_world(comm_hdl *chp)
{
#ifdef MPI
  *chp = MPI_COMM_WORLD;
#endif
}

void comm_free(comm_hdl ch)
{
#ifdef MPI
  MPI_Comm_free(&ch);
#endif
}


void comm_id(const comm_hdl ch, int *id)
{
#ifdef MPI
  MPI_Comm_rank(ch,id);
#endif
}

void comm_np(const comm_hdl ch, int *np)
{
#ifdef MPI
  MPI_Comm_size(ch,np);
#endif
}


void comm_type_int(comm_type *ct)
{
#ifdef MPI
  *ct = MPI_INTEGER;
#endif
}

void comm_type_int8(comm_type *ct)
{
#ifdef MPI
  *ct = MPI_INTEGER8;
#endif
}

void comm_type_real(comm_type *ct)
{
#ifdef MPI
  *ct = MPI_REAL;
#endif
}

void comm_type_dp(comm_type *ct)
{
#ifdef MPI
  *ct = MPI_DOUBLE_PRECISION;
#endif
}

void comm_tag_ub(const comm_hdl ch, int *ub)
{
#ifdef MPI
  int val=0;
  int flag=0;
  MPI_Attr_get(ch,MPI_TAG_UB,&val,&flag);
  *ub = val;
#endif
}

void comm_time(double *tm)
{
#ifdef MPI
  *tm = MPI_Wtime();
#else
  *tm = 0.0;
#endif
}

void comm_barrier(const comm_hdl ch)
{
#ifdef MPI
  MPI_Barrier(ch);
#endif
}

void comm_bcast(const comm_hdl ch, void *p, size_t n, uint root)
{
#ifdef MPI
  MPI_Bcast(p,n,MPI_BYTE,root,ch);
#endif
}


uint comm_gbl_id=0, comm_gbl_np=1;

GS_DEFINE_IDENTITIES()
GS_DEFINE_DOM_SIZES()

static void scan_imp(void *scan, const struct comm *com, gs_dom dom, gs_op op,
                     const void *v, uint vn, void *buffer)
{
  comm_req req[2];
  size_t vsize = vn*gs_dom_size[dom];
  const uint id=com->id, np=com->np;
  uint n = np, c=1, odd=0, base=0;
  void *buf[2];
  void *red = (char*)scan+vsize;
  buf[0]=buffer,buf[1]=(char*)buffer+vsize;
  while(n>1) {
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
  gs_init_array(scan,vn,dom,op);
  memcpy(red,v,vsize);
  while(n<np) {
    if(c&1) n-=(odd&1), base-=n;
    c>>=1, n<<=1, n+=(odd&1);
    odd>>=1;
    if(base==id) {
      comm_irecv(&req[0],com, buf[0],vsize, id+n/2,id+n/2);
      comm_isend(&req[1],com, red   ,vsize, id+n/2,id);
      comm_wait(req,2);
      gs_gather_array(red,buf[0],vn,dom,op);
    } else {
      comm_irecv(&req[0],com, scan,vsize, base,base);
      comm_isend(&req[1],com, red ,vsize, base,id);
      comm_wait(req,2);
      break;
    }
  }
  while(n>1) {
    if(base==id) {
      comm_send(com, scan  ,2*vsize, id+n/2,id);
    } else {
      comm_recv(com, buffer,2*vsize, base,base);
      gs_gather_array(scan,buf[0],vn,dom,op);
      memcpy(red,buf[1],vsize);
    }
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
}


static void allreduce_imp(const struct comm *com, gs_dom dom, gs_op op,
                          void *v, uint vn, void *buf)
{
  size_t total_size = vn*gs_dom_size[dom];
  const uint id=com->id, np=com->np;
  uint n = np, c=1, odd=0, base=0;
  while(n>1) {
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
  while(n<np) {
    if(c&1) n-=(odd&1), base-=n;
    c>>=1, n<<=1, n+=(odd&1);
    odd>>=1;
    if(base==id) {
      comm_recv(com, buf,total_size, id+n/2,id+n/2);
      gs_gather_array(v,buf,vn, dom,op);
    } else {
      comm_send(com, v,total_size, base,id);
      break;
    }
  }
  while(n>1) {
    if(base==id)
      comm_send(com, v,total_size, id+n/2,id);
    else
      comm_recv(com, v,total_size, base,base);
    odd=(odd<<1)|(n&1);
    c<<=1, n>>=1;
    if(id>=base+n) c|=1, base+=n, n+=(odd&1);
  }
}

void comm_scan(void *scan, const struct comm *com, gs_dom dom, gs_op op,
               const void *v, uint vn, void *buffer)
{
  scan_imp(scan, com,dom,op, v,vn, buffer);
}

void comm_allreduce_cdom(const comm_hdl ch, comm_type cdom, gs_op op,
                         void *v, uint vn, void *buf)
{
  gs_dom dom;
  int dom_ok = 1;
  switch(cdom) {
    case MPI_INTEGER:          dom = gs_int; break;
    case MPI_INTEGER8:         dom = gs_long; break;
    case MPI_REAL:             dom = gs_float; break;
    case MPI_DOUBLE_PRECISION: dom = gs_double; break;            
    default: dom_ok = 0;
  }

  if (dom_ok == 1) {
    struct comm c;
    int id, np;
    c.h = ch;
    comm_id(ch,&id);
    comm_np(ch,&np);
    c.id = id;
    c.np = np;
    comm_allreduce(&c,dom,op,v,vn,buf);
  }
  else {
    fail(1,__FILE__,__LINE__,
      "comm_allreduce_cdom: cannot identify cdom P=%u.",cdom);
  }
}


void comm_allreduce(const struct comm *com, gs_dom dom, gs_op op,
                    void *v, uint vn, void *buf)
{
  if(vn==0) return;
#ifdef MPI
  {
    MPI_Datatype mpitype;
    MPI_Op mpiop;
    #define DOMAIN_SWITCH() do { \
      switch(dom) { case gs_double:    mpitype=MPI_DOUBLE;    break; \
                    case gs_float:     mpitype=MPI_FLOAT;     break; \
                    case gs_int:       mpitype=MPI_INT;       break; \
                    case gs_long:      mpitype=MPI_LONG;      break; \
     WHEN_LONG_LONG(case gs_long_long: mpitype=MPI_LONG_LONG; break;) \
                  default:        goto comm_allreduce_byhand; \
      } \
    } while(0)
    DOMAIN_SWITCH();
    #undef DOMAIN_SWITCH
    switch(op) { case gs_add: mpiop=MPI_SUM;  break;
                 case gs_mul: mpiop=MPI_PROD; break;
                 case gs_min: mpiop=MPI_MIN;  break;
                 case gs_max: mpiop=MPI_MAX;  break;
                 default:        goto comm_allreduce_byhand;
    }
    MPI_Allreduce(v,buf,vn,mpitype,mpiop,com->h);
    memcpy(v,buf,vn*gs_dom_size[dom]);
    return;
  }
#endif
#ifdef MPI
comm_allreduce_byhand:
  allreduce_imp(com,dom,op, v,vn, buf);
#endif
}

double comm_dot(const struct comm *comm, double *v, double *w, uint n)
{
  double s=tensor_dot(v,w,n),b;
  comm_allreduce(comm,gs_double,gs_add, &s,1, &b);
  return s;
}

/* T comm_reduce__T(const struct comm *comm, gs_op op, const T *in, uint n) */

#define SWITCH_OP_CASE(T,OP) case gs_##OP: WITH_OP(T,OP); break;
#define SWITCH_OP(T,op) do switch(op) { \
    GS_FOR_EACH_OP(T,SWITCH_OP_CASE) case gs_op_n: break; } while(0)

#define WITH_OP(T,OP) \
  do { T v = *in++; GS_DO_##OP(accum,v); } while(--n)

#define DEFINE_REDUCE(T) \
T PREFIXED_NAME(comm_reduce__##T)( \
    const struct comm *comm, gs_op op, const T *in, uint n) \
{                                                           \
  T accum = gs_identity_##T[op], buf;                       \
  if(n!=0) SWITCH_OP(T,op);                                 \
  comm_allreduce(comm,gs_##T,op, &accum,1, &buf);           \
  return accum;                                             \
}

GS_FOR_EACH_DOMAIN(DEFINE_REDUCE)

#undef DEFINE_REDUCE
#undef WITH_OP
#undef SWITCH_OP
#undef SWITCH_OP_CASE



/*------------------------------------------------------------------------------
  FORTRAN interface
------------------------------------------------------------------------------*/

#undef comm_init
#undef comm_finalize
#undef comm_world
#undef comm_free
#undef comm_id
#undef comm_np
#undef comm_type_int
#undef comm_type_int8
#undef comm_type_real
#undef comm_type_dp
#undef comm_tag_ub
#undef comm_time
#undef comm_barrier
#undef comm_bcast
#undef comm_allreduce_add
#undef comm_allreduce_min
#undef comm_allreduce_max
#undef comm_allreduce_mul

#define ccomm_init      PREFIXED_NAME(comm_init     )
#define ccomm_finalize  PREFIXED_NAME(comm_finalize )
#define ccomm_world     PREFIXED_NAME(comm_world)
#define ccomm_free      PREFIXED_NAME(comm_free     )
#define ccomm_id        PREFIXED_NAME(comm_id       )
#define ccomm_np        PREFIXED_NAME(comm_np       )
#define ccomm_type_int  PREFIXED_NAME(comm_type_int )
#define ccomm_type_int8 PREFIXED_NAME(comm_type_int8)
#define ccomm_type_real PREFIXED_NAME(comm_type_real)
#define ccomm_type_dp   PREFIXED_NAME(comm_type_dp  )
#define ccomm_tag_ub    PREFIXED_NAME(comm_tag_ub   )
#define ccomm_time      PREFIXED_NAME(comm_time     )
#define ccomm_barrier   PREFIXED_NAME(comm_barrier  )
#define ccomm_bcast     PREFIXED_NAME(comm_bcast    )
#define ccomm_allreduce_add PREFIXED_NAME(comm_allreduce_add)
#define ccomm_allreduce_min PREFIXED_NAME(comm_allreduce_min)
#define ccomm_allreduce_max PREFIXED_NAME(comm_allreduce_max)
#define ccomm_allreduce_mul PREFIXED_NAME(comm_allreduce_mul)

#define fcomm_init      FORTRAN_NAME(comm_init     , COMM_INIT     )
#define fcomm_finalize  FORTRAN_NAME(comm_finalize , COMM_FINALIZE )
#define fcomm_world     FORTRAN_NAME(comm_world    , COMM_WORLD    )
#define fcomm_free      FORTRAN_NAME(comm_free     , COMM_FREE     )
#define fcomm_id        FORTRAN_NAME(comm_id       , COMM_ID       )
#define fcomm_np        FORTRAN_NAME(comm_np       , COMM_NP       )
#define fcomm_type_int  FORTRAN_NAME(comm_type_int , COMM_TYPE_INT )
#define fcomm_type_int8 FORTRAN_NAME(comm_type_int8, COMM_TYPE_INT8)
#define fcomm_type_real FORTRAN_NAME(comm_type_real, COMM_TYPE_REAL)
#define fcomm_type_dp   FORTRAN_NAME(comm_type_dp  , COMM_TYPE_DP  )
#define fcomm_tag_ub    FORTRAN_NAME(comm_tag_ub   , COMM_TAG_UB   )
#define fcomm_time      FORTRAN_NAME(comm_time     , COMM_TIME     )
#define fcomm_barrier   FORTRAN_NAME(comm_barrier  , COMM_BARRIER  )
#define fcomm_bcast     FORTRAN_NAME(comm_bcast    , COMM_BCAST    )
#define fcomm_allreduce_add FORTRAN_NAME(comm_allreduce_add, COMM_ALLREDUCE_ADD)
#define fcomm_allreduce_min FORTRAN_NAME(comm_allreduce_min, COMM_ALLREDUCE_MIN)
#define fcomm_allreduce_max FORTRAN_NAME(comm_allreduce_max, COMM_ALLREDUCE_MAX)
#define fcomm_allreduce_mul FORTRAN_NAME(comm_allreduce_mul, COMM_ALLREDUCE_PROD)


void fcomm_init(void)
{
  ccomm_init();
}

void fcomm_finalize(void)
{
  ccomm_finalize();
}


void fcomm_world(comm_hdl *chp)
{
  ccomm_world(chp);
}

void fcomm_free(comm_hdl *chp)
{
  ccomm_free(*chp);
}


void fcomm_id(const comm_hdl *chp, int *id)
{
  ccomm_id(*chp, id);
}

void fcomm_np(const comm_hdl *chp, int *np)
{
  ccomm_np(*chp, np);
}


void fcomm_type_int(comm_type *ct)
{
  ccomm_type_int(ct);
}

void fcomm_type_int8(comm_type *ct)
{
  ccomm_type_int8(ct);
}

void fcomm_type_real(comm_type *ct)
{
  ccomm_type_real(ct);
}

void fcomm_type_dp(comm_type *ct)
{
  ccomm_type_dp(ct);
}

void fcomm_tag_ub(const comm_hdl *chp, int *ub)
{
  ccomm_tag_ub(*chp, ub);
}

void fcomm_time(double *tm)
{
  ccomm_time(tm);
}


void fcomm_barrier(const comm_hdl *chp)
{
  ccomm_barrier(*chp);
}

void fcomm_bcast(const comm_hdl *chp, void *p, size_t *n, uint *root)
{
  ccomm_bcast(*chp,p,*n,*root);
}


void fcomm_allreduce_add(const comm_hdl *chp, comm_type *ct,
                         void *v, uint *vn, void *buf)
{
  comm_allreduce_cdom(*chp,*ct,gs_add,v,*vn,buf);
}

void fcomm_allreduce_min(const comm_hdl *chp, comm_type *ct,
                         void *v, uint *vn, void *buf)
{
  comm_allreduce_cdom(*chp,*ct,gs_min,v,*vn,buf);
}

void fcomm_allreduce_max(const comm_hdl *chp, comm_type *ct,
                         void *v, uint *vn, void *buf)
{
  comm_allreduce_cdom(*chp,*ct,gs_max,v,*vn,buf);
}

void fcomm_allreduce_mul(const comm_hdl *chp, comm_type *ct,
                         void *v, uint *vn, void *buf)
{
  comm_allreduce_cdom(*chp,*ct,gs_mul,v,*vn,buf);
}



