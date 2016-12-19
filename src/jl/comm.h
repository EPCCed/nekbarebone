#ifndef COMM_H
#define COMM_H

/* requires:
     <stddef.h>            for size_t
     <stdlib.h>            for exit
     "fail.h", "types.h"
     "gs_defs.h"           for comm_allreduce, comm_scan, comm_reduce_T
*/

#if !defined(FAIL_H) || !defined(TYPES_H)
#warning "comm.h" requires "fail.h" and "types.h"
#endif

/*
  When the preprocessor macro MPI is defined, defines (very) thin wrappers
  for the handful of used MPI routines. Alternatively, when MPI is not defined,
  these wrappers become dummy routines suitable for a single process run.
  No code outside of "comm.h" and "comm.c" makes use of MPI at all.

  Basic usage:
  
    struct comm c;
  
    comm_init(&c, MPI_COMM_WORLD);  // initializes c using MPI_Comm_dup

    comm_free(&c);
  
  Very thin MPI wrappers: (see below for implementation)

    comm_send,_recv,_isend,_irecv,_time,_barrier
    
  Additionally, some reduction and scan routines are provided making use
    of the definitions in "gs_defs.h" (provided this has been included first).

  Example comm_allreduce usage:
    
    double v[5], buf[5];
    comm_allreduce(&c, gs_double,gs_add, v,5,buf);
      // Computes the vector sum of v across all procs, using
      // buf as a scratch area. Delegates to MPI_Allreduce if possible.
    
  Example comm_scan usage:
    
    long in[5], out[2][5], buf[2][5];
    comm_scan(out, &c,gs_long,gs_add, in,5,buf);
      // out[0] will be the vector sum of "in" across procs with ids
           *strictly* less than this one (exclusive behavior),
         and out[1] will be the vector sum across all procs, as would
           be computed with comm_allreduce.
         Note: differs from MPI_Scan which has inclusive behavior
  
  Example comm_reduce_double, etc. usage:
  
    T out, in[10];
    out = comm_reduce_T(&c, gs_max, in, 10);
      // out will equal the largest element of "in",
         across all processors
      // T can be "double", "float", "int", "long", "slong", "sint", etc.
         as defined in "gs_defs.h"
         
*/

typedef struct comm * comm_ptr;

#ifdef MPI
#include <mpi.h>
typedef MPI_Comm comm_hdl;
typedef MPI_Datatype comm_type;
typedef MPI_Op comm_op;
typedef MPI_Request comm_req;
#else
typedef int comm_hdl;
typedef int comm_type;
typedef int comm_op;
typedef int comm_req;
#endif


#define comm_init      PREFIXED_NAME(comm_init     )
#define comm_finalize  PREFIXED_NAME(comm_finalize )

#define comm_world     PREFIXED_NAME(comm_world    )
#define comm_dup       PREFIXED_NAME(comm_dup      )
#define comm_free      PREFIXED_NAME(comm_free     )

#define comm_np        PREFIXED_NAME(comm_np       )
#define comm_id        PREFIXED_NAME(comm_id       )

#define comm_type_int  PREFIXED_NAME(comm_type_int )
#define comm_type_int8 PREFIXED_NAME(comm_type_int8)
#define comm_type_real PREFIXED_NAME(comm_type_real)
#define comm_type_dp   PREFIXED_NAME(comm_type_dp  )
#define comm_tag_ub    PREFIXED_NAME(comm_tag_ub   )
#define comm_time      PREFIXED_NAME(comm_time     )

#define comm_barrier   PREFIXED_NAME(comm_barrier  )
#define comm_bcast     PREFIXED_NAME(comm_bcast    )

#define comm_allreduce_cdom PREFIXED_NAME(comm_allreduce_cdom)
#define comm_allreduce PREFIXED_NAME(comm_allreduce)
#define comm_scan      PREFIXED_NAME(comm_scan     )
#define comm_dot       PREFIXED_NAME(comm_dot      )

/* global id, np vars strictly for diagnostic messages (fail.c) */
#ifndef comm_gbl_id
#define comm_gbl_id PREFIXED_NAME(comm_gbl_id)
#define comm_gbl_np PREFIXED_NAME(comm_gbl_np)
extern uint comm_gbl_id, comm_gbl_np;
#endif


struct comm {
  comm_hdl h;
  int np, id;
#ifdef __UPC__
  shared[] char *shared *buf_dir; /* Global directory of buffers */
  shared int *flgs;
  char *buf;			  /* Local part of buffers */
  size_t buf_len;		  /* Shared buffer size */
#endif
};


void comm_init(void);
void comm_finalize(void);

void comm_world(comm_ptr *cpp);
void comm_dup(comm_ptr *cpp, const comm_ptr cp);
void comm_free(comm_ptr *cpp);

void comm_np(const comm_ptr cp, int *np);
void comm_id(const comm_ptr cp, int *id);

void comm_type_int(comm_type *ct);
void comm_type_int8(comm_type *ct);
void comm_type_real(comm_type *ct);
void comm_type_dp(comm_type *ct);
void comm_tag_ub(const comm_ptr cp, int *ub);
void comm_time(double *tm);

void comm_barrier(const comm_ptr cp);
void comm_bcast(const comm_ptr cp, void *p, size_t n, uint root);


#ifdef GS_DEFS_H
void comm_allreduce_cdom(const comm_ptr cp, comm_type cdom, gs_op op,
                         void *v, uint vn, void *buf);

void comm_allreduce(const comm_ptr cp, gs_dom dom, gs_op op,
                    void *v, uint vn, void *buf);

void comm_scan(void *scan, const comm_ptr cp, gs_dom dom, gs_op op,
               const void *v, uint vn, void *buffer);

double comm_dot(const comm_ptr cp, double *v, double *w, uint n);

#define DEFINE_REDUCE(T) \
T PREFIXED_NAME(comm_reduce__##T)( \
    const comm_ptr cp, gs_op op, const T *in, uint n); \
static T comm_reduce_##T(const comm_ptr cp, gs_op op, const T *v, uint vn) \
{ return PREFIXED_NAME(comm_reduce__##T)(cp,op,v,vn); }
GS_FOR_EACH_DOMAIN(DEFINE_REDUCE)
#undef DEFINE_REDUCE

#define comm_reduce_sint \
    TYPE_LOCAL(comm_reduce_int,comm_reduce_long,comm_reduce_long_long)
#define comm_reduce_slong \
   TYPE_GLOBAL(comm_reduce_int,comm_reduce_long,comm_reduce_long_long)

#endif

static void comm_recv(const comm_ptr cp, void *p, size_t n,
                      uint src, int tag)
{
#ifdef MPI
# ifndef MPI_STATUS_IGNORE
  MPI_Status stat;
  MPI_Recv(p,n,MPI_UNSIGNED_CHAR,src,tag,cp->h,&stat);
# else  
  MPI_Recv(p,n,MPI_UNSIGNED_CHAR,src,tag,cp->h,MPI_STATUS_IGNORE);
# endif
#endif
}

static void comm_send(const comm_ptr cp, void *p, size_t n,
                      uint dst, int tag)
{
#ifdef MPI
  MPI_Send(p,n,MPI_UNSIGNED_CHAR,dst,tag,cp->h);
#endif
}

static void comm_irecv(comm_req *req, const comm_ptr cp,
                       void *p, size_t n, uint src, int tag)
{
#ifdef MPI
  MPI_Irecv(p,n,MPI_UNSIGNED_CHAR,src,tag,cp->h,req);
#endif
}

static void comm_isend(comm_req *req, const comm_ptr cp,
                       void *p, size_t n, uint dst, int tag)
{
#ifdef MPI
  MPI_Isend(p,n,MPI_UNSIGNED_CHAR,dst,tag,cp->h,req);
#endif
}

static void comm_wait(comm_req *req, int n)
{
#ifdef MPI
# ifndef MPI_STATUSES_IGNORE
  MPI_Status status[8];
  while(n>=8) MPI_Waitall(8,req,status), req+=8, n-=8;
  if(n>0) MPI_Waitall(n,req,status);
# else
  MPI_Waitall(n,req,MPI_STATUSES_IGNORE);
# endif  
#endif
}

#endif
