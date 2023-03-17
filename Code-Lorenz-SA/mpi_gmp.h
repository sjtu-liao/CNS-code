/**********************************************/
/* mpi_gmp.h:                                 */
/* Copyright (C) 2003 Tomonori Kouya          */
/*                                            */
/* Version 0.0: 2003.04/01                    */
/* Version 0.1: 2003.05/26 append _GET_PREC macro */
/*                                            */
/* This library is free software; you can re- */
/* distribute it and/or modify it under the   */
/* terms of the GNU Lesser General Public     */
/* License as published by the Free Software  */
/* Foundation; either version 2.1 of the      */
/* License, or (at your option) any later     */
/* version.                                   */
/*                                            */
/* This library is distributed in the hope    */
/* that it will be useful, but WITHOUT ANY    */
/* WARRANTY; without even the implied         */
/* warranty of MERCHANTABILITY or FITNESS FOR */
/* A PARTICULAR PURPOSE.  See the GNU Lesser  */
/* General Public License for more details.   */
/**********************************************/
#ifndef __MPI_GMP_H
#define __MPI_GMP_H
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmp.h"


#ifndef __MPI_INCLUDE
	#include "mpi.h"
#endif

#ifdef USE_MPFR
	#include "mpfr.h"
	#include "mpf2mpfr.h"
#endif

#ifdef __MPFR_H
#if MPFR_VERSION_MAJOR >= 2
  #if MPFR_VERSION_MINOR >= 1
  //#define _NUM_LIMB(a) (mpfr_get_prec((mpfr_ptr)(a)))
  #define _NUM_LIMB(a) ((a->_mpfr_prec - 1) / GMP_LIMB_BITS + 1)
  //#define _NUM_LIMB(a) ((unsigned long)ceil((double)((a->_mpfr_prec) / GMP_LIMB_BITS)))
  #else
    #define _NUM_LIMB(a) ((a->_mpfr_prec / GMP_LIMB_BITS) + 1)
  #endif
#else
  #define _NUM_LIMB(a) ((a->_mpfr_prec / GMP_LIMB_BITS) + 1)
#endif
#define _GET_PREC(a_prec) (a_prec)
#else
#define _GET_PREC(a_prec) ((a_prec-1) * GMP_LIMB_BITS)
#endif

#define MPI_GMP_MAXPROCS 128 /* because of rsh! */

#if __DEBUGMODE
#define MPF_PRINT(str, a) {printf("%s", str); mpf_out_str(stdout, 10, 0, a);printf("\n");fflush(stdout);}
#else
#define MPF_PRINT(str, a)
#endif

/* Datatype of mpf_t */
MPI_Datatype gmp_mpf_array[128];
MPI_Datatype gmp_mpf;
MPI_Datatype gmp_mpf_b128;	/* 38 decimal digits */
MPI_Datatype gmp_mpf_b256;	/* 77 decimal digts */
MPI_Datatype gmp_mpf_b512;	/* 154 decimal digits */
MPI_Datatype gmp_mpf_b1024;	/* 308 decimal digits */
MPI_Datatype gmp_mpf_b2048;	/* 616 decimal digits */
MPI_Datatype gmp_mpf_b4096;	/* 1233 decimal digits */
MPI_Datatype gmp_mpf_b8192;	/* 2466 decimal digits */
MPI_Datatype gmp_mpf_b16384;	/* 4932 decimal digits */
MPI_Datatype gmp_mpf_b32768;	/* 9864 decimal digits */
MPI_Datatype gmp_mpf_b65536;	/* 19728 decimal digits */
MPI_Datatype gmp_mpf_b131072;	/* 39456 decimal digits */

MPI_Datatype gmp_mpf_b167;	/* 50 decimal digits */
MPI_Datatype gmp_mpf_b333;	/* 100 decimal digits */
MPI_Datatype gmp_mpf_b667;	/* 200 decimal digits */
MPI_Datatype gmp_mpf_b1332;	/* 400 decimal digits */
MPI_Datatype gmp_mpf_b1661;	/* 500 decimal digits */
MPI_Datatype gmp_mpf_b3322;	/* 1000 decimal digits */
MPI_Datatype gmp_mpf_b33220;	/* 10000 decimal digits */
MPI_Datatype gmp_mpf_b332193;	/* 100000 decimal digits */

#define MPI_MPF gmp_mpf
#define MPI_MPF_B128 gmp_mpf_b128
#define MPI_MPF_B256 gmp_mpf_b256
#define MPI_MPF_B512 gmp_mpf_b512
#define MPI_MPF_B1024 gmp_mpf_b1024
#define MPI_MPF_B2048 gmp_mpf_b2048
#define MPI_MPF_B4096 gmp_mpf_b4096
#define MPI_MPF_B8192 gmp_mpf_b8192
#define MPI_MPF_B16384 gmp_mpf_b16384
#define MPI_MPF_B32768 gmp_mpf_b32768
#define MPI_MPF_B65536 gmp_mpf_b65536
#define MPI_MPF_B131072 gmp_mpf_b131072

#define MPI_MPF_D38 gmp_mpf_b128
#define MPI_MPF_D77 gmp_mpf_b256
#define MPI_MPF_D154 gmp_mpf_b512
#define MPI_MPF_D308 gmp_mpf_b1024
#define MPI_MPF_D616 gmp_mpf_b2048
#define MPI_MPF_D1233 gmp_mpf_b4096
#define MPI_MPF_D2466 gmp_mpf_b8192
#define MPI_MPF_D4932 gmp_mpf_b16384
#define MPI_MPF_D9864 gmp_mpf_b32768
#define MPI_MPF_D19728 gmp_mpf_b65536
#define MPI_MPF_D39456 gmp_mpf_b131072

#define MPF_D50 (167)
#define MPF_D100 (333)
#define MPF_D200 (667)
#define MPF_D400 (1332)
#define MPF_D500 (1661)
#define MPF_D1000 (3322)
#define MPF_D10000 (33220)
#define MPF_D100000 (332193)

#define MPI_MPF_D50 gmp_mpf_b167
#define MPI_MPF_D100 gmp_mpf_b333
#define MPI_MPF_D200 gmp_mpf_b667
#define MPI_MPF_D400 gmp_mpf_b1332
#define MPI_MPF_D500 gmp_mpf_b1661
#define MPI_MPF_D1000 gmp_mpf_b3322
#define MPI_MPF_D10000 gmp_mpf_b33220
#define MPI_MPF_D100000 gmp_mpf_b332193

/* Operation for mpf_ts */
MPI_Op gmp_mpf_add;
MPI_Op gmp_mpf2_add;
#define MPI_MPF_SUM gmp_mpf_add
#define MPI_MPF2_SUM gmp_mpf2_add

/* divide number of dimension */
long int _mpi_divide_dim(long int d_dim[], long int dim, int num_procs);

/* get bufsize for mpf_t */
size_t get_bufsize_mpf(mpf_t, int);
size_t get_bufsize_mpf2(mpf_t, int);
size_t get_bufsize_mpz(mpz_t);
size_t get_bufsize_mpq(mpq_t);

/* allocate buf for mpf_t */
void *allocbuf_mpf(unsigned long prec, int incount);
void *allocbuf_mpf2(unsigned long prec, int incount);
void *allocbuf_mpz(mpz_t);
void *allocbuf_mpq(mpq_t);

/* free buf */
void freebuf_mpf(void *);
void freebuf_mpf2(void *);
void freebuf_mpz(void *);
void freebuf_mpq(void *);

/* pack mpf */
#ifdef __USE_MPI_PACK
void pack_mpf(mpf_t a, int incount, void *buf, int *bufsize, int *pos, MPI_Comm comm);
#else
void pack_mpf(mpf_t a, int incount, void *buf);
#endif
void pack_mpz(mpz_t a, void *buf);
void pack_mpq(mpq_t a, void *buf);

/* unpack mpf */
#ifdef __USE_MPI_PACK
void unpack_mpf(void *buf, int bufsize, int *pos, mpf_t ret, int count, MPI_Comm comm);
#else
void unpack_mpf(void *buf, mpf_t ret, int count);
#endif
void unpack_mpz(void *buf, mpz_t ret);
void unpack_mpq(void *buf, mpq_t ret);

void pack_mpf2(mpf_t a, int incount, void *buf, unsigned long prec);
void unpack_mpf2(void *buf, mpf_t ret, int count, unsigned long *prec);

/* typedef and commit to mpich */
void commit_mpf(MPI_Datatype *mpi_mpf_t, unsigned long prec, MPI_Comm comm);

/* create op */
void create_mpf_op(MPI_Op *mpi_mpf_op, void (*func)(void *, void *, int *, MPI_Datatype *), MPI_Comm comm);

/* clear type */
void free_mpf(MPI_Datatype *mpi_mpf_t);

/* clear op */
void free_mpf_op(MPI_Op *mpi_mpf_op);

/* macros at each precision */
void use_mpi_mpf_b128(MPI_Comm comm);
void use_mpi_mpf_b256(MPI_Comm comm);
void use_mpi_mpf_b512(MPI_Comm comm);
void use_mpi_mpf_b1024(MPI_Comm comm);
void use_mpi_mpf_b2048(MPI_Comm comm);
void use_mpi_mpf_b4096(MPI_Comm comm);
void use_mpi_mpf_b8192(MPI_Comm comm);
void use_mpi_mpf_b16384(MPI_Comm comm);
void use_mpi_mpf_b32768(MPI_Comm comm);
void use_mpi_mpf_b65536(MPI_Comm comm);
void use_mpi_mpf_b131072(MPI_Comm comm);
void use_mpi_mpf_d50(MPI_Comm comm);
void use_mpi_mpf_d100(MPI_Comm comm);
void use_mpi_mpf_d1000(MPI_Comm comm);
void use_mpi_mpf_d10000(MPI_Comm comm);
void use_mpi_mpf_d100000(MPI_Comm comm);
void free_mpi_mpf_b128(void);
void free_mpi_mpf_b256(void);
void free_mpi_mpf_b512(void);
void free_mpi_mpf_b1024(void);
void free_mpi_mpf_b2048(void);
void free_mpi_mpf_b4096(void);
void free_mpi_mpf_b8192(void);
void free_mpi_mpf_b16384(void);
void free_mpi_mpf_b32768(void);
void free_mpi_mpf_b65536(void);
void free_mpi_mpf_b131072(void);
void free_mpi_mpf_d50(void);
void free_mpi_mpf_d100(void);
void free_mpi_mpf_d1000(void);
void free_mpi_mpf_d10000(void);
void free_mpi_mpf_d100000(void);

/* Operations for mpf_t */

/* mpf_add for MPI */
void _mpi_mpf_add(void *in, void *ret, int *len, MPI_Datatype *datatype);
void _mpi_mpf2_add(void *in, void *ret, int *len, MPI_Datatype *datatype);
