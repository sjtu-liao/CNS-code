#include "mpi_gmp.h"
char *readme="write by wpf IAP/LASG 2005.3.11, add taylor method 2010.3\n";
size_t get_bufsize_mpf(mpf_t a, int incount)
{
	size_t bufsize;

#ifdef __MPFR_H
	/* Definition of the main structure */
	//
	//	typedef struct {
	//	  mpfr_prec_t  _mpfr_prec;
	//	  mpfr_sign_t  _mpfr_sign;
	//	  mp_exp_t     _mpfr_exp;
	//	  mp_limb_t   *_mpfr_d;
	//	} __mpfr_struct;
	//	typedef __mpfr_struct mpfr_t[1];
	//	typedef __mpfr_struct *mpfr_ptr;
	//	typedef __gmp_const __mpfr_struct *mpfr_srcptr;
	//	above def. from mpfr.h in mpfr-2.2.0
	//
	//bufsize = sizeof(mpfr_sign_t) + sizeof(mpfr_prec_t) + sizeof(mp_exp_t) + sizeof(mp_limb_t) * _NUM_LIMB(a);
#if MPFR_VERSION_MAJOR >= 2	
#if MPFR_VERSION_MINOR >= 1
	bufsize = sizeof(mpfr_sign_t) + sizeof(mpfr_prec_t);
#else	
	bufsize = sizeof(mpfr_size_t) + sizeof(mpfr_prec_t);
#endif
#else
	bufsize = sizeof(mp_size_t) + sizeof(mp_prec_t);
#endif
	bufsize += sizeof(mp_exp_t);
	bufsize += sizeof(mp_limb_t) * _NUM_LIMB(a);
	//bufsize = sizeof(int) * 2 + sizeof(mp_exp_t) + sizeof(mp_limb_t) * _NUM_LIMB(a);
	// typedef struct
	// {
	//		int _mp_prec;
	//		int _mp_size;
	//		mp_exp_t _mp_exp;		 Exponent, in the base of `mp_limb_t'.  */
	//		mp_limb_t *_mp_d;		 Pointer to the limbs.  */
	//	} __mpf_struct;
	// typedef __mpf_struct MP_FLOAT; 
	// typedef __mpf_struct mpf_t[1];
	// above def. from gmp.h in gmp-4.1.4
	//
#else
	bufsize = sizeof(int) * 2 + sizeof(mp_exp_t) + sizeof(mp_limb_t) * (a->_mp_prec + 1);
#endif

	return (size_t)(bufsize * incount);
}
/* mpf_add for MPI */
void _mpi_mpf_add(void *in, void *ret, int *len, MPI_Datatype *datatype)
{
	int i, itmp;
	unsigned long prec;
	void *ptr_in, *ptr_ret;
	mpf_t tmp_ret, tmp_in;
	// printf("use _mpf_add\n");
	
	ptr_in = (void *)in;
	memcpy(&itmp, ptr_in, sizeof(int));
	prec = _GET_PREC(itmp);


	if(*datatype == MPI_MPF)
	{
		mpf_init(tmp_ret);
		mpf_init(tmp_in);
	}
	else
	{	
		mpf_init2(tmp_ret, prec);
		mpf_init2(tmp_in, prec);
	}

	ptr_in = (void *)in;
	ptr_ret = (void *)ret;

	for(i = 0; i < *len; i++)
	{
		/* unpack */
		unpack_mpf(ptr_ret, tmp_ret, 1);
		unpack_mpf(ptr_in, tmp_in, 1);

		/* add */
		mpf_add(tmp_ret, tmp_in, tmp_ret);

		/* pack */
		pack_mpf(tmp_ret, 1, (void *)ptr_ret);

		ptr_ret += get_bufsize_mpf(tmp_ret, 1);
		ptr_in += get_bufsize_mpf(tmp_in, 1);
	}
	//	printf("add over\n");

	mpf_clear(tmp_in);
	mpf_clear(tmp_ret);
}

/* unpack mpf */
#ifdef __USE_MPI_PACK
void unpack_mpf(void *buf, int bufsize, int *pos, mpf_t ret, int count, MPI_Comm comm)
#else
void unpack_mpf(void *buf, mpf_t ret, int count)
#endif
{
	unsigned long prec;
	int i;
	//unsigned char *tmp_buf;
	unsigned char *tmp_buf;
	//mpf_ptr ptr_ret;
	unsigned char *ptr_ret;
	//mpf_ptr tmp_ret;
	mpf_ptr tmp_ret;
	mp_limb_t *limb_d;

	//tmp_buf = (unsigned char *)buf;
	tmp_buf = (unsigned char *)buf;
	//ptr_ret = (mpf_ptr)ret;
	ptr_ret = (unsigned char *)ret;
	tmp_ret = (mpf_ptr)ret;

	for(i = 0; i < count; i++)
	{
#ifndef __USE_MPI_PACK
#ifdef __MPFR_H
#if MPFR_VERSION_MAJOR >= 2
#if MPFR_VERSION_MINOR >= 1
		memcpy(&(ret->_mpfr_sign), tmp_buf, sizeof(mpfr_sign_t));
		tmp_buf += sizeof(mpfr_sign_t); ptr_ret += sizeof(mpfr_sign_t);
		memcpy(&(ret->_mpfr_prec), tmp_buf, sizeof(mpfr_prec_t));
		tmp_buf += sizeof(mpfr_prec_t); ptr_ret += sizeof(mpfr_prec_t);
#else
		memcpy(&(ret->_mpfr_size), tmp_buf, sizeof(int));
		tmp_buf += sizeof(int); ptr_ret += sizeof(int);
		memcpy(&(ret->_mpfr_prec), tmp_buf, sizeof(mpfr_prec_t));
		tmp_buf += sizeof(mpfr_prec_t); ptr_ret += sizeof(mpfr_prec_t);
#endif
#else
		memcpy(&(ret->_mpfr_size), tmp_buf, sizeof(int));
		tmp_buf += sizeof(int); ptr_ret += sizeof(int);
		memcpy(&(ret->_mpfr_prec), tmp_buf, sizeof(mp_prec_t));
		tmp_buf += sizeof(mp_prec_t); ptr_ret += sizeof(mp_prec_t);
#endif
		memcpy(&(ret->_mpfr_exp), tmp_buf, sizeof(mp_exp_t));
		tmp_buf += sizeof(mp_exp_t); ptr_ret += sizeof(mp_exp_t);

		memcpy(ret->_mpfr_d, tmp_buf, (size_t)(sizeof(mp_limb_t) * _NUM_LIMB(ret)));
		tmp_buf += sizeof(mp_limb_t) * _NUM_LIMB(ret); ptr_ret += sizeof(mp_limb_t *);
#else
		memcpy(&(ret->_mp_size), tmp_buf, sizeof(int));
		tmp_buf += sizeof(int); ptr_ret += sizeof(int);
		memcpy(&(ret->_mp_prec), tmp_buf, sizeof(int));
		tmp_buf += sizeof(int); ptr_ret += sizeof(int);
		memcpy(&(ret->_mp_exp), tmp_buf, sizeof(mp_exp_t));
		tmp_buf += sizeof(mp_exp_t); ptr_ret += sizeof(mp_exp_t);
		memcpy(ret->_mp_d, tmp_buf, sizeof(mp_limb_t) * (ret->_mp_prec + 1));
		tmp_buf += sizeof(mp_limb_t) * (ret->_mp_prec + 1); ptr_ret += sizeof(mp_limb_t *);
#endif
#else
#ifdef __MPFR_H
		MPI_Unpack(buf, bufsize, pos, &ret->_mpfr_size, 1, MPI_INT, comm); 
		MPI_Unpack(buf, bufsize, pos, &ret->_mpfr_prec, 1, MPI_INT, comm); 
		MPI_Unpack(buf, bufsize, pos, &ret->_mpfr_exp, 1, MPI_LONG, comm); 
		MPI_Unpack(buf, bufsize, pos, ret->_mpfr_d, _NUM_LIMB(ret), MPI_UNSIGNED_LONG, comm);
#else
		MPI_Unpack(buf, bufsize, pos, &ret->_mp_size, 1, MPI_INT, comm); 
		MPI_Unpack(buf, bufsize, pos, &ret->_mp_prec, 1, MPI_INT, comm); 
		MPI_Unpack(buf, bufsize, pos, &ret->_mp_exp, 1, MPI_LONG, comm); 
		MPI_Unpack(buf, bufsize, pos, ret->_mp_d, ret->_mp_prec + 1, MPI_UNSIGNED_LONG, comm);
#endif
#endif
		//		(unsigned char *)ret = ptr_ret;
		//ret = (unsigned char *)ptr_ret;
		//printf("unpack_mpf->"); mpf_out_str(stdout, 10, 0, ret); printf("\n");
		ret++;
	}
#if __GNUC__ >= 4
	ret = tmp_ret;
#else
	(mpf_ptr)ret = tmp_ret;
#endif
}

/* create op */
void create_mpf_op(MPI_Op *mpi_mpf_op, void (*func)(void *, void *, int *, MPI_Datatype *), MPI_Comm comm)
{
	/* operation create */
//	MPI_Op_create(func, MPI_COMM_WORLD, mpi_mpf_op);
	MPI_Op_create(func, comm, mpi_mpf_op);
#ifdef __MPFR_H
	mpfr_set_default_rounding_mode(GMP_RNDN);
#endif
}


/* clear type */
void free_mpf(MPI_Datatype *mpi_mpf_t)
{
	MPI_Type_free(mpi_mpf_t);
}

/* clear op */
void free_mpf_op(MPI_Op *mpi_mpf_op)
{
	MPI_Op_free(mpi_mpf_op);
}




/* typedef and commit to mpich */
void commit_mpf(MPI_Datatype *mpi_mpf_t, unsigned long prec, MPI_Comm comm)
{
	int i, blockcounts[5];
	MPI_Datatype types[5];
	MPI_Aint displacements[5];

	mpf_t a;
	void *buf;
	int pos, bufsize;

	mpf_init2(a, prec);
	mpf_set_ui(a, 2);
	mpf_sqrt(a, a);

	blockcounts[0] = 1;
#if MPFR_VERSION_MAJOR >= 2
	/* struct of mpfr */
	/*   mpfr_prec_t(long) _mpfr_prec */
	/*   mpfr_sign_t(int)  _mpfr_sign */
	/*   mp_exp_t(long)  _mpfr_exp */
	/*   mp_limb_t(unsigned long) *_mpfr_d */
	types[0] = MPI_LONG;
#else
	types[0] = MPI_INT;
#endif

	blockcounts[1] = 1;
	types[1] = MPI_INT;

	blockcounts[2] = 1;
	types[2] = MPI_LONG;

#ifdef __MPFR_H
	blockcounts[3] = _NUM_LIMB(a);
#else
	blockcounts[3] = a->_mp_prec + 1;
#endif
	types[3] = MPI_UNSIGNED_LONG;


	blockcounts[4] = 1;
	types[4] = MPI_UB;

	/* Pack */
	pos = 0;
	//buf = NULL;
	buf = allocbuf_mpf(prec, 1);
	pack_mpf(a, 1, buf);

	MPI_Address(buf, &displacements[0]);
	MPI_Address(buf + sizeof(long), &displacements[1]);
	MPI_Address(buf + sizeof(long) + sizeof(int), &displacements[2]);
	MPI_Address(buf + sizeof(long) + sizeof(int) + sizeof(mp_exp_t), &displacements[3]);
	displacements[4] = get_bufsize_mpf(a, 1);
	displacements[3] -= displacements[0];
	displacements[2] -= displacements[0];
	displacements[1] -= displacements[0];
	displacements[0] = 0;

	/* print typemap */
//	for(i = 0; i < 4; i++)
//		printf("%d: (%d, %d) %d\n", i, types[i], blockcounts[i], displacements[i]);

	MPI_Type_struct(5, blockcounts, displacements, types, mpi_mpf_t);
	MPI_Type_commit(mpi_mpf_t);

	/* get size */
	MPI_Pack_size(1, *mpi_mpf_t, comm, &pos);
	//printf("size = %d\n", pos);

	mpf_clear(a);
	free(buf);
}

/* allocate buf for mpf_t                     */
/* |<---------get_bufsize_mpf()------------>| */
/* +--------+--------+---+------------------+ */
/* |mpf_t[0]|mpf_t[1]|...|mpf_t[incount - 1]| */
/* +--------+--------+---+------------------+ */
void *allocbuf_mpf(unsigned long prec, int incount)
{
	void *buf;
	int bufsize;
	mpf_t a;
	//	printf("use allocbuf");
	mpf_init2(a, prec);
	
	buf = (void *)malloc(get_bufsize_mpf(a, incount));

	mpf_clear(a);

	return buf;
}

/* pack mpf */
#ifdef __USE_MPI_PACK
void pack_mpf(mpf_t a, int incount, void *buf, int *bufsize, int *pos, MPI_Comm comm)
#else
void pack_mpf(mpf_t a, int incount, void *buf)
#endif
{
	int i;
	unsigned long prec;
	unsigned char *tmp_buf;
	unsigned char *ptr_a;
	mpf_ptr tmp_a;
	mp_limb_t *limb_d;

	prec = mpf_get_prec(a);

	if(buf == NULL)
	{
		buf = allocbuf_mpf(prec, incount);
	}

	tmp_buf = (unsigned char *)buf;
	ptr_a = (unsigned char *)a;
	tmp_a = (mpf_ptr)a;

	for(i = 0; i < incount; i++)
	{
	  //   printf("pack_mpf->"); mpf_out_str(stdout, 10, 0, a); printf("\n");
/* 2. int a->_mp_size; */
/* 3. int a->_mp_prec; */
/* 4. mp_exp_t(long?) a->_mp_exp; */
/* 5. mp_limb_t(unsigned long?) a->_mp_limb; */

#ifndef __USE_MPI_PACK
#ifdef __MPFR_H
#if MPFR_VERSION_MAJOR >= 2
#if MPFR_VERSION_MINOR >= 1
		memcpy(tmp_buf, &(a->_mpfr_sign), sizeof(mpfr_sign_t));
		tmp_buf += sizeof(mpfr_sign_t); ptr_a += sizeof(mpfr_sign_t);
		memcpy(tmp_buf, &(a->_mpfr_prec), sizeof(mpfr_prec_t));
		tmp_buf += sizeof(mpfr_prec_t); ptr_a += sizeof(mpfr_prec_t);
#else
		memcpy(tmp_buf, &(a->_mpfr_size), sizeof(int));
		tmp_buf += sizeof(int); ptr_a += sizeof(int);
		memcpy(tmp_buf, &(a->_mpfr_prec), sizeof(mpfr_prec_t));
		tmp_buf += sizeof(mpfr_prec_t); ptr_a += sizeof(mpfr_prec_t);
#endif
#else
		memcpy(tmp_buf, &(a->_mpfr_size), sizeof(int));
		tmp_buf += sizeof(int); ptr_a += sizeof(int);
		memcpy(tmp_buf, &(a->_mpfr_prec), sizeof(mp_prec_t));
		tmp_buf += sizeof(mp_prec_t); ptr_a += sizeof(mp_prec_t);
#endif
		memcpy(tmp_buf, &(a->_mpfr_exp), sizeof(mp_exp_t));
		tmp_buf += sizeof(mp_exp_t); ptr_a += sizeof(mp_exp_t);
		memcpy(tmp_buf, a->_mpfr_d, (size_t)(sizeof(mp_limb_t) * _NUM_LIMB(a)));
		tmp_buf += sizeof(mp_limb_t) * _NUM_LIMB(a); ptr_a += sizeof(mp_limb_t *);
#else
		memcpy(tmp_buf, &(a->_mp_size), sizeof(int));
		tmp_buf += sizeof(int); ptr_a += sizeof(int);
		memcpy(tmp_buf, &(a->_mp_prec), sizeof(int));
		tmp_buf += sizeof(int); ptr_a += sizeof(int);
		memcpy(tmp_buf, &(a->_mp_exp), sizeof(mp_exp_t));
		tmp_buf += sizeof(mp_exp_t); ptr_a += sizeof(mp_exp_t);
		memcpy(tmp_buf, a->_mp_d, (size_t)(sizeof(mp_limb_t) * (a->_mp_prec + 1))); 
		tmp_buf += sizeof(mp_limb_t) * (a->_mp_prec + 1); ptr_a += sizeof(mp_limb_t *);
#endif
#else
#ifdef __MPFR_H
		MPI_Pack(&(a->_mpfr_size), 1, MPI_INT, buf, *bufsize, pos, comm);
		MPI_Pack(&(a->_mpfr_prec), 1, MPI_INT, buf, *bufsize, pos, comm);
		MPI_Pack(&(a->_mpfr_exp), 1, MPI_LONG, buf, *bufsize, pos, comm);
		MPI_Pack(a->_mpfr_d, _NUM_LIMB(a), MPI_UNSIGNED_LONG, buf, *bufsize, pos, comm);
#else
		MPI_Pack(&(a->_mp_size), 1, MPI_INT, buf, *bufsize, pos, comm);
		MPI_Pack(&(a->_mp_prec), 1, MPI_INT, buf, *bufsize, pos, comm);
		MPI_Pack(&(a->_mp_exp), 1, MPI_LONG, buf, *bufsize, pos, comm);
		MPI_Pack(a->_mp_d, a->_mp_prec + 1, MPI_UNSIGNED_LONG, buf, *bufsize, pos, comm);
#endif
#endif
		//(unsigned char *)a = (void *)ptr_a;
		//a = (mpf_ptr)ptr_a;
		a++;
		//printf("%d\n", tmp_buf);
	}

	//(mpf_ptr)a = tmp_a;
	a = tmp_a;
}


//---------------------------------------------
double my_atod(char *in_str)
{
  char st[100],st1[100],st2[100];
  int d0,d1,d2,d3;
  double f0,f1,f2,f3;
  double out_f,out_f1;
  int i,l;
	  strcpy(st,in_str);
	  if(st[0]=='-'){
		  d0=-1;
		  strcpy(st,(in_str+1));
	  }
	  else{
		  d0=1;
	  }
	  l=strlen(st);
	  for(i=0;i<l;i++){
		  if(st[i]=='/'){
			  strcpy(st1,&st[i+1]);
			  st[i]='\x0';
			  out_f=my_atod(st);
			  out_f1=my_atod(st1);
			  f0=d0;
			  return f0*out_f/out_f1;
		  }
	  }
		  sscanf(&st[0],"%d.%s",&d1,&st2);
		  sscanf(&st[0],"%d.%d",&d1,&d2);
		  f0=d0;
		  f1=d1;
		  f2=d2;
		  d3=strlen(st2);
		  f3=1.0;
		  for(i=0;i<d3;i++){
			  f3=f3*10.0;
		  }
	out_f=f0*(f1+f2/f3);
	return out_f;
}

//---------------------------------------------
void mp_my_atod(char *in_str,mpfr_t mp_out,int prec)
{
  char st[100],st1[100],st2[100];
  int d0,d1,d2,d3;
  double f0,f1,f2,f3;
  double out_f,out_f1;
  mpfr_t mp_f1,mp_f2,mp_f3;
  mpfr_t mp_out_f,mp_out_f1;
  mpfr_t mp_out_f2,mp_out_f3,mp_f0;
  int i,l;
  mpfr_init2(mp_out_f, prec);
  mpfr_init2(mp_out_f1, prec);
  mpfr_init2(mp_out_f2, prec);
  mpfr_init2(mp_out_f3, prec);
  mpfr_init2(mp_f0, prec);
  mpfr_init2(mp_f1, prec);
  mpfr_init2(mp_f2, prec);
  mpfr_init2(mp_f3, prec);
	  strcpy(st,in_str);
	  if(st[0]=='-'){
		  d0=-1;
		  strcpy(st,(in_str+1));
	  }
	  else{
		  d0=1;
	  }
	  l=strlen(st);
	  for(i=0;i<l;i++){
		  if(st[i]=='/'){
			  strcpy(st1,&st[i+1]);
			  st[i]='\x0';
			  mp_my_atod(st,mp_out_f,prec);
			  mp_my_atod(st1,mp_out_f1,prec);
			  f0=d0;
			  mpfr_set_d(mp_f0,f0, GMP_RNDD);
			  mpfr_div(mp_out_f2,mp_out_f,mp_out_f1,GMP_RNDD);
			  mpfr_mul(mp_out_f3,mp_out_f2,mp_f0,GMP_RNDD);
			  mpfr_set(mp_out,mp_out_f3,GMP_RNDD);
			  return;
		  }
	  }
		  sscanf(&st[0],"%d.%s",&d1,&st2);
		  sscanf(&st[0],"%d.%d",&d1,&d2);
		  f0=d0;
		  f1=d1;
		  f2=d2;
		  d3=strlen(st2);
		    mpfr_set_d(mp_f3,1.0,GMP_RNDD);
		  for(i=0;i<d3;i++){
		    mpfr_mul_ui(mp_f3,mp_f3,10,GMP_RNDD);
		  }
	mpfr_set_d(mp_f1,f1,GMP_RNDD);
	mpfr_set_d(mp_f2,f2,GMP_RNDD);
	mpfr_set_d(mp_f0,f0,GMP_RNDD);

	mpfr_div(mp_f2,mp_f2,mp_f3,GMP_RNDD);
	mpfr_add(mp_f1,mp_f1,mp_f2,GMP_RNDD);
	mpfr_mul(mp_f0,mp_f0,mp_f1,GMP_RNDD);
	mpfr_set(mp_out,mp_f0, GMP_RNDD);
}
