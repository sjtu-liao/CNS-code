#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#include"mpi.h"
#include"gmp.h"
#include"mpfr.h"
#include"mpf2mpfr.h"
#define _GMP_H
#define _MPFR_H
#include"mpi_gmp.h"
#include"mpi_mpfr.h"

#define T_c 10000.0    /* critical predictable time */
#define gamma 1.1    /* safety factor */
#define delta_T 50.0    /* temporal interval of changing precision*/
#define out_T 1.0    /* temporal interval of output & delta_T/out_T must be an integer*/
#define remaining_T 500.0    /* precision is stopped decreasing after the time of T_c - remaining_T */

#define prec 20000
#define M 10000
mpfr_t X[M+1];
mpfr_t Y[M+1];
mpfr_t Z[M+1];

int myid, numprocs;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    mpfr_set_default_prec(prec);
    commit_mpf(&(MPI_MPF), prec, MPI_COMM_WORLD);
    create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);
    int MM, PRE;
    MM=M;
    double TT;
    double st, et;
    st=MPI_Wtime();
    FILE *fp;
    char filename[64];
    char Char[10000];
    double Tc, outT;
    mpfr_t Ym, YM, z, zz;
    mpfr_inits2(prec, Ym, YM, z, zz, (mpfr_ptr) 0);
    mpfr_t tol, CM, VSk1, VSk2, VSY1, VSY2, VSY, VSa, VSb;
    mpfr_inits2(prec, tol, CM, VSk1, VSk2, VSY1, VSY2, VSY, VSa, VSb, (mpfr_ptr) 0);
    mpfr_t temp, h, hi, t, T, PD, PONE, TEN;
    mpfr_inits2(prec, temp, h, hi, t, T, PD, PONE, TEN, (mpfr_ptr) 0);
    mpfr_set_str(TEN, "10.0", 10, GMP_RNDN);
    mpfr_t Mx, My, Mz, M3[3];
    mpfr_inits2(prec, Mx, My, Mz, M3[0], M3[1], M3[2], (mpfr_ptr) 0);
    mpfr_t Dxz, Dxy, D2[2];
    mpfr_inits2(prec, Dxz, Dxy, D2[0], D2[1], (mpfr_ptr) 0);
    mpfr_set_str(t, "0.0", 10, GMP_RNDN);
    Tc=T_c;
    mpfr_set_d(T, Tc, GMP_RNDN);
    mpfr_sprintf(Char, "%.10Re", T);
    mpfr_set_str(T, Char, 10, GMP_RNDN);
    outT=out_T;
    mpfr_set_d(PD, outT, GMP_RNDN);
    mpfr_sprintf(Char, "%.10Re", PD);
    mpfr_set_str(PD, Char, 10, GMP_RNDN);
    mpfr_set(PONE, PD, GMP_RNDN);
    int m, n;
    int i, j, k, Zint;
    int a, b1, b2, c, Compare, mm, PDD;
    Compare=0;
    PDD=1;
    a=10;
    b1=8;
    b2=3;
    c=28;
    void *packed_temp, *packed_TEMP;
    void *packed_temp2, *packed_temp3;
    packed_temp = allocbuf_mpf(prec, 1);
    packed_TEMP = allocbuf_mpf(prec, 1);
    packed_temp2 = allocbuf_mpf(prec, 2);
    packed_temp3 = allocbuf_mpf(prec, 3);
    for(m=0; m<=M; m++)
    {
        mpfr_inits2(prec, X[m], Y[m], Z[m], (mpfr_ptr) 0);
    }
    
    mpfr_set_str(X[0], "-15.8", 10, GMP_RNDN);    /* initial condition */
    mpfr_set_str(Y[0], "-17.48", 10, GMP_RNDN);
    mpfr_set_str(Z[0], "35.64", 10, GMP_RNDN);

m=0;
TT=m*out_T;
for(; Compare<=0; )
{
    if(myid==0)
    {
        if(PDD==1)
        {
            sprintf(filename, "XYZ_out-MP.dat");
            fp=fopen(filename, "a");
            fprintf(fp, "%lf\t", TT);
            mpfr_sprintf(Char, "%.5000Re", X[0]);
            fprintf(fp, "%.5010s\t", Char);
            mpfr_sprintf(Char, "%.5000Re", Y[0]);
            fprintf(fp, "%.5010s\t", Char);
            mpfr_sprintf(Char, "%.5000Re", Z[0]);
            fprintf(fp, "%.5010s\n", Char);
            fclose(fp);
        }
    }
    
    if(myid==0)
    {
        if(PDD==1)
        {
            sprintf(filename, "XYZ_out-double.dat");
            fp=fopen(filename, "a");
            fprintf(fp, "%lf\t", TT);
            mpfr_sprintf(Char, "%.16Re", X[0]);
            fprintf(fp, "%.26s\t", Char);
            mpfr_sprintf(Char, "%.16Re", Y[0]);
            fprintf(fp, "%.26s\t", Char);
            mpfr_sprintf(Char, "%.16Re", Z[0]);
            fprintf(fp, "%.26s\n", Char);
            fclose(fp);
        }
    }
    
    if( PDD==1 && m%(int)(delta_T/out_T)==0 && m<=(int)((T_c-remaining_T)/out_T) )
    {
        if(myid==0)
        {
            et=MPI_Wtime();
            printf("t=%lf: CPU's time is %fs\n", TT, et-st);
        }
        
        PRE=ceil(gamma*(T_c-TT)/2.53);    /* significant digit number of multiple-precision */
        
        MM=ceil(1.5*PRE);    /* order of Taylor expansion */
        
        Zint=(-1.0)*PRE;
        mpfr_pow_si(tol, TEN, Zint, GMP_RNDN);    /* allowed tolerance for the optimal time step */
        
        mpfr_ui_pow_ui(temp, 10, PRE, GMP_RNDN);
        mpfr_log2(temp, temp, GMP_RNDN);
        PRE=mpfr_get_ui(temp, GMP_RNDN)+1;
        
        mpfr_set_default_prec(PRE);
        commit_mpf(&(MPI_MPF), PRE, MPI_COMM_WORLD);
        create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);
        
        for(i=0; i<=MM; i++)
        {
            mpfr_prec_round(X[i], PRE, GMP_RNDN);
            mpfr_prec_round(Y[i], PRE, GMP_RNDN);
            mpfr_prec_round(Z[i], PRE, GMP_RNDN);
        }
        
        mpfr_prec_round(temp, PRE, GMP_RNDN);
        mpfr_prec_round(h, PRE, GMP_RNDN);
        mpfr_prec_round(hi, PRE, GMP_RNDN);
        mpfr_prec_round(t, PRE, GMP_RNDN);
        mpfr_prec_round(T, PRE, GMP_RNDN);
        mpfr_prec_round(PD, PRE, GMP_RNDN);
        mpfr_prec_round(PONE, PRE, GMP_RNDN);
        mpfr_prec_round(TEN, PRE, GMP_RNDN);
        mpfr_prec_round(Mx, PRE, GMP_RNDN);
        mpfr_prec_round(My, PRE, GMP_RNDN);
        mpfr_prec_round(Mz, PRE, GMP_RNDN);
        mpfr_prec_round(Dxz, PRE, GMP_RNDN);
        mpfr_prec_round(Dxy, PRE, GMP_RNDN);
        mpfr_prec_round(D2[0], PRE, GMP_RNDN);
        mpfr_prec_round(D2[1], PRE, GMP_RNDN);
        mpfr_prec_round(M3[0], PRE, GMP_RNDN);
        mpfr_prec_round(M3[1], PRE, GMP_RNDN);
        mpfr_prec_round(M3[2], PRE, GMP_RNDN);
        
        mpfr_prec_round(Ym, PRE, GMP_RNDN);
        mpfr_prec_round(YM, PRE, GMP_RNDN);
        mpfr_prec_round(z, PRE, GMP_RNDN);
        mpfr_prec_round(zz, PRE, GMP_RNDN);
        mpfr_prec_round(tol, PRE, GMP_RNDN);
        mpfr_prec_round(CM, PRE, GMP_RNDN);
        mpfr_prec_round(VSk1, PRE, GMP_RNDN);
        mpfr_prec_round(VSk2, PRE, GMP_RNDN);
        mpfr_prec_round(VSY1, PRE, GMP_RNDN);
        mpfr_prec_round(VSY2, PRE, GMP_RNDN);
        mpfr_prec_round(VSY, PRE, GMP_RNDN);
        mpfr_prec_round(VSa, PRE, GMP_RNDN);
        mpfr_prec_round(VSb, PRE, GMP_RNDN);
        
        packed_temp = allocbuf_mpf(PRE, 1);
        packed_TEMP = allocbuf_mpf(PRE, 1);
        packed_temp2 = allocbuf_mpf(PRE, 2);
        packed_temp3 = allocbuf_mpf(PRE, 3);
    }
    
        mm=MM;
        mpfr_set_str(CM, "1.0", 10, GMP_RNDN);
        mpfr_div_si(CM, CM, mm, GMP_RNDN);
        mpfr_pow(VSk1, tol, CM, GMP_RNDN);
        mm=MM-1;
        mpfr_set_str(CM, "-1.0", 10, GMP_RNDN);
        mpfr_div_si(VSY1, CM, mm, GMP_RNDN);
        mm=MM+1;
        mpfr_set_str(CM, "1.0", 10, GMP_RNDN);
        mpfr_div_si(CM, CM, mm, GMP_RNDN);
        mpfr_pow(VSk2, tol, CM, GMP_RNDN);
        mm=MM;
        mpfr_set_str(CM, "-1.0", 10, GMP_RNDN);
        mpfr_div_si(VSY2, CM, mm, GMP_RNDN);
    
    for(i=1; i<=MM; i++)
    {
        mpfr_set_str(Dxz, "0.0", 10, GMP_RNDN);
        mpfr_set_str(Dxy, "0.0", 10, GMP_RNDN);
        for(j=myid; j<=i-1; j+=numprocs)
        {
            mpfr_mul(temp, Z[i-1-j], X[j], GMP_RNDN);
            mpfr_add(Dxz, Dxz, temp, GMP_RNDN);
            mpfr_mul(temp, Y[i-1-j], X[j], GMP_RNDN);
            mpfr_add(Dxy, Dxy, temp, GMP_RNDN);
        }
        pack_mpf(Dxz, 1, packed_temp);
        MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
        /*MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);*/
        if(myid==0)
            unpack_mpf(packed_TEMP, D2[0], 1);
        
        pack_mpf(Dxy, 1, packed_temp);
        MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
        /*MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);*/
        if(myid==0)
            unpack_mpf(packed_TEMP, D2[1], 1);
        
        if(myid==0)
        {
            pack_mpf(D2[0], 2, packed_temp2);
        }
        MPI_Bcast(packed_temp2, 2, MPI_MPF, 0, MPI_COMM_WORLD);
        unpack_mpf(packed_temp2, D2[0], 2);
        
        mpfr_sub(X[i], Y[i-1], X[i-1], GMP_RNDN);
        mpfr_mul_ui(X[i], X[i], a, GMP_RNDN);
        mpfr_div_ui(X[i], X[i], i, GMP_RNDN);
        mpfr_mul_ui(Y[i], X[i-1], c, GMP_RNDN);
        mpfr_sub(Y[i], Y[i], Y[i-1], GMP_RNDN);
        mpfr_sub(Y[i], Y[i], D2[0], GMP_RNDN);
        mpfr_div_ui(Y[i], Y[i], i, GMP_RNDN);
        mpfr_mul_ui(Z[i], Z[i-1], b1, GMP_RNDN);
        mpfr_div_ui(Z[i], Z[i], b2, GMP_RNDN);
        mpfr_sub(Z[i], D2[1], Z[i], GMP_RNDN);
        mpfr_div_ui(Z[i], Z[i], i, GMP_RNDN);
    }
    
            mpfr_abs(Ym, X[MM-1], GMP_RNDN);
            mpfr_abs(YM, X[MM], GMP_RNDN);
            mpfr_abs(z, Y[MM-1], GMP_RNDN);
            mpfr_abs(zz, Y[MM], GMP_RNDN);
            Compare=mpfr_cmp(z, Ym);
            if(Compare>0)
            {
                mpfr_set(Ym, z, GMP_RNDN);
            }
            Compare=mpfr_cmp(zz, YM);
            if(Compare>0)
            {
                mpfr_set(YM, zz, GMP_RNDN);
            }
            mpfr_abs(z, Z[MM-1], GMP_RNDN);
            mpfr_abs(zz, Z[MM], GMP_RNDN);
            Compare=mpfr_cmp(z, Ym);
            if(Compare>0)
            {
                mpfr_set(Ym, z, GMP_RNDN);
            }
            Compare=mpfr_cmp(zz, YM);
            if(Compare>0)
            {
                mpfr_set(YM, zz, GMP_RNDN);
            }
        mpfr_pow(VSY, Ym, VSY1, GMP_RNDN);
        mpfr_mul(VSa, VSk1, VSY, GMP_RNDN);
        mpfr_pow(VSY, YM, VSY2, GMP_RNDN);
        mpfr_mul(VSb, VSk2, VSY, GMP_RNDN);
        Compare=mpfr_cmp(VSa, VSb);
        if(Compare>=0)
        {
            mpfr_set(h, VSb, GMP_RNDN);
        }
        if(Compare<0)
        {
            mpfr_set(h, VSa, GMP_RNDN);
        }
        
        PDD=0;
        mpfr_add(z, t, h, GMP_RNDN);
        Compare=mpfr_cmp(z, PD);
        if(Compare>=0)
        {
            mpfr_sub(h, PD, t, GMP_RNDN);
            mpfr_add(PD, PD, PONE, GMP_RNDN);
            PDD++;
            m++;
            TT=m*out_T;
        }
    
    mpfr_set_str(Mx, "0.0", 10, GMP_RNDN);
    mpfr_set_str(My, "0.0", 10, GMP_RNDN);
    mpfr_set_str(Mz, "0.0", 10, GMP_RNDN);
    for(i=(myid+1); i<=MM; i+=numprocs)
    {
        mpfr_set_str(hi, "1.0", 10, GMP_RNDN);
        for(k=1; k<=i; k++)
        {
            mpfr_mul(hi, hi, h, GMP_RNDN);
        }
        mpfr_mul(temp, hi, X[i], GMP_RNDN);
        mpfr_add(Mx, Mx, temp, GMP_RNDN);
        mpfr_mul(temp, hi, Y[i], GMP_RNDN);
        mpfr_add(My, My, temp, GMP_RNDN);
        mpfr_mul(temp, hi, Z[i], GMP_RNDN);
        mpfr_add(Mz, Mz, temp, GMP_RNDN);
    }
    
    pack_mpf(Mx, 1, packed_temp);
    MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
    /*MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);*/
    if(myid==0)
        unpack_mpf(packed_TEMP, M3[0], 1);
    
    pack_mpf(My, 1, packed_temp);
    MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
    /*MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);*/
    if(myid==0)
        unpack_mpf(packed_TEMP, M3[1], 1);
    
    pack_mpf(Mz, 1, packed_temp);
    MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
    /*MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);*/
    if(myid==0)
        unpack_mpf(packed_TEMP, M3[2], 1);
    
    if(myid==0)
    {
        pack_mpf(M3[0], 3, packed_temp3);
    }
    MPI_Bcast(packed_temp3, 3, MPI_MPF, 0, MPI_COMM_WORLD);
    unpack_mpf(packed_temp3, M3[0], 3);
    
    mpfr_add(X[0], X[0], M3[0], GMP_RNDN);
    mpfr_add(Y[0], Y[0], M3[1], GMP_RNDN);
    mpfr_add(Z[0], Z[0], M3[2], GMP_RNDN);
    
    mpfr_add(t, t, h, GMP_RNDN);
    Compare=mpfr_cmp(t, T);
}
    if(myid==0)
    {
        et=MPI_Wtime();
        printf("Total CPU's time is %fs\n", et-st);
    }
    free_mpf_op(&(MPI_MPF_SUM));
    free_mpf(&(MPI_MPF));
    MPI_Finalize();
    return 0;
}
