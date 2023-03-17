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

#define prec 13886    /* number of significant digits */
#define M 3500    /* order of Taylor expansion */

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
    double st, et;
    st=MPI_Wtime();
    FILE *fp;
    char filename[64];
    char Char[5000];
    mpfr_t temp, h, hi;
    mpfr_inits2(prec, temp, h, hi, (mpfr_ptr) 0);
    mpfr_t Mx, My, Mz, MX, MY, MZ;
    mpfr_inits2(prec, Mx, My, Mz, MX, MY, MZ, (mpfr_ptr) 0);
    mpfr_t Dxz, Dxy, DXZ, DXY;
    mpfr_inits2(prec, Dxz, Dxy, DXZ, DXY, (mpfr_ptr) 0);
    int H, m, n;
    double T, TT, out_T, h_double;
    mpfr_set_str(h, "0.01", 10, GMP_RNDN);    /* time-step */
    T=10000.0;    /* computation time (Lorenz unit time) */
    out_T=1.0;    /* temporal interval of output */
    h_double=mpfr_get_d(h, GMP_RNDN);
    H=T/h_double;
    int i, j, k, Zint;
    int a, b1, b2, c;
    a=10;    /* sigma */
    b1=8; b2=3;    /* b=-b1/b2 */
    c=28;    /* R */
    void *packed_temp, *packed_TEMP;
    packed_temp = allocbuf_mpf(prec, 1);
    packed_TEMP = allocbuf_mpf(prec, 1);
    for(m=0; m<=M; m++)
    {
        mpfr_inits2(prec, X[m], Y[m], Z[m], (mpfr_ptr) 0);
    }
    mpfr_set_str(X[0], "-15.8", 10, GMP_RNDN);    /* initial condition */
    mpfr_set_str(Y[0], "-17.48", 10, GMP_RNDN);
    mpfr_set_str(Z[0], "35.64", 10, GMP_RNDN);
    
for(m=0; m<=H; m++)
{
    TT=m*h_double;
    if(myid==0)
    {
        if(m%(int)(out_T/h_double)==0)
        {
            sprintf(filename, "XYZ_out-MP.dat");
            fp=fopen(filename, "a");
            fprintf(fp, "%lf\t", TT);
            mpfr_sprintf(Char, "%.4180Re", X[0]);
            fprintf(fp, "%.4190s\t", Char);
            mpfr_sprintf(Char, "%.4180Re", Y[0]);
            fprintf(fp, "%.4190s\t", Char);
            mpfr_sprintf(Char, "%.4180Re", Z[0]);
            fprintf(fp, "%.4190s\n", Char);
            fclose(fp);
        }
    }
    
    if(myid==1)
    {
        if(m%(int)(out_T/h_double)==0)
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
    
    for(i=1; i<=M; i++)
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
        MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);
        unpack_mpf(packed_TEMP, DXZ, 1);
        
        pack_mpf(Dxy, 1, packed_temp);
        MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);
        unpack_mpf(packed_TEMP, DXY, 1);
        
        mpfr_sub(X[i], Y[i-1], X[i-1], GMP_RNDN);
        mpfr_mul_ui(X[i], X[i], a, GMP_RNDN);
        mpfr_div_ui(X[i], X[i], i, GMP_RNDN);
        mpfr_mul_ui(Y[i], X[i-1], c, GMP_RNDN);
        mpfr_sub(Y[i], Y[i], Y[i-1], GMP_RNDN);
        mpfr_sub(Y[i], Y[i], DXZ, GMP_RNDN);
        mpfr_div_ui(Y[i], Y[i], i, GMP_RNDN);
        mpfr_mul_ui(Z[i], Z[i-1], b1, GMP_RNDN);
        mpfr_div_ui(Z[i], Z[i], b2, GMP_RNDN);
        mpfr_sub(Z[i], DXY, Z[i], GMP_RNDN);
        mpfr_div_ui(Z[i], Z[i], i, GMP_RNDN);
    }
    mpfr_set_str(Mx, "0.0", 10, GMP_RNDN);
    mpfr_set_str(My, "0.0", 10, GMP_RNDN);
    mpfr_set_str(Mz, "0.0", 10, GMP_RNDN);
    for(i=(myid+1); i<=M; i+=numprocs)
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
    MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);
    unpack_mpf(packed_TEMP, MX, 1);
    
    pack_mpf(My, 1, packed_temp);
    MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);
    unpack_mpf(packed_TEMP, MY, 1);
    
    pack_mpf(Mz, 1, packed_temp);
    MPI_Reduce(packed_temp, packed_TEMP, 1, MPI_MPF, MPI_MPF_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(packed_TEMP, 1, MPI_MPF, 0, MPI_COMM_WORLD);
    unpack_mpf(packed_TEMP, MZ, 1);
    
    mpfr_add(X[0], X[0], MX, GMP_RNDN);
    mpfr_add(Y[0], Y[0], MY, GMP_RNDN);
    mpfr_add(Z[0], Z[0], MZ, GMP_RNDN);
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
