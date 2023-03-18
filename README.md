# CNS code

To gain a reproducible/reliable numerical simulation of chaotic system, Liao [1] proposed a brand-new numerical strategy, namely the "Clean Numerical Simulation" (CNS) [2-4], to control the background numerical noise, say, truncation and round-off errors, during a temporal interval [0, T<sub>c</sub>], where T<sub>c</sub> is the so-called "critical predictable time" and this temporal interval should be long enough for calculating statistics. In the frame of the CNS [1-10], the temporal and spatial truncation errors can be decreased to a required small level via using the Taylor expansion method with a high enough order in the temporal dimension and adopting a fine enough discretization method in the spatial dimension (such as the high-order spatial Fourier expansion), respectively. Significantly, all of the physical and numerical variables/parameters should be represented by means of the multiple precision with a large enough number of significant digits, and thus the round-off error can also be decreased to a required small level. In this way, different from some other general numerical algorithms, the CNS is able to give the reproducible/reliable numerical simulation of a chaotic dynamical system within a finite but long enough temporal interval.

Up to now, the CNS has been applied to many chaotic dynamical systems successfully with the corresponding computer-generated simulations being reproducible and of course reliable.

## Lorenz system

In 2014, Liao & Wang [4] obtained a reproducible/convergent numerical simulation of the Lorenz system in quite a long temporal interval [0, 10000] (Lorenz unit time) by means of a parallel algorithm of the CNS using the 3500th-order Taylor expansion with the constant time step ∆t = 0.01 and 4180-digit multiple precision for all physical and numerical variables/parameters. It took 220.9 hours (i.e. about 9 days and 5 hours) using 1200 CPUs of the National Supercomputer TH-1A at Tianjin, China.

For the code, please refer to the website:

[Code-Lorenz](https://github.com/sjtu-liao/CNS-code/blob/master/Code-Lorenz)

**"Lorenz.c" is the main program written in C language using the MPFR library and MPI parallel technique, which needs two header files "mpi_gmp.h" and "mpi_mpfr.h": the number of significant digits and the order of Taylor expansion can be set at Lines 15 & 16, and the settings of other numerical parameters are annotated. One can compile the main program via the command *mpicc* such as that written in "makefile", and run the compiled executable file via the command *mpirun* such as that written in "run.sh".**

In 2023, according to the exponentially growing property of noise in numerical simulations of chaos, Qin and Liao [11] proposed a modified strategy of the CNS, called the "self-adaptive CNS", to significantly increase the computational efficiency of the CNS algorithm. Using this self-adaptive CNS strategy, a reproducible/convergent numerical simulation of the above-mentioned Lorenz system in the temporal interval [0, 10000] with the critical predictable time T<sub>c</sub>=10000 is obtained successfully, by means of the parallel CNS algorithm using 1200 CPUs of the National Supercomputer TH-2 at Guangzhou, China. Note that it took only 37.2 hours (i.e. about 1 day and 13 hours), just about 17% of the CPU time of the previous CNS algorithm applied by Liao & Wang [4] who likewise used a supercomputer.

For the code, please refer to the website:

[Code-Lorenz-SA](https://github.com/sjtu-liao/CNS-code/blob/master/Code-Lorenz-SA)

**"Lorenz.c" is the main program written in C language using the MPFR library and MPI parallel technique, which needs two header files "mpi_gmp.h" and "mpi_mpfr.h": critical predictable time, safety factor, temporal interval of changing precision, and so on can be set at Lines 15-19. One can compile the main program via the command *mpicc* such as that written in "makefile", and run the compiled executable file via the command *mpirun* such as that written in "run.sh".**

---

### References
[1] S. Liao, On the reliability of computed chaotic solutions of non-linear differential equations, Tellus Ser. A-Dyn. Meteorol. Oceanol. 61 (4) (2009) 550–564.

[2] S. Liao, On the numerical simulation of propagation of micro-level inherent uncertainty for chaotic dynamic systems, Chaos Solitons Fractals 47 (2013) 1–12.

[3] S. Liao, Physical limit of prediction for chaotic motion of three-body problem, Commun. Nonlinear Sci. Numer. Simul. 19 (3) (2014) 601–616.

[4] S. Liao, P. Wang, On the mathematically reliable long-term simulation of chaotic solutions of lorenz equation in the interval [0, 10000], Sci. China-Phys. Mech. Astron. 57 (2) (2014) 330–335.

[5] T. Hu, S. Liao, On the risks of using double precision in numerical simulations of spatio-temporal chaos, J. Comput. Phys. 418 (2020) 109629.

[6] S. Qin, S. Liao, Influence of numerical noises on computer-generated simulation of spatio-temporal chaos, Chaos Solitons Fractals 136 (2020) 109790.

[7] T. Xu, J. Li, Z. Li, S. Liao, Accurate predictions of chaotic motion of a free fall disk, Phys. Fluids 33 (3) (2021) 037111.

[8] S. Liao, S. Qin, Ultra-chaos: an insurmountable objective obstacle of reproducibility and replication, Adv. Appl. Math. Mech. 14 (4) (2022) 799–815.

[9] S. Qin, S. Liao, Large-scale influence of numerical noises as artificial stochastic disturbances on a sustained turbulence, J. Fluid Mech. 948 (2022) A7.

[10] Y. Yang, S. Qin, S. Liao, Ultra-chaos of a mobile robot: a higher disorder than normal-chaos, Chaos Solitons Fractals 167 (2023) 113037.

[11] S. Qin, S. Liao, A self-adaptive algorithm of the clean numerical simulation (CNS) for chaos, Adv. Appl. Math. Mech. (Accepted).
