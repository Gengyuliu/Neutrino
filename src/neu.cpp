#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <iostream>
#include <fstream>
#include "neu.hpp"

#define sgn(v) ((v>0) - (v<0)) //sign function in branchless implementation


void Neu::parsing(int argc, char *const *argv, phys_sys_t *phys_sys, num_sys_t *num_sys){
    int opt;
    const char *optstr = "";
    const struct option long_opts[] = {
        {"flavor", required_argument, NULL, 1},
        {"theta", required_argument, NULL, 2},
        {"mu", required_argument, NULL, 3},
        {"f_0", required_argument, NULL, 4},
        {"alpha", required_argument, NULL, 5},
        {"beta", required_argument, NULL, 6},
        {"dim", required_argument, NULL, 7},
        {"Nt", required_argument, NULL, 8},
        {"CFL", required_argument, NULL, 9},
        {"Nx", required_argument, NULL, 10},
        {"x1", required_argument, NULL, 11},
        {"x2", required_argument, NULL, 12},
        {"Nvx", required_argument, NULL, 13},
        {"Ny", required_argument, NULL, 14},
        {"y1", required_argument, NULL, 15},
        {"y2", required_argument, NULL, 16},
        {"Nvy", required_argument, NULL, 17},
        {"Nz", required_argument, NULL, 18},
        {"z1", required_argument, NULL, 19},
        {"z2", required_argument, NULL, 20},
        {"Nvz", required_argument, NULL, 21},
        {"init", required_argument, NULL, 22}, 
        {"A", required_argument, NULL, 23},
        {"z_0", required_argument, NULL, 24},
        {"sigma", required_argument, NULL, 25},
        {0,0,0,0}
    };
    double val;
    while ((opt = getopt_long(argc, argv, optstr, long_opts, NULL)) != -1){
        switch(opt){
            case 1:
                sscanf(optarg, "%d", &phys_sys->flavor);
                break;
            case 2:
                sscanf(optarg, "%lf", &val);
                phys_sys->theta = val*M_PI/180.0;
                break;
            case 3:
                sscanf(optarg, "%lf", &phys_sys->mu);
                break;
            case 4:
                sscanf(optarg, "%lf", &phys_sys->f_0);
                break;
            case 5:
                sscanf(optarg, "%lf", &phys_sys->alpha);
                break;
            case 6:
                sscanf(optarg, "%lf", &phys_sys->beta);
                break;
            case 7:
                sscanf(optarg, "%d", &num_sys->dim);
                break;
            case 8:
                sscanf(optarg, "%d", &num_sys->Nt);
                break;
            case 9:
                sscanf(optarg, "%lf", &num_sys->CFL);
                break;
            case 10:
                sscanf(optarg, "%d", &num_sys->Nx);
                break;
            case 11:
                sscanf(optarg, "%lf", &num_sys->x1);
                break; 
            case 12:
                sscanf(optarg, "%lf", &num_sys->x2);
                break;
            case 13:
                sscanf(optarg, "%d", &num_sys->Nvx);
                break;
            case 14:
                sscanf(optarg, "%d", &num_sys->Ny);
                break;
            case 15:
                sscanf(optarg, "%lf", &num_sys->y1);
                break; 
            case 16:
                sscanf(optarg, "%lf", &num_sys->y2);
                break;
            case 17:
                sscanf(optarg, "%d", &num_sys->Nvy);
                break;
            case 18:
                sscanf(optarg, "%d", &num_sys->Nz);
                break;
            case 19:
                sscanf(optarg, "%lf", &num_sys->z1);
                break; 
            case 20:
                sscanf(optarg, "%lf", &num_sys->z2);
                break;
            case 21:
                sscanf(optarg, "%d", &num_sys->Nvz);
                break;
            case 22:
                sscanf(optarg, "%d", &phys_sys->init);
                break;
            case 23:
                sscanf(optarg, "%lf", &phys_sys->A);
                //osc->A = 1e-7;
                //printf("this is A : %.8lf\n", phys_sys->A);
                break;
            case 24:
                sscanf(optarg, "%lf", &phys_sys->z_0);
                break;
            case 25:
                sscanf(optarg, "%lf", &phys_sys->sigma);
                break;
        }
    }
}

#define MAT2_ADD(C, a, A, b, B) {   \
    C._11 = a*A._11 + b*B._11;      \
    C._12 = a*A._12 + b*B._12;      \
    C._21 = a*A._21 + b*B._21;      \
    C._22 = a*A._22 + b*B._22;      \
}

//#define ADV2(y, a, x, N) {      \
//    y._11 += a*(8.0*(x[1*N]._11 - x[-1*N]._11) - (x[2*N]._11 - x[-2*N]._11)); \
//    y._12 += a*(8.0*(x[1*N]._12 - x[-1*N]._12) - (x[2*N]._12 - x[-2*N]._12)); \
//    y._21 += a*(8.0*(x[1*N]._21 - x[-1*N]._21) - (x[2*N]._21 - x[-2*N]._21)); \
//    y._22 += a*(8.0*(x[1*N]._22 - x[-1*N]._22) - (x[2*N]._22 - x[-2*N]._22)); \
//}

static inline void ADV2(op2_t *y, double a, op2_t *x, size_t N){
    y->_11 += a*(8.0*(x[1*N]._11 - x[-1*N]._11) - (x[2*N]._11 - x[-2*N]._11));
    y->_12 += a*(8.0*(x[1*N]._12 - x[-1*N]._12) - (x[2*N]._12 - x[-2*N]._12));
    y->_21 += a*(8.0*(x[1*N]._21 - x[-1*N]._21) - (x[2*N]._21 - x[-2*N]._21));
    y->_22 += a*(8.0*(x[1*N]._22 - x[-1*N]._22) - (x[2*N]._22 - x[-2*N]._22));

}

static inline void KO2(op2_t *y, double a, op2_t *x, size_t N){
    y->_11 += a*(x[2*N]._11+x[-2*N]._11-4*(x[1*N]._11+x[-1*N]._11)+6*x[0]._11);
    y->_12 += a*(x[2*N]._12+x[-2*N]._12-4*(x[1*N]._12+x[-1*N]._12)+6*x[0]._12); 
    y->_21 += a*(x[2*N]._21+x[-2*N]._21-4*(x[1*N]._21+x[-1*N]._21)+6*x[0]._21); 
    y->_22 += a*(x[2*N]._22+x[-2*N]._22-4*(x[1*N]._22+x[-1*N]._22)+6*x[0]._22); 
}

#define VAC2(H, theta){ \
    H.v._11 = -cos(2.0*theta);  \
    H.v._12 = sin(2.0*theta);   \
    H.v._21 = sin(2.0*theta);   \
    H.v._22 = cos(2.0*theta);   \
    H.v_bar = H.v;              \
}

#define CONJ2(A, B){ \
    A._11 = B._11;          \
    A._12 = conj(B._12);    \
    A._21 = conj(B._21);    \
    A._22 = B._22;          \
}

#define INT2(A, a, B, C){  \
    A._11 += a*(B._11 - C._11); \
    A._12 += a*(B._12 - C._12); \
    A._21 += a*(B._21 - C._21); \
    A._22 += a*(B._22 - C._22); \
}

#define COMMUTATOR2(C, a, A, B){ \
    C._11 = a*(A._12*B._21 - A._21*B._12);  \
    C._12 = a*((A._11 - A._22)*B._12 + A._12*(B._22 - B._11));  \
    C._21 = a*((A._22 - A._11)*B._21 + A._21*(B._11 - B._22));  \
    C._22 = a*(A._21*B._12 - A._12*B._21);  \
}

//#define KO2(y, a, x, N) { \
//    y._11 += a*(x[2*N]._11+x[-2*N]._11 - 4*(x[1*N]._11+x[-1*N]._11) + 6*x[0]._11); \
//    y._12 += a*(x[2*N]._12+x[-2*N]._12 - 4*(x[1*N]._12+x[-1*N]._12) + 6*x[0]._12); \
//    y._21 += a*(x[2*N]._21+x[-2*N]._21 - 4*(x[1*N]._21+x[-1*N]._21) + 6*x[0]._21); \
//    y._22 += a*(x[2*N]._22+x[-2*N]._22 - 4*(x[1*N]._22+x[-1*N]._22) + 6*x[0]._22); \
//}
//
void Neu::Qke2::rk4(field2_t *curr, field2_t *next) {
    field2_t f0, f1, f2, f3;
    field2_t curr1, curr2, curr3; 
    for (size_t i = 0; i < SIZE; ++i){
        bzero(&next->v[i], sizeof(op2_t));
        bzero(&next->v_bar[i], sizeof(op2_t));
    }
    f0.v = new op2_t[SIZE]; f0.v_bar = new op2_t[SIZE];
    f1.v = new op2_t[SIZE]; f1.v_bar = new op2_t[SIZE];
    f2.v = new op2_t[SIZE]; f2.v_bar = new op2_t[SIZE];
    f3.v = new op2_t[SIZE]; f3.v_bar = new op2_t[SIZE];
    curr1.v = new op2_t[SIZE]; curr1.v_bar = new op2_t[SIZE];
    curr2.v = new op2_t[SIZE]; curr2.v_bar = new op2_t[SIZE];
    curr3.v = new op2_t[SIZE]; curr3.v_bar = new op2_t[SIZE];

    /*	Step 1	*/    
    periodic_ghost_zone(curr);
    compute_F(&f0, curr);
    fieldadd(&curr1, 1.0, curr, 0.5*dt, &f0);

    /*	Step 2	*/
    periodic_ghost_zone(&curr1);
    compute_F(&f1, &curr1);
    fieldadd(&curr2, 1.0, curr, 0.5*dt, &f1);
    
    /*	Step 3	*/
    periodic_ghost_zone(&curr2);
    compute_F(&f2, &curr2);
    fieldadd(&curr3, 1.0, curr, dt, &f2);

    /*	Step 4	*/
    periodic_ghost_zone(&curr3);
    compute_F(&f3, &curr3);

    fieldadd(next, 1.0, next, 1.0, curr);
    fieldadd(next, 1.0, next, dt/6.0, &f0);
    fieldadd(next, 1.0, next, dt/3.0, &f1);
    fieldadd(next, 1.0, next, dt/3.0, &f2);
    fieldadd(next, 1.0, next, dt/6.0, &f3);

    
    delete[] f0.v; delete[] f0.v_bar;
    delete[] f1.v; delete[] f1.v_bar;
    delete[] f2.v; delete[] f2.v_bar;
    delete[] f3.v; delete[] f3.v_bar;
    delete[] curr1.v; delete[] curr1.v_bar;
    delete[] curr2.v; delete[] curr2.v_bar;
    delete[] curr3.v; delete[] curr3.v_bar;

}


void Neu::Qke2::fieldadd(field2_t *C, const cnum a, const field2_t *A, const cnum b, const field2_t *B) {
    switch (sys.dim){
        case 1:
#pragma omp parallel for
           for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    //int idx = _index(i2, k, osc);
                    //int _j = j%sys.Nz;
                    //if (_j < 0)
                    //    _j += sys.Nz;
                    size_t idx = index(sys.dim, j, k);

                    MAT2_ADD(C->v[idx], a, A->v[idx], b, B->v[idx]);
                    MAT2_ADD(C->v_bar[idx], a, A->v_bar[idx], b, B->v_bar[idx]);
                }
            }
            break;
        case 2:
#pragma omp parallel for
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue; 
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            size_t idx = index(sys.dim, jx, jz, kx, kz);
                            MAT2_ADD(C->v[idx], a, A->v[idx], b, B->v[idx]);
                            MAT2_ADD(C->v_bar[idx], a, A->v_bar[idx], b, B->v_bar[idx]);
                        }
                    }
                }
            
            } 
            break;
        case 3:
            break;
    }
}

void Neu::Qke2::periodic_ghost_zone(field2_t *A) {
    switch (sys.dim){
        case 1:
#pragma omp parallel for
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < gz; ++j){
                    /*	lower side	*/ 
                    size_t idx = index(sys.dim, -j-1, k);
                    A->v[idx] = A->v[index(sys.dim, sys.Nz-j-1,k)];
                    A->v_bar[idx] = A->v_bar[index(sys.dim, sys.Nz-j-1,k)];
                    /*	lower side	*/
                    idx = index(sys.dim, sys.Nz+j,k);
                    A->v[idx] = A->v[index(sys.dim, j, k)];
                    A->v_bar[idx] = A->v_bar[index(sys.dim, j, k)];
                }
            }
            break;
        case 2:
#pragma omp parallel for 
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2), + pow(vx[kx],2) > 1)
                        continue;
                    for (int jz = 0; jz < gz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            /*	z	*/ 
                            size_t idx = index(sys.dim, jx, -jz-1, kx, kz);
                            size_t id1 = index(sys.dim, jx,sys.Nz-jz-1,kx,kz); 
                            size_t id2 = index(sys.dim, jx,jz, kx, kz);
                            A->v[idx] = A->v[id1];
                            A->v_bar[idx] = A->v_bar[id1];
                            
                            idx = index(sys.dim, jx, sys.Nz+jz, kx, kz);
                            A->v[idx] = A->v[id2];
                            A->v_bar[idx] = A->v_bar[id2];
                        }
                    }
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < gx; ++jx){
                            /*	x	*/
                            size_t idx = index(sys.dim, -jx-1, jz, kx, kz);
                            size_t id1 = index(sys.dim, sys.Nx-jx-1,jz,kx,kz);
                            size_t id2 = index(sys.dim, jx, jz, kx, kz);
                            A->v[idx] = A->v[id1];
                            A->v_bar[idx] = A->v_bar[id1];
                            
                            idx = index(sys.dim, sys.Nx+jx, jz, kx, kz);
                            A->v[idx] = A->v[id2];
                            A->v_bar[idx] = A->v_bar[id2];

                        }
                    }
                }
            }
            break;
        case 3:
            break;
    } 

}
void Neu::Qke2::compute_F(field2_t *f, const field2_t *a) {
    switch (sys.dim){
        case 1:
#pragma omp parallel for
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    size_t idx = index(sys.dim, j, k);
                    double v = vz[k];
                    op2_t *a_ptr    = &a->v[idx];  
                    op2_t *a_bptr   = &a->v_bar[idx];
                    pair2_t f_adv, f_osc;
                    pair2_t H, H_vac, H_vv;
                    
                    memset(&f_adv.v, 0, sizeof(op2_t));
                    memset(&f_adv.v_bar, 0, sizeof(op2_t));
                    memset(&f_osc.v, 0, sizeof(op2_t));
                    memset(&f_osc.v_bar, 0, sizeof(op2_t));
                    memset(&H.v, 0, sizeof(op2_t));
                    memset(&H.v_bar, 0, sizeof(op2_t));
#ifdef ADV 
                    /*	advection	*/
                    double factor = -v/(12.0*dz);
                    ADV2(&f_adv.v, factor, a_ptr, 1);
                    ADV2(&f_adv.v_bar, factor, a_bptr, 1);
#endif

#ifdef KO
/*	Kreiss-Oliger dissipation (3-rd order)	*/
                    double ko_eps = -1e-4/(dz*16);
                    KO2(&f_adv.v, ko_eps, a_ptr, 1);
                    KO2(&f_adv.v_bar, ko_eps, a_bptr, 1);
#endif
                    /*	oscillation	*/
#ifdef VACUUM
                    VAC2(H, osc.theta);
#elif INT_VV
                    for (int kk = 0; kk < sys.Nvz; ++kk){
                        size_t idxx = index(sys.dim, j, kk);
                        double vv = vz[kk];
                        op2_t a_conj, a_bar_conj;
                        CONJ2(a_bar_conj, a->v_bar[idxx]);
                        INT2(H.v, osc.mu*dvz*(1.0 - v*vv), a->v[idxx], a_bar_conj);

                        CONJ2(a_conj, a->v[idxx]);
                        INT2(H.v_bar, -osc.mu*dvz*(1.0-v*vv), a_conj, a->v_bar[idxx]);
                    } 

#elif VAC_VV
                    VAC2(H_vac, osc.theta);
                    memset(&H_vv.v, 0, sizeof(op2_t));
                    memset(&H_vv.v_bar, 0, sizeof(op2_t));
                    for (int kk = 0; kk < sys.Nvz; ++kk){
                        size_t idxx = index(sys.dim, j, kk);
                        double vv = vz[kk];
                        op2_t a_conj, a_bar_conj;
                        CONJ2(a_bar_conj, a->v_bar[idxx]);
                        INT2(H_vv.v, osc.mu*dvz*(1.0 - v*vv), a->v[idxx], a_bar_conj);

                        CONJ2(a_conj, a->v[idxx]);
                        INT2(H_vv.v_bar, -osc.mu*dvz*(1.0-v*vv), a_conj, a->v_bar[idxx]);
                    }
                    MAT2_ADD(H.v, 1.0, H_vac.v, 1.0, H_vv.v);
                    MAT2_ADD(H.v_bar, 1.0, H_vac.v_bar, 1.0, H_vv.v_bar);
#endif                
                    
                    COMMUTATOR2(f_osc.v, I, a->v[idx], H.v); 
                    COMMUTATOR2(f_osc.v_bar, I, a->v_bar[idx], H.v_bar); 
                    
                    MAT2_ADD(f->v[idx], 1.0, f_adv.v, 1.0, f_osc.v);
                    MAT2_ADD(f->v_bar[idx], 1.0, f_adv.v_bar, 1.0, f_osc.v_bar);

                    assert((f->v[idx]._12 == conj(f->v[idx]._21) ));
                }
            }
            break;
        case 2:
#pragma omp parallel for
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue;
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                        
                            size_t idx = index(sys.dim, jx, jz, kx, kz);
                            op2_t *a_ptr    = &a->v[idx];  
                            op2_t *a_bptr   = &a->v_bar[idx];
                            pair2_t f_adv, f_osc;
                            pair2_t H, H_vac, H_vv;
                    
                            memset(&f_adv.v, 0, sizeof(op2_t));
                            memset(&f_adv.v_bar, 0, sizeof(op2_t));
                            memset(&f_osc.v, 0, sizeof(op2_t));
                            memset(&f_osc.v_bar, 0, sizeof(op2_t));
                            memset(&H.v, 0, sizeof(op2_t));
                            memset(&H.v_bar, 0, sizeof(op2_t));
#ifdef ADV 
                            /*	advection	*/
                            double factor1 = -vx[kx]/(12.0*dx);
                            double factor2 = -vz[kz]/(12.0*dz);
                            ADV2(&f_adv.v, factor1, a_ptr, 1);
                            ADV2(&f_adv.v, factor2, a_ptr, sys.Nx+2*gx );
                            ADV2(&f_adv.v_bar, factor1, a_bptr, 1);
                            ADV2(&f_adv.v_bar, factor2, a_bptr, sys.Nx + 2*gx);
#endif

#ifdef KO
/*	Kreiss-Oliger dissipation (3-rd order)	*/
                            double ko_eps = -1e-4/(dz*16);
                            KO2(&f_adv.v, ko_eps, a_ptr, 1);
                            KO2(&f_adv.v, ko_eps, a_ptr, sys.Nx + 2*gx );
                            KO2(&f_adv.v_bar, ko_eps, a_bptr, 1);
                            KO2(&f_adv.v_bar, ko_eps, a_bptr, sys.Nx + 2*gx);
#endif

                            /*	oscillation	*/
#ifdef VACUUM
                            VAC2(H, osc.theta);
#elif INT_VV
                            for (int kkz = 0; kkz < sys.Nvz; ++kkz){
                                for (int kkx = 0; kkx < sys.Nvx; ++kkx){
                                    size_t idxx = index(sys.dim, jx, jz, kkx, kkz);
                                    op2_t a_conj, a_bar_conj;
                                    CONJ2(a_bar_conj, a->v_bar[idxx]);
                                    INT2(H.v, osc.mu*dvz*dvx*(1.0 - (vz[kz]*vz[kkz] + vx[kx]*vx[kkx])), a->v[idxx], a_bar_conj);

                                    CONJ2(a_conj, a->v[idxx]);
                                    INT2(H.v_bar, -osc.mu*dvz*dvx*(1.0-(vz[kz]*vz[kkz] + vx[kx]*vx[kkx])), a_conj, a->v_bar[idxx]);
                                }
                            } 
#elif VAC_VV
                            VAC2(H_vac, osc.theta);
                            memset(&H_vv.v, 0, sizeof(op2_t));
                            memset(&H_vv.v_bar, 0, sizeof(op2_t));
                            for (int kkz = 0; kkz < sys.Nvz; ++kkz){
                                for (int kkx = 0; kkx < sys.Nvx; ++kkx){
                                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                                        continue;
                                    size_t idxx = index(sys.dim, jx, jz, kkx, kkz);
                                    op2_t a_conj, a_bar_conj;
                                    CONJ2(a_bar_conj, a->v_bar[idxx]);
                                    INT2(H_vv.v, osc.mu*dvz*dvx*(1.0 - (vz[kz]*vz[kkz] + vx[kx]*vx[kkx])), a->v[idxx], a_bar_conj);

                                    CONJ2(a_conj, a->v[idxx]);
                                    INT2(H_vv.v_bar, -osc.mu*dvz*dvx*(1.0-(vz[kz]*vz[kkz] + vx[kx]*vx[kkx])), a_conj, a->v_bar[idxx]);
                                }
                            } 

                            MAT2_ADD(H.v, 1.0, H_vac.v, 1.0, H_vv.v);
                            MAT2_ADD(H.v_bar, 1.0, H_vac.v_bar, 1.0, H_vv.v_bar);
#endif                
                    
                            COMMUTATOR2(f_osc.v, I, a->v[idx], H.v); 
                            COMMUTATOR2(f_osc.v_bar, I, a->v_bar[idx], H.v_bar); 
                    
                            MAT2_ADD(f->v[idx], 1.0, f_adv.v, 1.0, f_osc.v);
                            MAT2_ADD(f->v_bar[idx], 1.0, f_adv.v_bar, 1.0, f_osc.v_bar);

                            assert((f->v[idx]._12 == conj(f->v[idx]._21) ));
                        }
                    }
                }
            }
            break;
        case 3:
            break;
    
    }   
 
}

void Neu::Qke2::set_init(field2_t *start) {
    switch (sys.dim) {
        case 1:
#pragma omp parallel for
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    size_t idx = index(sys.dim, j, k);
                    double e = Gaussian(1, z[j]);
                    //double g_ = g(vz[k], 0.6)/2.0;
                    //double g_bar = g(vz[k], 0.53)/2.0;
                    bzero(&start->v[idx], sizeof(op2_t));
                    bzero(&start->v_bar[idx], sizeof(op2_t));
                    //start->v[idx]._11 = g_*2*sqrt(1 - e*e);
                    //start->v[idx]._22 = 0;
                    //start->v[idx]._12 = g_*e;
                    //start->v[idx]._21 = g_*e;

                    //start->v_bar[idx]._11 = osc.alpha*g_bar*2*sqrt(1 - e*e);
                    //start->v_bar[idx]._22 = 0;
                    //start->v_bar[idx]._12 = osc.alpha*g_bar*e;
                    //start->v_bar[idx]._21 = osc.alpha*g_bar*e;
                    start->v[idx]._11 = e;
                    start->v[idx]._22 = 0;
                    start->v[idx]._12 = 0;
                    start->v[idx]._21 = start->v[idx]._12;
                      
                    start->v_bar[idx]._11 = e;
                    start->v_bar[idx]._22 = 0;
                    start->v_bar[idx]._12 = 0;
                    start->v_bar[idx]._21 = start->v_bar[idx]._12;

                }
            }
            break;
        case 2:
#pragma omp parallel for
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue;
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            size_t idx = index(sys.dim, jx, jz, kx, kz);
                            double e = Gaussian(2, x[jx], z[jz]);
                            //if (jz == 50 && jx == 50){
                            //    printf("%ld z[%d] = %lf, x[%d] = %lf\n", idx, jz, z[jz], jx, x[jx]);
                            //    printf("%lf \n", e);
                            //}
                            
                            //double gz = g(vz[kz], 0.6)/2.0;
                            //double gx = g(vx[kx], 0.6)/2.0;
                            //double gz_bar = g(vz[kz], 0.53)/2.0;
                            //double gx_bar = g(vx[kx], 0.53)/2.0;
                           
                            bzero(&start->v[idx], sizeof(op2_t));  
                            bzero(&start->v_bar[idx], sizeof(op2_t));
                            
                            //start->v[idx]._11 = gz*gx*2*sqrt(1-pow(osc.A*ex*ez,2));
                            //start->v[idx]._22 = 0;
                            //start->v[idx]._12  = gz*gx*osc.A*ex*ez;
                            //start->v[idx]._21 = start->v[idx]._12;
                            //
                            //start->v_bar[idx]._11 = osc.alpha*gz_bar*gx_bar*2*sqrt(1-pow(osc.A*ex*ez,2));
                            //start->v_bar[idx]._22 = 0;
                            //start->v_bar[idx]._12  = osc.alpha*gz_bar*gx_bar*osc.A*ex*ez;
                            //start->v_bar[idx]._21 = start->v_bar[idx]._12;
                            
                            start->v[idx]._11 = e;
                            start->v[idx]._22 = 0;
                            start->v[idx]._12 = 0;
                            start->v[idx]._21 = start->v[idx]._12;
                            
                            start->v_bar[idx]._11 = e;
                            start->v_bar[idx]._22 = 0;
                            start->v_bar[idx]._12 = 0;
                            start->v_bar[idx]._21 = start->v_bar[idx]._12;
                        }
                    } 
                }
            }
            break;
        case 3:
            break;
    }

}


void Neu::Qke2::field_to_pol(pol_t *P, const field2_t *fd) {
    switch (sys.dim) {
        case 1:
#pragma omp parallel for
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    int idx = index(sys.dim, j, k);
                    double g_ = g(vz[k], 0.6)/2.0;
                    double g_bar = g(vz[k], 0.53)/2.0;
                    
                    P->v[idx]._1_ = creal(fd->v[idx]._12)/g_;
                    P->v[idx]._2_ = cimag(fd->v[idx]._12)/g_;
                    P->v[idx]._3_ = creal(fd->v[idx]._11 - fd->v[idx]._22)/(2*g_);
                    P->v_bar[idx]._1_ = creal(fd->v_bar[idx]._12)/g_;
                    P->v_bar[idx]._2_ = cimag(fd->v_bar[idx]._12)/g_;
                    P->v_bar[idx]._3_ = creal(fd->v_bar[idx]._11 - fd->v_bar[idx]._22)/(2*g_);
                }
            }
            break;
        case 2:
#pragma omp parallel for 
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    /*	spherical symmetry ??	*/
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue;
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            int idx = index(sys.dim, jx, jz, kx, kz);
                            double gz = g(vz[kz], 0.6)/2.0;
                            double gx = g(vx[kx], 0.6)/2.0;
                            double gz_bar = g(vz[kz], 0.53)/2.0;
                            double gx_bar = g(vx[kx], 0.53)/2.0;
                           

                            P->v[idx]._1_ = creal(fd->v[idx]._12)/(gz*gx);
                            P->v[idx]._2_ = cimag(fd->v[idx]._12)/(gz*gx);
                            P->v[idx]._3_ = creal(fd->v[idx]._11 - fd->v[idx]._22)/(2*gz*gx);
                            P->v_bar[idx]._1_ = creal(fd->v_bar[idx]._12)/(gz*gx);
                            P->v_bar[idx]._2_ = cimag(fd->v_bar[idx]._12)/(gz*gx);
                            P->v_bar[idx]._3_ = creal(fd->v_bar[idx]._11 - fd->v_bar[idx]._22)/(2*gz*gx);
                        }
                    }
                }
            }
            break;
        case 3:
            break;
    }
    

}

double Neu::Qke2::Gaussian(int dim, ...){
    double x, y, z;
    va_list valist;
    va_start(valist, dim);
    if (dim == 1){
        z = va_arg(valist, double);
        va_end(valist);
        return exp(-pow(z-osc.z_0, 2)/(2*pow(osc.sigma,2)));
    }
    else if (dim == 2){
        x = va_arg(valist, double);
        z = va_arg(valist, double);
        va_end(valist); 
        return exp(-(pow(x-osc.z_0,2)+pow(z-osc.z_0,2))/(2*pow(osc.sigma,2)));
    }
    else
       return 0;
}

void Neu::Qke2::vac_exact(){
    double a = pow(cos(osc.theta), 4) + 0.5*pow(sin(2*osc.theta), 2)*cos(2*dt*osc.n) + pow(sin(osc.theta), 4);
    cnum b = I*sin(osc.n*dt)*sin(2*osc.theta)*(pow(cos(osc.theta),2)*(cos(osc.n*dt) + I*sin(osc.n*dt)) + pow(sin(osc.theta), 2)*(cos(osc.n*dt)-I*sin(osc.n*dt)));
    double c = pow(sin(osc.n*dt)*sin(2*osc.theta), 2);
    switch (sys.dim){
        case 1:
#pragma omp parallel for shared(a, b, c)
            for (int k = 0; k < sys.Nvz; ++k){
                double z_cent = vz[k]*(dt*(osc.n%NTz[k]));
                double z_prd  = z_cent - sgn(vz[k])*(sys.z2-sys.z1);
                double val;
                for (int j = 0; j < sys.Nz; ++j){
                    int idx = index(sys.dim, j, k);
                    //rho_ee = a;
                    if (fabs(z[j] - z_cent) < fabs(z[j] - z_prd)){
                        val = Gaussian(1, z[j] - z_cent);
                        rho_exact.v[idx]._11 = a*val;
                        rho_exact.v[idx]._12 = b*val;
                        rho_exact.v[idx]._21 = conj(rho_exact.v[idx]._12);
                        rho_exact.v[idx]._22 = c*val;
                    }
                    else{
                        val = Gaussian(1, z[j] - z_prd);
                        rho_exact.v[idx]._11 = a*val;
                        rho_exact.v[idx]._12 = b*val;
                        rho_exact.v[idx]._21 = conj(rho_exact.v[idx]._12);
                        rho_exact.v[idx]._22 = c*val;
                    }
                }
            }
            break;
        case 2:
#pragma omp parallel for shared(a, b, c)
           for (int kz = 0; kz < sys.Nvz; ++kz){
                double z_cent = vz[kz]*(dt*(osc.n%NTz[kz]));
                double z_prd  = z_cent - sgn(vz[kz])*(sys.z2-sys.z1);
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue;
                    double x_cent = vx[kx]*(dt*(osc.n%NTx[kx]));
                    double x_prd  = x_cent -sgn(vx[kx])*(sys.x2-sys.x1);
                    for (int jz = 0; jz < sys.Nz; jz+=1){
                        for (int jx = 0; jx < sys.Nx; jx+=1){
                            int idx = index(sys.dim, jx, jz, kx, kz);
                            //rho_ee = a;
                            double r1 = fabs(x[jx] - x_cent);
                            double r2 = fabs(x[jx] - x_prd);
                            double r3 = fabs(z[jz] - z_cent);
                            double r4 = fabs(z[jz] - z_prd);
                            if (r1 < r2 && r3 < r4){
                                rho_exact.v[idx]._11=a*Gaussian(2,x[jx]-x_cent,z[jz]-z_cent);
                            }
                            else if (r1 >= r2 && r3 < r4){
                                rho_exact.v[idx]._11=a*Gaussian(2,x[jx]-x_prd,z[jz]-z_cent);
                            }
                            else if (r1 < r2 && r3 >= r4){
                                rho_exact.v[idx]._11=a*Gaussian(2,x[jx]-x_cent,z[jz]-z_prd);
                            }
                            else{
                                rho_exact.v[idx]._11=a*Gaussian(2,x[jx]-x_prd,z[jz]-z_prd);
                            } 
                        }
                    }
                }
            }
            break;
        case 3:
            break; 
    }
}


void Neu::Qke2::output_rho(const char *filename, field2_t *rho){
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (!file){
       std::cout << "*** Open fails: " << filename << std::endl;
    }

 
    double error;
    switch (sys.dim){
        case 1:
            for (int k = 0; k < sys.Nvz; ++k){
                error  = 0.0;
                for (int j = 0; j < sys.Nz; j+=1){
                    int idx = index(sys.dim, j, k);
                    
                    //file << vz[k] << " " << z[j] << " " << creal(rho->v[idx]._11) << " " << creal(rho->v[idx]._12) << " " << cimag(rho->v[idx]._12) << " " << creal(rho->v[idx]._22) << std::endl;
                    error += pow(creal(rho->v[idx]._11)- creal(rho_exact.v[idx]._11), 2);
                    
                    file << vz[k] << " " << z[j] << " " << creal(rho->v[idx]._11) << " " << creal(rho_exact.v[idx]._11) << std::endl;
                }
                error = sqrt(error);
                file << error << std::endl;
            } 
            break;
        case 2:
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue;
                    error = 0.0;
                    for (int jz = 0; jz < sys.Nz; jz+=1){
                        for (int jx = 0; jx < sys.Nx; jx+=1){
                            int idx = index(sys.dim, jx, jz, kx, kz);
                            
                            error += dz*dx*pow(creal(rho->v[idx]._11-rho_exact.v[idx]._11),2);
                            //file << vz[kz] << " " << vx[kx] << dt*" " << z[jz] << " " << x[jx] << " " <<  creal(rho->v[idx]._11) << " " << creal(rho->v[idx]._12) << " " << cimag(rho->v[idx]._12) << " " << creal(rho->v[idx]._22) << " " << rho_ee << std::endl;
                            file << vz[kz] << " " << vx[kx] << " " << z[jz] << " " << x[jx] << " " <<  creal(rho->v[idx]._11) << " " << creal(rho_exact.v[idx]._11) << " " << std::endl;
                        }
                    }
                    error = sqrt(error);
                    file << error << std::endl;
                }
            }
            break;
        case 3:
            break;
    }
    file.close();
}

void Neu::Qke2::output_P(const char *filename, pol_t *P){
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (!file){
       std::cout << "*** Open fails: " << filename << std::endl;
    }

    switch (sys.dim){
        case 1:
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; j+=1){
                    int idx = index(sys.dim, j, k);
                    file << vz[k] << " " << z[j] << " " << P->v[idx]._1_ << " " << P->v[idx]._2_ << " " << P->v[idx]._3_ << std::endl;
                }
            } 
            break;
        case 2:
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int kx = 0; kx < sys.Nvx; ++kx){
                    if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                        continue;
                    for (int jz = 0; jz < sys.Nz; jz+=1){
                        for (int jx = 0; jx < sys.Nx; jx+=1){
                            int idx = index(sys.dim, jx, jz, kx, kz);
                            file << vz[kz] << " " << vx[kx] << " " << z[jz] << " " << x[jx] << " " << P->v[idx]._1_ << " " << P->v[idx]._2_ << " " << P->v[idx]._3_ << std::endl;

                        }
                    }
                }
            }
            break;
        case 3:
            break;
    }
    file.close();
}


void Neu::Qke2::run(){
    char file[0x100];
    printf("dt = %lf\n", dt);
    switch (sys.dim){
        case 1:
           printf("dz = %lf\n", dz); 
           break;
        case 2:
            printf("dx = %lf dz = %lf\n", dx, dz);
            break;
    }   

    set_init(&rho);
    for (int i = 0; i <= sys.Nt; ++i){
        printf("%d\n", i);
        osc.n = i;
        if (!(i & 1)) {
            if (i%100 == 0){
                //field_to_pol(&P, &rho);
                //sprintf(file, "P_%05d.txt", i); 
                //output_P(file, &P);
                vac_exact();
                sprintf(file, "rho_%05d.txt", i);
                output_rho(file, &rho); 
            }  
            rk4(&rho, &rho_next);
            
        }
        else {
            rk4(&rho_next, &rho);
        }
    } 
}

double Neu::g(double u, double xi){
    return sqrt(2.0/M_PI)*exp(-pow(u-1.0,2)/(2.0*xi*xi))/(xi*erf(sqrt(2)/xi));
}









