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
#include "objects.hpp"
#include "neu.hpp"

#define sgn(v) ((v>0) - (v<0)) //sign function in branchless implementation



#define DECLARE_UTIL(flavor)    \
    static inline void ASSIGN##flavor(op##flavor##_t *y, op##flavor##_t *x) { \
        FOR(m, flavor*flavor){                \
            y->element[m]= x->element[m];    \
        }   \
    }\
    static inline void ADV##flavor(op##flavor##_t *y, cnum a, op##flavor##_t *x, size_t N){                   \
        FOR(m, flavor*flavor)             \
            y->element[m] += a*(8.0*(x[1*N].element[m] - x[-1*N].element[m]) - (x[2*N].element[m]- x[-2*N].element[m]));    \
    }                                                       \
                                                            \
    static inline void KO##flavor(op##flavor##_t *y, cnum a, op##flavor##_t *x, size_t N){                    \
        FOR(m, flavor*flavor)             \
            y->element[m] += a*((x[-2*N].element[m] + x[2*N].element[m]) - 4*(x[-1*N].element[m] + x[-1*N].element[m]) + 6*x[0].element[m]);\
    }                           \
                                \
    static inline void MAT##flavor##_ADD(op##flavor##_t *C, cnum a, op##flavor##_t *A, cnum b, op##flavor##_t *B){           \
        FOR(m, flavor*flavor)             \
            C->element[m] = a*A->element[m] + b*B->element[m]; \
    }                                                       \
    static inline void INT##flavor(op##flavor##_t *A, cnum a, const op##flavor##_t *B, const op##flavor##_t *C){                          \
        FOR(m, flavor*flavor)                                         \
            A->element[m] += a*(B->element[m] - C->element[m]);    \
    }




DECLARE_UTIL(2)
DECLARE_UTIL(3) 
#undef DECLARE_UTIL

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
        {"Nphi", required_argument, NULL, 13},
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
                sscanf(optarg, "%d", &num_sys->Nphi);
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
                //sscanf(optarg, "%d", &num_sys->Nvy);
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
                printf("this is A : %.8lf\n", phys_sys->A);
                break;
            case 24:
                sscanf(optarg, "%lf", &phys_sys->z_0);
                break;
            case 25:
                sscanf(optarg, "%lf", &phys_sys->sigma);
                printf("This is sigma = %lf\n", phys_sys->sigma);
                break;
        }
    }
}



#define VAC2(H, theta){ \
    H.v.element[0] += -cos(2.0*theta);  \
    H.v.element[1] += sin(2.0*theta);   \
    H.v.element[2] += sin(2.0*theta);   \
    H.v.element[3] += cos(2.0*theta);   \
    H.v_bar = H.v;              \
}

#define CONJ2(A, B){ \
    A.element[0] = B.element[0];          \
    A.element[1] = conj(B.element[1]);    \
    A.element[2] = conj(B.element[2]);    \
    A.element[3] = B.element[3];          \
}

#define COMMUTATOR2(C, a, A, B){ \
    C.element[0] += a*(A.element[1]*B.element[2] - A.element[2]*B.element[1]);\
    C.element[1] += a*((A.element[0] - A.element[3])*B.element[1] + A.element[1]*(B.element[3] - B.element[0]));  \
    C.element[2] += a*((A.element[3] - A.element[0])*B.element[2] + A.element[2]*(B.element[0] - B.element[3]));  \
    C.element[3] += a*(A.element[2]*B.element[1] - A.element[1]*B.element[2]);\
}

void Neu::Qke2::rk4(field2_t *curr, field2_t *next) {
    field2_t *f0, *f1, *f2, *f3;
    field2_t *curr1, *curr2, *curr3; 
    for (size_t i = 0; i < SIZE; ++i){
        FOR(m, 4){
            next->v[i].element[m] = 0.0;
            next->v_bar[i].element[m] = 0.0;
        }
    }
    f0 = make_field2(SIZE);
    f1 = make_field2(SIZE);
    f2 = make_field2(SIZE);
    f3 = make_field2(SIZE);
    curr1 = make_field2(SIZE);
    curr2 = make_field2(SIZE);
    curr3 = make_field2(SIZE);
    

    /*	Step 1	*/    
    periodic_ghost_zone(curr);
    compute_F(f0, curr);
    fieldadd(curr1, 1.0, curr, 0.5*dt, f0);

    /*	Step 2	*/
    periodic_ghost_zone(curr1);
    compute_F(f1, curr1);
    fieldadd(curr2, 1.0, curr, 0.5*dt, f1);
    
    /*	Step 3	*/
    periodic_ghost_zone(curr2);
    compute_F(f2, curr2);
    fieldadd(curr3, 1.0, curr, dt, f2);

    /*	Step 4	*/
    periodic_ghost_zone(curr3);
    compute_F(f3, curr3);

    fieldadd(next, 1.0, next, 1.0, curr);
    fieldadd(next, 1.0, next, dt/6.0, f0);
    fieldadd(next, 1.0, next, dt/3.0, f1);
    fieldadd(next, 1.0, next, dt/3.0, f2);
    fieldadd(next, 1.0, next, dt/6.0, f3);

   
    free_field2(f0);
    free_field2(f1);
    free_field2(f2);
    free_field2(f3);
    free_field2(curr1);
    free_field2(curr2);
    free_field2(curr3);


}


void Neu::Qke2::fieldadd(field2_t *C, const cnum a, field2_t *A, const cnum b, field2_t *B) {
    switch (sys.dim){
        case 1:
#pragma omp parallel for 
           for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    size_t idx = index(sys.dim, j, k);

                    MAT2_ADD(&C->v[idx], a, &A->v[idx], b, &B->v[idx]);
                    MAT2_ADD(&C->v_bar[idx], a, &A->v_bar[idx], b, &B->v_bar[idx]);
                }
            }
            break;
        case 2:
#pragma omp parallel for
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            size_t idx = index(sys.dim, jx, jz, f, kz);
                            MAT2_ADD(&C->v[idx], a, &A->v[idx], b, &B->v[idx]);
                            MAT2_ADD(&C->v_bar[idx], a, &A->v_bar[idx], b, &B->v_bar[idx]);
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
                    size_t idx_l = index(sys.dim, -j-1, k);
                    size_t idx_r = index(sys.dim, sys.Nz+j, k);
                    size_t idx;
                    idx = index(sys.dim, sys.Nz-j-1, k);
                    ASSIGN2(&A->v[idx_l], &A->v[idx]);
                    ASSIGN2(&A->v_bar[idx_l], &A->v_bar[idx]);
                    //A->v[idx] = A->v[index(sys.dim, sys.Nz-j-1,k)];
                    //A->v_bar[idx] = A->v_bar[index(sys.dim, sys.Nz-j-1,k)];
                    /*	lower side	*/
                    idx = index(sys.dim, j,k);
                    ASSIGN2(&A->v[idx_r], &A->v[idx]);
                    ASSIGN2(&A->v_bar[idx_r], &A->v_bar[idx]);

                    //A->v[idx] = A->v[index(sys.dim, j, k)];
                    //A->v_bar[idx] = A->v_bar[index(sys.dim, j, k)];
                }
            }
            break;
        case 2:
#pragma omp parallel for 
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2), + pow(phi[kx],2) > 1)
                    //    continue;
                    for (int jz = 0; jz < gz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            /*	z	*/ 
                            size_t idx = index(sys.dim, jx, -jz-1, f, kz);
                            size_t id1 = index(sys.dim, jx,sys.Nz-jz-1,f,kz); 
                            size_t id2 = index(sys.dim, jx,jz, f, kz);
                            ASSIGN2(&A->v[idx], &A->v[id1]);
                            ASSIGN2(&A->v_bar[idx], &A->v_bar[id1]);
                            //A->v[idx] = A->v[id1];
                            //A->v_bar[idx] = A->v_bar[id1];
                            
                            idx = index(sys.dim, jx, sys.Nz+jz, f, kz);
                            ASSIGN2(&A->v[idx], &A->v[id2]);
                            ASSIGN2(&A->v_bar[idx], &A->v_bar[id2]);
                            //A->v[idx] = A->v[id2];
                            //A->v_bar[idx] = A->v_bar[id2];
                        }
                    }
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < gx; ++jx){
                            /*	x	*/
                            size_t idx = index(sys.dim, -jx-1, jz, f, kz);
                            size_t id1 = index(sys.dim, sys.Nx-jx-1,jz,f,kz);
                            size_t id2 = index(sys.dim, jx, jz, f, kz);
                            ASSIGN2(&A->v[idx], &A->v[id1]);
                            ASSIGN2(&A->v_bar[idx], &A->v_bar[id1]);

                            //A->v[idx] = A->v[id1];
                            //A->v_bar[idx] = A->v_bar[id1];
                            
                            idx = index(sys.dim, sys.Nx+jx, jz, f, kz);
                            ASSIGN2(&A->v[idx], &A->v[id2]);
                            ASSIGN2(&A->v_bar[idx], &A->v_bar[id2]);
                            //A->v[idx] = A->v[id2];
                            //A->v_bar[idx] = A->v_bar[id2];

                        }
                    }
                }
            }
            break;
        case 3:
            break;
    } 

}
void Neu::Qke2::compute_F(field2_t *fld, const field2_t *a) {
    switch (sys.dim){
        case 1:
//#pragma omp parallel for
#pragma omp parallel for
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    size_t idx = index(sys.dim, j, k);
                    double v = vz[k];
                    op2_t *a_ptr    = &a->v[idx];  
                    op2_t *a_bptr   = &a->v_bar[idx];

                    pair2_t f, H;
                    zero_pair2(&f);
                    zero_pair2(&H);
#ifdef ADV 
/*	Advection	*/
                    double factor = -v/(12.0*dz);
                    ADV2(&f.v, factor, a_ptr, 1);
                    ADV2(&f.v_bar, factor, a_bptr, 1);

#endif

#ifdef KO
/*	Kreiss-Oliger dissipation (3-rd order)	*/
                    double ko_eps = -1e-4/(dz*16);
                    KO2(&f.v, ko_eps, a_ptr, 1);
                    KO2(&f.v_bar, ko_eps, a_bptr, 1);

#endif
                    /*	oscillation	*/
#ifdef VACUUM
                    VAC2(H, osc.theta);
#elif INT_VV
                    for (int kk = 0; kk < sys.Nvz; ++kk){
                        size_t idxx = index(sys.dim, j, kk);
                        double vv = vz[kk];
                        pair2_t a_conj;
                        CONJ2(a_conj.v_bar, a->v_bar[idxx]);
                        INT2(&H.v, osc.mu*dvz*(1.0 - v*vv), &a->v[idxx], &a_conj.v_bar);

                        CONJ2(a_conj.v, a->v[idxx]);
                        INT2(&H.v_bar, -osc.mu*dvz*(1.0-v*vv), &a_conj.v, &a->v_bar[idxx]);
                    } 

#elif VAC_VV
                    VAC2(H, osc.theta);
                    for (int kk = 0; kk < sys.Nvz; ++kk){
                        size_t idxx = index(sys.dim, j, kk);
                        double vv = vz[kk];
                        pair2_t a_conj;
                        CONJ2(a_conj->v_bar, a->v_bar[idxx]);
                        INT2(&H.v, osc.mu*dvz*(1.0 - v*vv), &a->v[idxx], &a_conj->v_bar);

                        CONJ2(a_conj->v, a->v[idxx]);
                        INT2(&H.v_bar, -osc.mu*dvz*(1.0-v*vv), &a_conj->v, &a->v_bar[idxx]);
                    }
#endif                
                    
                    COMMUTATOR2(f.v, I, a->v[idx], H.v); 
                    COMMUTATOR2(f.v_bar, I, a->v_bar[idx], H.v_bar); 

                    ASSIGN2(&fld->v[idx], &f.v);
                    ASSIGN2(&fld->v_bar[idx], &f.v_bar);
                    
                    assert((fld->v[idx].element[1] == conj(fld->v[idx].element[2]) ));
                }
            }
            break;
        case 2:
            
#pragma omp parallel for
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                    //    continue;
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                        
                            size_t idx = index(sys.dim, jx, jz, f, kz);
                            op2_t *a_ptr    = &a->v[idx];  
                            op2_t *a_bptr   = &a->v_bar[idx];
                            pair2_t _f, H;                  
                            zero_pair2(&_f);
                            zero_pair2(&H);
                             
#ifdef ADV 
/*	Advection	*/
                            double factor1 = -vx[kz*sys.Nphi+f]/(12.0*dx);
                            double factor2 = -vz[kz]/(12.0*dz);
                            ADV2(&_f.v, factor1, a_ptr, 1);
                            ADV2(&_f.v, factor2, a_ptr, sys.Nx + 2*gx );
                            ADV2(&_f.v_bar, factor1, a_bptr, 1);
                            ADV2(&_f.v_bar, factor2, a_bptr, sys.Nx + 2*gx);
#endif

#ifdef KO
/*	Kreiss-Oliger dissipation (3-rd order)	*/
                            double ko_eps = -1e-4/(dz*16);
                            KO2(&_f.v, ko_eps, a_ptr, 1);
                            KO2(&_f.v, ko_eps, a_ptr, sys.Nx + 2*gx );
                            KO2(&_f.v_bar, ko_eps, a_bptr, 1);
                            KO2(&_f.v_bar, ko_eps, a_bptr, sys.Nx + 2*gx);
#endif

                            /*	oscillation	*/
#ifdef VACUUM
                            VAC2(H, osc.theta);
#elif INT_VV
                            for (int kkz = 0; kkz < sys.Nvz; ++kkz){
                                for (int ff = 0; ff < sys.Nphi; ++ff){
                                    size_t idxx = index(sys.dim,jx,jz,ff,kkz);
                                    pair2_t a_conj;
                                    CONJ2(a_conj.v_bar, a->v_bar[idxx]);
                                    INT2(&H.v, osc.mu*dvz*dphi/(2*M_PI)*(1.0 - (vz[kz]*vz[kkz] + vx[kz*sys.Nphi+f]*vx[kkz*sys.Nphi+ff] + vy[kz*sys.Nphi+f]*vy[kkz*sys.Nphi+ff])),&a->v[idxx], &a_conj.v_bar);
                                    //INT2(H.v, osc.mu*dvz*dphi/(2*M_PI)*(1.0 - (vz[kz]*vz[kkz] )), a->v[idxx], a_bar_conj);


                                    CONJ2(a_conj.v, a->v[idxx]);
                                    INT2(&H.v_bar, -osc.mu*dvz*dphi/(2*M_PI)*(1.0 - (vz[kz]*vz[kkz] + vx[kz*sys.Nphi+f]*vx[kkz*sys.Nphi+ff] + vy[kz*sys.Nphi+f]*vy[kkz*sys.Nphi+ff])),&a_conj.v, &a->v_bar[idxx]);
                                    //INT2(H.v_bar, -osc.mu*dvz*dphi/(2*M_PI)*(1.0-(vz[kz]*vz[kkz] )), a_conj, a->v_bar[idxx]);
                                }
                            } 
#elif VAC_VV
                            VAC2(H, osc.theta);
                            //memset(&H.v, 0, sizeof(op2_t));
                            //memset(&H_vv.v_bar, 0, sizeof(op2_t));
                            for (int kkz = 0; kkz < sys.Nvz; ++kkz){
                                for (int ff = 0; ff < sys.Nphi; ++ff){
                                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                                    //    continue;
                                    size_t idxx = index(sys.dim,jx,jz,ff,kkz);
                                    pair2_t a_conj;
                                    CONJ2(a_conj.v_bar, a->v_bar[idxx]);
                                    INT2(&H->v, osc.mu*dvz*dphi/(2*M_PI)*(1.0 - (vz[kz]*vz[kkz] + vx[kz*sys.Nphi+f]*vx[kkz*sys.Nphi+ff] + vy[kz*sys.Nphi+f]*vy[kkz*sys.Nphi+ff])), &a->v[idxx], &a_conj.v_bar);


                                    CONJ2(a_conj.v, a->v[idxx]);
                                    INT2(&H->v_bar, -osc.mu*dvz*dphi/(2*M_PI)*(1.0-(vz[kz]*vz[kkz] + vx[kz*sys.Nphi+f]*vx[kkz*sys.Nphi+ff] + vy[kz*sys.Nphi+f]*vy[kkz*sys.Nphi+ff])), &a_conj.v, &a->v_bar[idxx]);
                                }
                            } 

#endif                
                    
                            COMMUTATOR2(_f.v, I, a->v[idx], H.v); 
                            COMMUTATOR2(_f.v_bar, I, a->v_bar[idx], H.v_bar); 
                    
                            ASSIGN2(&fld->v[idx], &_f.v);
                            ASSIGN2(&fld->v_bar[idx], &_f.v_bar);

                            assert((fld->v[idx].element[1] == conj(fld->v[idx].element[2]) ));
                        
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
                    double g_ = g(vz[k], 0.6, 1)/2.0;
                    double g_bar = g(vz[k], 0.53, 1)/2.0;
                    //bzero(&start->v[idx], sizeof(op2_t));
                    //bzero(&start->v_bar[idx], sizeof(op2_t));
                    start->v[idx].element[0] = g_*(1 + sqrt(1-pow(osc.A*e,2)));
                    start->v[idx].element[3] = g_*(1 - sqrt(1-pow(osc.A*e,2)));
                    start->v[idx].element[1] = g_*osc.A*e;
                    start->v[idx].element[2] = g_*osc.A*e;

                    start->v_bar[idx].element[0] = osc.alpha*g_bar*(1+sqrt(1-pow(osc.A*e,2)));
                    start->v_bar[idx].element[3] = osc.alpha*g_bar*(1-sqrt(1-pow(osc.A*e,2)));

                    start->v_bar[idx].element[1] = osc.alpha*g_bar*osc.A*e;
                    start->v_bar[idx].element[2] = osc.alpha*g_bar*osc.A*e;
                    //start->v[idx]._11 = e;
                    //start->v[idx]._22 = 0;
                    //start->v[idx]._12 = 0;
                    //start->v[idx]._21 = start->v[idx]._12;
                    //  
                    //start->v_bar[idx]._11 = e;
                    //start->v_bar[idx]._22 = 0;
                    //start->v_bar[idx]._12 = 0;
                    //start->v_bar[idx]._21 = start->v_bar[idx]._12;

                }
            }
            break;
        case 2:
#pragma omp parallel for
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                    //    continue;
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            size_t idx = index(sys.dim, jx, jz, f, kz);
                            double e = Gaussian(1, z[jz]);
                            //if (jz == 50 && jx == 50){
                            //    printf("%ld z[%d] = %lf, x[%d] = %lf\n", idx, jz, z[jz], jx, x[jx]);
                            //    printf("%lf \n", e);
                            //}
                            
                            double gz = g(vz[kz], 0.6, 1)/2.0;
                            //double gx = g(vx[kz*sys.Nphi+f], 0.6, 0.8);
                            double gz_bar = g(vz[kz], 0.53, 1)/2.0;
                            //double gx_bar = g(vx[kz*sys.Nphi+f], 0.53, 0.8);
                           
                            //bzero(&start->v[idx], sizeof(op2_t));  
                            //bzero(&start->v_bar[idx], sizeof(op2_t));
                            //
                            start->v[idx].element[0] = gz*(1+sqrt(1-pow(osc.A*e,2)));
                            start->v[idx].element[3] = gz*(1-sqrt(1-pow(osc.A*e,2)));
                            start->v[idx].element[1]  = gz*osc.A*e;
                            start->v[idx].element[2] = start->v[idx].element[1];
                            //
                            start->v_bar[idx].element[0] = osc.alpha*gz_bar*(1+sqrt(1-pow(osc.A*e,2)));
                            start->v_bar[idx].element[3] = osc.alpha*gz_bar*(1-sqrt(1-pow(osc.A*e,2)));
                            start->v_bar[idx].element[1]  = osc.alpha*gz_bar*osc.A*e;
                            start->v_bar[idx].element[2] = start->v_bar[idx].element[1];
                            
                            //start->v[idx]._11 = e;
                            //start->v[idx]._22 = 0;
                            //start->v[idx]._12 = 0;
                            //start->v[idx]._21 = start->v[idx]._12;
                            
                            //start->v_bar[idx]._11 = e;
                            //start->v_bar[idx]._22 = 0;
                            //start->v_bar[idx]._12 = 0;
                            //start->v_bar[idx]._21 = start->v_bar[idx]._12;
                        }
                    } 
                }
            }
            break;
        case 3:
            break;
    }

}


void Neu::Qke2::field_to_pol(pol2_t *P, const field2_t *fd) {
    switch (sys.dim) {
        case 1:
#pragma omp parallel for
            for (int k = 0; k < sys.Nvz; ++k){
                for (int j = 0; j < sys.Nz; ++j){
                    int idx = index(sys.dim, j, k);
                    double g_ = g(vz[k], 0.6, 1)/2.0;
                    double g_bar = g(vz[k], 0.53, 1)/2.0;                    
                    
                    P->v[idx].data[0] = creal(fd->v[idx].element[1])/g_;
                    P->v[idx].data[1] = cimag(fd->v[idx].element[1])/g_;
                    P->v[idx].data[2] = creal(fd->v[idx].element[0] - fd->v[idx].element[3])/(2*g_);
                    P->v_bar[idx].data[0] = creal(fd->v_bar[idx].element[1])/g_;
                    P->v_bar[idx].data[1] = cimag(fd->v_bar[idx].element[1])/g_;
                    P->v_bar[idx].data[2] = creal(fd->v_bar[idx].element[0] - fd->v_bar[idx].element[3])/(2*g_);
                }
            }
            break;
        case 2:
#pragma omp parallel for 
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                    //    continue;
                    for (int jz = 0; jz < sys.Nz; ++jz){
                        for (int jx = 0; jx < sys.Nx; ++jx){
                            int idx = index(sys.dim, jx, jz, f, kz);
                            double gz = g(vz[kz], 0.6, 1)/2.0;
                            //double gx = g(vx[f], 0.6, 0.8);
                            double gz_bar = g(vz[kz], 0.53, 1)/2.0;
                            //double gx_bar = g(vx[f], 0.53, 0.8);
                           

                            P->v[idx].data[0] = creal(fd->v[idx].element[1])/(gz);
                            P->v[idx].data[1] = cimag(fd->v[idx].element[1])/(gz);
                            P->v[idx].data[2] = creal(fd->v[idx].element[0] - fd->v[idx].element[3])/(2*gz);
                            P->v_bar[idx].data[0] = creal(fd->v_bar[idx].element[1])/(gz);
                            P->v_bar[idx].data[1] = cimag(fd->v_bar[idx].element[1])/(gz);
                            P->v_bar[idx].data[2] = creal(fd->v_bar[idx].element[0] - fd->v_bar[idx].element[3])/(2*gz);
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
                        rho_exact->v[idx].element[0] = a*val;
                        rho_exact->v[idx].element[1] = b*val;
                        rho_exact->v[idx].element[2] = conj(rho_exact->v[idx].element[1]);
                        rho_exact->v[idx].element[3] = c*val;
                    }
                    else{
                        val = Gaussian(1, z[j] - z_prd);
                        rho_exact->v[idx].element[0] = a*val;
                        rho_exact->v[idx].element[1] = b*val;
                        rho_exact->v[idx].element[2] = conj(rho_exact->v[idx].element[1]);
                        rho_exact->v[idx].element[3] = c*val;
                    }
                }
            }
            break;
        case 2:
#pragma omp parallel for shared(a, b, c)
           for (int kz = 0; kz < sys.Nvz; ++kz){
                double z_cent = vz[kz]*(dt*(osc.n%NTz[kz]));
                double z_prd  = z_cent - sgn(vz[kz])*(sys.z2-sys.z1);
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                    //    continue;
                    double x_cent = vx[kz*sys.Nphi+f]*(dt*(osc.n%NTx[kz*sys.Nphi+f]));
                    double x_prd  = x_cent -sgn(vx[kz*sys.Nphi+f])*(sys.x2-sys.x1);
                    for (int jz = 0; jz < sys.Nz; jz+=1){
                        for (int jx = 0; jx < sys.Nx; jx+=1){
                            int idx = index(sys.dim, jx, jz, f, kz);
                            //rho_ee = a;
                            double r1 = fabs(x[jx] - x_cent);
                            double r2 = fabs(x[jx] - x_prd);
                            double r3 = fabs(z[jz] - z_cent);
                            double r4 = fabs(z[jz] - z_prd);
                            if (r1 < r2 && r3 < r4){
                                rho_exact->v[idx].element[0]=a*Gaussian(2,x[jx]-x_cent,z[jz]-z_cent);
                            }
                            else if (r1 >= r2 && r3 < r4){
                                rho_exact->v[idx].element[0]=a*Gaussian(2,x[jx]-x_prd,z[jz]-z_cent);
                            }
                            else if (r1 < r2 && r3 >= r4){
                                rho_exact->v[idx].element[0]=a*Gaussian(2,x[jx]-x_cent,z[jz]-z_prd);
                            }
                            else{
                                rho_exact->v[idx].element[0]=a*Gaussian(2,x[jx]-x_prd,z[jz]-z_prd);
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
                for (int j = 0; j < sys.Nz; j+=10){
                    int idx = index(sys.dim, j, k);
                    
                    //file << vz[k] << " " << z[j] << " " << creal(rho->v[idx]._11) << " " << creal(rho->v[idx]._12) << " " << cimag(rho->v[idx]._12) << " " << creal(rho->v[idx]._22) << std::endl;
                    error += pow(creal(rho->v[idx].element[0])- creal(rho_exact->v[idx].element[0]), 2);
                    
                    file << vz[k] << " " << z[j] << " " << creal(rho->v[idx].element[0]) << " " << creal(rho_exact->v[idx].element[0]) << std::endl;
                }
                error = sqrt(error);
                file << error << std::endl;
            } 
            break;
        case 2:
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                    //    continue;
                    error = 0.0;
                    for (int jz = 0; jz < sys.Nz; jz+=1){
                        for (int jx = 0; jx < sys.Nx; jx+=1){
                            int idx = index(sys.dim, jx, jz, f, kz);
                            
                            error += dz*dx*pow(creal(rho->v[idx].element[0]-rho_exact->v[idx].element[0]),2);
                            //file << vz[kz] << " " << vx[kx] << dt*" " << z[jz] << " " << x[jx] << " " <<  creal(rho->v[idx]._11) << " " << creal(rho->v[idx]._12) << " " << cimag(rho->v[idx]._12) << " " << creal(rho->v[idx]._22) << " " << rho_ee << std::endl;
                            file << vz[kz] << " " << vx[kz*sys.Nphi+f] << " " << z[jz] << " " << x[jx] << " " <<  creal(rho->v[idx].element[0]) << " " << creal(rho_exact->v[idx].element[0]) << " " << std::endl;
                 
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

void Neu::Qke2::output_P(const char *filename, const pol2_t *P){
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
                    file << vz[k] << " " << z[j] << " " << P->v[idx].data[0] << " " << P->v[idx].data[1] << " " << P->v[idx].data[2] << std::endl;
                }
            } 
            break;
        case 2:
            for (int kz = 0; kz < sys.Nvz; ++kz){
                for (int f = 0; f < sys.Nphi; ++f){
                    //if (pow(vz[kz],2) + pow(vx[kx],2) > 1)
                    //    continue;
                    for (int jz = 0; jz < sys.Nz; jz+=1){
                        for (int jx = 0; jx < sys.Nx; jx+=5){
                            int idx = index(sys.dim, jx, jz, f, kz);
                            file << vz[kz] << " " << vx[kz*sys.Nphi+f] << " " << z[jz] << " " << x[jx] << " " << P->v[idx].data[0] << " " << P->v[idx].data[1] << " " << P->v[idx].data[2] << std::endl;
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

    set_init(rho);
    for (int i = 0; i <= sys.Nt; ++i){
        printf("%d\n", i);
        osc.n = i;
        if (!(i & 1)) {
            if (i%10 == 0){
#ifdef VACUUM
                //vac_exact();
                //sprintf(file, "rho_%05d.txt", i);
                //output_rho(file, &rho);
                field_to_pol(P, rho);
                sprintf(file, "P_%05d_1dim.txt", i);
                output_P(file, P);
#elif INT_VV
                field_to_pol(P, rho);
                sprintf(file, "P_%05d_%ddim.txt", i, sys.dim);
                output_P(file, P);
#elif ADV
                field_to_pol(P, rho);
                sprintf(file, "P_%05d_%ddim.txt", i, sys.dim);
                output_P(file, P);
#endif 
            }  
            rk4(rho, rho_next);
            
        }
        else {
            rk4(rho_next, rho);
        }
    } 
}

double Neu::g(double u, double xi, double v){
    double N = sqrt(M_PI/2.0)*xi*(erf((v+1)/(sqrt(2)*xi))-erf((v-1)/(sqrt(2)*xi)));
    return exp(-pow(u-v,2)/(2.0*xi*xi))/N;
}
