#ifndef NEU_H
#define NEU_H
#include <stdio.h>
#include <stdarg.h>
#include <complex.h>
#include "objects.hpp"

typedef struct phys_sys_t {
    int     flavor;
    double  pmo = 1;
    double  theta = 33.3;      /*	mixing angle	*/
    double  mu = 1.0;         /*	interaction strength	*/
 
    int     n;                  /*	n-th time-iteration	*/ 
    int     init; 
    double  f_0; 
    double  alpha;
    double  beta;
    double  A, z_0, sigma;
} phys_sys_t;

typedef struct num_sys_t {
    int     dim;        /*	1 + dim. space-time dimensions	*/
    
    double  CFL; 
    int     Nt;
    
    double  x1, x2, y1, y2, z1, z2;
    int     Nx, Ny, Nz;
    
    int     Nphi,  Nvz;
} num_sys_t;

///*	2-flavor
// *	SU(2)	*/
//typedef struct operator2_t {
//    cnum _11; cnum _12;
//    cnum _21; cnum _22;
//} op2_t;
//typedef struct vector2_t {
//    double _1_;
//    double _2_;
//    double _3_;
//} vec2_t;
//
//typedef struct field2_t {
//    op2_t *v;
//    op2_t *v_bar;
//} field2_t;
//
//typedef struct pair2_t{
//    op2_t v;
//    op2_t v_bar;
//} pair2_t;
//
//typedef struct polarization_t {
//    vec2_t *v;
//    vec2_t *v_bar;
//} pol_t;
//
///*	3-flavor
// *	SU(3)	*/
//typedef struct operator3_t {
//    cnum _11; cnum _12; cnum _13;
//    cnum _21; cnum _22; cnum _23;
//    cnum _31; cnum _32; cnum _33;
//} op3_t;
//typedef struct vector3_t {
//    double _1_;
//    double _2_; 
//    double _3_; 
//    double _4_; 
//    double _5_; 
//    double _6_; 
//    double _7_; 
//    double _8_; 
//} vec3_t;
//
//typedef struct field3_t {
//    op3_t *v;
//    op3_t *v_bar;
//} field3_t;
//
//typedef struct pair3_t {
//    op3_t v;
//    op3_t v_bar;
//} pair3_t;
//
//
namespace Neu {
    void parsing(int argc, char *const *argv, phys_sys_t *a, num_sys_t *b);
    class Qke2 {
        protected:
            phys_sys_t  osc; 
            num_sys_t   sys;
            double      t;
            double      dt;
            double      v1, v2;
            size_t      SIZE;
            field2_t    *rho_exact; 
            field2_t    *rho;      
            field2_t    *rho_next;
            pol2_t      *P; 
            double  dx, dy, dz;
            double  dvz, dphi;
            double  *x, *z;
            double  *vz, *phi;
            double  *vx, *vy;
            int     gx, gy, gz;

            int     *NTx, *NTz;    //period of iterations for PBC
        public:    
            Qke2(const phys_sys_t _phys_sys, const num_sys_t _num_sys){
                osc    = _phys_sys;
                sys    = _num_sys;
                t = 0;
                v1 = -1.0;
                v2 = 1.0;
                gx = gy = gz = 2; 
                switch (sys.dim){
                    case 1:
                        SIZE = (size_t)(sys.Nz+2*gz)*sys.Nvz;
                        dz  = (sys.z2 - sys.z1)/sys.Nz;
                        dvz = (v2 - v1)/sys.Nvz;
                        z   = new double[sys.Nz];
                        vz  = new double[sys.Nvz];
                        NTz  = new int[sys.Nvz]; 

                        dt = _num_sys.CFL*dz/v2;
                        for (int j = 0; j < sys.Nz; ++j) 
                            z[j] = sys.z1 + (j+0.5)*dz;
                        for (int k = 0; k < sys.Nvz; ++k){ 
                            vz[k] = v1 + (k+0.5)*dvz;
                            NTz[k] = (int)(ceil((sys.z2 - sys.z1)/(dt*fabs(vz[k]))));
                        }

                        break;
                    case 2:
                        SIZE = (size_t) (sys.Nx+2*gx)*(sys.Nz+2*gz)*sys.Nphi*sys.Nvz;
                        dx  = (sys.x2 - sys.x1)/sys.Nx;
                        dz  = (sys.z2 - sys.z1)/sys.Nz;
                        dvz = (v2 - v1)/sys.Nvz;
                        dphi= (2*M_PI)/sys.Nphi; 
                        x   = new double[sys.Nx];
                        z   = new double[sys.Nz];
                        phi = new double[sys.Nphi];
                    
                        vz  = new double[sys.Nvz];
                        vx  = new double[sys.Nvz*sys.Nphi];
                        vy  = new double[sys.Nvz*sys.Nphi];
                        NTx  = new int[sys.Nvz*sys.Nphi];
                        //NTz  = new int[sys.Nvz];

                        dt = _num_sys.CFL*dz/v2;
                        for (int j = 0; j < sys.Nx; ++j)
                            x[j] = sys.x1 + (j+0.5)*dx;
                        for (int j = 0; j < sys.Nz; ++j)
                            z[j] = sys.z1 + (j+0.5)*dz;
                        for (int kz = 0; kz < sys.Nvz; ++kz){
                            vz[kz] = v1 + (kz+0.5)*dvz;
                            NTz[kz] = (int)(ceil((sys.z2-sys.z1)/(vz[kz]*dt)));
                            for (int f = 0; f < sys.Nphi; ++f){
                                phi[f] = (f + 0.5)*dphi;
                                vx[kz*sys.Nphi+f] = sqrt(1-vz[kz]*vz[kz])*cos(phi[f]);
                                vy[kz*sys.Nphi+f] = sqrt(1-vz[kz]*vz[kz])*sin(phi[f]); 
                                NTx[kz*sys.Nphi + f] = (int)(ceil((sys.x2-sys.x1)/(vx[kz*sys.Nphi + f]*dt)));
                            }
                        }
                        //for (int kz = 0; kz < sys.Nvz; ++kz){
                        //    for (int kx = 0; kx < sys.Nvx; ++kx){
                        //        double v = sqrt(pow(vx[kx],2) + pow(vz[kz],2));
                        //        NT[kz*sys.Nvx+kx] = (int) (ceil(L/(v*dt)));
                        //    }
                        //}
                        break;
                    case 3:
                        //SIZE = (size_t) (sys.Nx+2*gx)*(sys.Ny+2*gy)*(sys.Nz+2*gz)*sys.Nphi*sys.Nvy*sys.Nvz;
                        break;
                }
                rho_exact   = make_field2(SIZE);
                rho         = make_field2(SIZE);
                rho_next    = make_field2(SIZE);
                P           = make_pol2(SIZE);
                //rho_exact.v     = new op2_t[SIZE];
                //rho_exact.v_bar = new op2_t[SIZE];
                //rho.v           = new op2_t[SIZE];
                //rho.v_bar       = new op2_t[SIZE];
                //rho_next.v      = new op2_t[SIZE];
                //rho_next.v_bar  = new op2_t[SIZE];
                //P.v             = new vec2_t[SIZE];
                //P.v_bar         = new vec2_t[SIZE];
            }            
            ~Qke2() {
                switch (sys.dim) {
                    case 1:
                        delete[] z; delete[] vz;
                        delete[] NTz;
                        break;
                    case 2:
                        delete[] x; delete[] z;
                        delete[] vz; delete[] phi; 
                        delete[] vx; delete[] vy;
                        delete[] NTx; //delete[] NTz;
                        break;
                    case 3:
                        break;
                
                } 
                free_field2(rho_exact);
                free_field2(rho);
                free_field2(rho_next);
                free_pol2(P);
                //delete[] rho_exact.v; delete [] rho_exact.v_bar;
                //delete[] rho.v; delete [] rho.v_bar;
                //delete[] rho_next.v; delete[] rho_next.v_bar;
                //delete[] P.v;   delete[] P.v_bar;
            }
            size_t index(int dim, ...) {
                //dim = sys.dim;
                int jx, jy, jz, f, kz;
                va_list valist;
                va_start(valist, dim); 
                if (dim == 1){
                    jz = va_arg(valist, int);
                    kz = va_arg(valist, int);
                    va_end(valist);
                    return (size_t)(kz*(sys.Nz+2*gz) + (jz+gz));
                }
                else if (dim == 2){   
                    jx = va_arg(valist, int);
                    jz = va_arg(valist, int);
                    f = va_arg(valist, int);
                    kz = va_arg(valist, int);
                    va_end(valist);
                    return (size_t)(kz*sys.Nphi*(sys.Nz+2*gz)*(sys.Nx+2*gx) + f*(sys.Nz+2*gz)*(sys.Nx+2*gx) + (jz+gz)*(sys.Nx+2*gx) + (jx+gx));
                } 
                else if (dim == 3){    
                    jx = va_arg(valist, int);
                    jy = va_arg(valist, int);
                    jz = va_arg(valist, int);
                    f = va_arg(valist, int);
                    //ky = va_arg(valist, int);
                    kz = va_arg(valist, int);
                    va_end(valist);
                    return 0;
                }
                else 
                   return 0;
            } 
            //void ccommutator(const cnum a, const op2_t *A, const op2_t *B, op2_t *C);

            void rk4(field2_t *curr, field2_t *next);
            
            void fieldadd(field2_t *C, const cnum a, field2_t *A, const cnum b, field2_t *B);
            void periodic_ghost_zone(field2_t *A);
            
            void injetopen_ghost_zone(field2_t *A);
            
            void compute_F(field2_t *f, const field2_t *A);
            
            void set_init(field2_t *start);

            void field_to_pol(pol2_t *P, const field2_t *fd);
            
            double Gaussian(int dim, ...);

            void vac_exact();

            void output_rho(const char *filename, field2_t *rho);

            void output_P(const char *filename, const pol2_t *P);
            void run();
    };
    double g(double u, double xi, double v);
    
    class Qke3;

}

#define DECLARE_FUNC(flavor)    \
    static inline void ADV##flavor(op##flavor##_t *y, cnum a, const op##flavor##_t *x, size_t N);                   \
    static inline void KO##flavor(op##flavor##_t *y, cnum a, const op##flavor##_t *x, size_t N);                    \
    static inline void MAT##flavor##_ADD(op##flavor##_t *C, cnum a, op##flavor##_t *A, cnum b, op##flavor##_t *B);  \
   static inline void INT##flavor(op##flavor##_t *A, cnum a, const op##flavor##_t *B, const op##flavor##_t *C);     \
   static inline void ASSIGN##flavor(op##flavor##_t *y, op##flavor##_t *x);


DECLARE_FUNC(2)
DECLARE_FUNC(3)
#undef DECLARE_FUNC

#endif
