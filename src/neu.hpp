#ifndef NEU_H
#define NEU_H
#include <stdio.h>
#include <stdarg.h>
#include <complex.h>
#define cnum _Complex double

typedef struct phys_sys_t {
    int     flavor;
    double  pmo = 1;
    double  theta = 33.3;      /*	mixing angle	*/
    double  mu = 1.0;         /*	interaction strength	*/
 
    double  t; 
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
    
    int     Nvx, Nvy, Nvz;
} num_sys_t;

/*	2-flavor
 *	SU(2)	*/
typedef struct operator2_t {
    cnum _11; cnum _12;
    cnum _21; cnum _22;
} op2_t;
typedef struct vector2_t {
    double _1_;
    double _2_;
    double _3_;
} vec2_t;

typedef struct field2_t {
    op2_t *v;
    op2_t *v_bar;
} field2_t;

typedef struct pair2_t{
    op2_t v;
    op2_t v_bar;
} pair2_t;

typedef struct polarization_t {
    vec2_t *v;
    vec2_t *v_bar;
} pol_t;

/*	3-flavor
 *	SU(3)	*/
typedef struct operator3_t {
    cnum _11; cnum _12; cnum _13;
    cnum _21; cnum _22; cnum _23;
    cnum _31; cnum _32; cnum _33;
} op3_t;
typedef struct vector3_t {
    double _1_;
    double _2_; 
    double _3_; 
    double _4_; 
    double _5_; 
    double _6_; 
    double _7_; 
    double _8_; 
} vec3_t;

typedef struct field3_t {
    op3_t *v;
    op3_t *v_bar;
} field3_t;

typedef struct pair3_t {
    op3_t v;
    op3_t v_bar;
} pair3_t;


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
            field2_t     rho;      
            field2_t     rho_next;
            pol_t        P; 
            double  dx, dy, dz;
            double  dvx, dvy, dvz;
            double  *x, *y, *z;
            double  *vx, *vy, *vz;
            int     gx, gy, gz;
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
                        
                        dt = _num_sys.CFL*dz/v2;
                        for (int j = 0; j < sys.Nz; ++j) 
                           z[j] = sys.z1 + (j+0.5)*dz;
                        for (int k = 0; k < sys.Nvz; ++k) 
                           vz[k] = v1 + (k+0.5)*dvz;

                        break;
                    case 2:
                        SIZE = (size_t) (sys.Nx+2*gx)*(sys.Nz+2*gz)*sys.Nvx*sys.Nvz;
                        dx  = (sys.x2 - sys.x1)/sys.Nx;
                        dz  = (sys.z2 - sys.z1)/sys.Nz;
                        dvx = (v2 - v1)/sys.Nvx;
                        dvz = (v2 - v1)/sys.Nvz;
                        x   = new double[sys.Nx];
                        z   = new double[sys.Nz];
                        vx  = new double[sys.Nvx];
                        vz  = new double[sys.Nvz];

                        dt = _num_sys.CFL*dz/v2;
                        for (int j = 0; j < sys.Nx; ++j)
                            x[j] = sys.x1 + (j+0.5)*dx;
                        for (int j = 0; j < sys.Nz; ++j)
                            z[j] = sys.z1 + (j+0.5)*dz;
                        for (int k = 0; k < sys.Nvz; ++k)
                            vz[k] = v1 + (k+0.5)*dvz;
                        for (int k = 0; k < sys.Nvx; ++k)
                            vx[k] = v1 + (k+0.5)*dvx;
                        break;
                    case 3:
                        SIZE = (size_t) (sys.Nx+2*gx)*(sys.Ny+2*gy)*(sys.Nz+2*gz)*sys.Nvx*sys.Nvy*sys.Nvz;
                        break;
                }
                rho.v           = new op2_t[SIZE];
                rho.v_bar       = new op2_t[SIZE];
                rho_next.v      = new op2_t[SIZE];
                rho_next.v_bar  = new op2_t[SIZE];
                P.v             = new vec2_t[SIZE];
                P.v_bar         = new vec2_t[SIZE];
            }            
            ~Qke2() {
                switch (sys.dim) {
                    case 1:
                        delete[] z; delete[] vz;
                        break;
                    case 2:
                        delete[] x; delete[] z;
                        delete[] vx; delete[] vz;
                        break;
                    case 3:
                        break;
                
                } 
                delete[] rho.v; delete [] rho.v_bar;
                delete[] rho_next.v; delete[] rho_next.v_bar;
                delete[] P.v;   delete[] P.v_bar;
            }
            size_t index(int dim, ...) {
                //dim = sys.dim;
                int jx, jy, jz, kx, ky, kz;
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
                    kx = va_arg(valist, int);
                    kz = va_arg(valist, int);
                    va_end(valist);
                    return (size_t)(kz*sys.Nvx*(sys.Nz+2*gz)*(sys.Nx+2*gx) + kx*(sys.Nz+2*gz)*(sys.Nx+2*gx) + (jz+gz)*(sys.Nx+2*gx) + (jx+gx));
                } 
                else if (dim == 3){    
                    jx = va_arg(valist, int);
                    jy = va_arg(valist, int);
                    jz = va_arg(valist, int);
                    kx = va_arg(valist, int);
                    ky = va_arg(valist, int);
                    kz = va_arg(valist, int);
                    va_end(valist);
                    return (size_t)(kz*sys.Nvx*sys.Nvy*(sys.Nz+2*gz)*(sys.Ny+2*gy)*(sys.Nx+2*gx) + ky*sys.Nvx*(sys.Nz+2*gz)*(sys.Ny+2*gy)*(sys.Nx+2*gx) + kx*(sys.Nz+2*gz)*(sys.Ny+2*gy)*(sys.Nx+2*gx) + (jz+gz)*(sys.Nx+2*gx)*(sys.Ny+2*gy) + (jy+gy)*(sys.Nx+2*gx) + (jx+gx));
                }
                else 
                   return 0;
            } 
            void ccommutator(const cnum a, const op2_t *A, const op2_t *B, op2_t *C);

            void rk4(field2_t *curr, field2_t *next);
            
            void fieldadd(field2_t *C, const cnum a, const field2_t *A, const cnum b, const field2_t *B);
            void periodic_ghost_zone(field2_t *A);
            
            void injetopen_ghost_zone(field2_t *A);
            
            void compute_F(field2_t *f, const field2_t *A);
            
            void set_init(field2_t *start);

            void field_to_pol(pol_t *P, const field2_t *fd);
            
            double Gaussian(int dim, ...);

            void output_rho(const char *filename, field2_t *rho);

            void output_P(const char *filename, pol_t *P);
            void run();
    };
    double g(double u, double xi);
    
    class Qke3;
    


}
#endif
