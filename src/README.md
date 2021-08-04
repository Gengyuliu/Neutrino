# Neutrino Oscillation Simulation

## Objects in Simulation

We declare required objects to suit the need of different flavor system. 

For example, `DECLARE_OBJECT(2)` means that we declare objects for two-flavor system.

The type of density matrix we used is `op##flavor##_t`, which is all defined inside the macro `DECLARE_OBJECT`:
```c=
#define DECLARE_OBJECTS(flavor)                 \
...                                             \
   typedef struct {                             \            
           cnum element[flavor*flavor];         \                
   } op##flavor##_t;                            \            
   typedef struct {                             \            
           int dim;                             \               
           double *data;                        \                
   } vec##flavor##_t;                           \
   typedef struct {                             \            
           int size;                            \                
           op##flavor##_t *v;                   \                
           op##flavor##_t *v_bar;               \              
   } field##flavor##_t;                         \
```

## Algorithm

We used Runge-Kutta methods along with 4-th order finite difference in space discretization. The detailed information is included in `doc`.
