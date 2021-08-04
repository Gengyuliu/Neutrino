#ifndef OBJECTS_H
#define OBJECTS_H


typedef _Complex double cnum;

//void _make_Matrix(Matrix *m, int size);

//void _free_Matrix(Matrix *m);

/*	object type depends on the neutrino flavor system 	*/
#define DECLARE_OBJECTS(flavor)                             \
    typedef struct {                                        \
        cnum element[flavor*flavor];                        \
    } op##flavor##_t;                                       \
    typedef struct {                                        \
        int dim;                                            \
        double *data;                                       \
    } vec##flavor##_t;                                      \
    void _make_vec##flavor(vec##flavor##_t *v, int flvr);   \
    void _free_vec##flavor(vec##flavor##_t *v);             \
                                                            \
    typedef struct {                                        \
        int size;                                           \
        op##flavor##_t *v;                                  \
        op##flavor##_t *v_bar;                              \
    } field##flavor##_t;                                    \
    field##flavor##_t* make_field##flavor(int size);        \
    void free_field##flavor(field##flavor##_t *fd);         \
                                                            \
    typedef struct {                                        \
        op##flavor##_t v;                                   \
        op##flavor##_t v_bar;                               \
    } pair##flavor##_t;                                     \
                                                            \
    /*void init_pair##flavor(pair##flavor##_t *pr);           \
    void end_pair##flavor(pair##flavor##_t *pr);            */\
    void zero_pair##flavor(pair##flavor##_t *pr);           \
                                                            \
    typedef struct {                                        \
        int size;                                           \
        vec##flavor##_t *v;                                 \
        vec##flavor##_t *v_bar;                             \
    } pol##flavor##_t;                                      \
    pol##flavor##_t* make_pol##flavor(int size);            \
    void free_pol##flavor(pol##flavor##_t *pol);            

DECLARE_OBJECTS(2);
DECLARE_OBJECTS(3);


#undef DECLARE_OBJECTS

#define FOR(m, M) \
    for (int m = 0; m < M; ++m) 

#endif
