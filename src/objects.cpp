#include <stdio.h>
#include <stdlib.h>

#include "./objects.hpp"


#define DECLARE_FUNC(flavor)                                                    \
    void _make_vec##flavor(vec##flavor##_t *v, int flvr){                       \
        v->dim = flvr*flvr - 1;                                                 \
        v->data = (double*)malloc(sizeof(double)*v->dim);                       \
    }                                                                           \
                                                                                \
    void _free_vec##flavor(vec##flavor##_t *v){                                 \
        free(v->data);                                                          \
    }                                                                           \
                                                                                \
    field##flavor##_t* make_field##flavor(int size) {                           \
        field##flavor##_t *fd = (field##flavor##_t*)malloc(sizeof(field##flavor##_t));\
        fd->size = size;                                                        \
        fd->v = (op##flavor##_t*)malloc(sizeof(op##flavor##_t)*size);           \
        fd->v_bar = (op##flavor##_t*)malloc(sizeof(op##flavor##_t)*size);       \
        /*for (int i = 0; i < size; ++i){                                         \
            _make_Matrix(&fd->v[i], flavor*flavor);                             \
            _make_Matrix(&fd->v_bar[i], flavor*flavor);                         \
        }*/\
       return fd;                                                              \
    }                                                                           \
                                                                                \
    void free_field##flavor(field##flavor##_t *fd){                             \
        free(fd->v);                                                            \
        free(fd->v_bar);                                                        \
        free(fd);                                                               \
    }                                                                           \
                                                                                \
    /*void init_pair##flavor(pair##flavor##_t *pr){                               \
        _make_Matrix(&pr->v, flavor*flavor);                                    \
        _make_Matrix(&pr->v_bar, flavor*flavor);                                \
    }                                                                           \
    void end_pair##flavor(pair##flavor##_t *pr){                                \
        _free_Matrix(&pr->v);                                                   \
        _free_Matrix(&pr->v_bar);                                               \
    }*/                                                                         \
    void zero_pair##flavor(pair##flavor##_t *pr){                               \
        FOR(m, flavor*flavor){                                                  \
            pr->v.element[m] = 0.0;                                             \
            pr->v_bar.element[m] = 0.0;                                         \
        }                                                                       \
    }                                                                           \
    pol##flavor##_t* make_pol##flavor(int size){                                \
        pol##flavor##_t *pol = (pol##flavor##_t*)malloc(sizeof(pol##flavor##_t));   \
        pol->size = size;                                                       \
        pol->v = (vec##flavor##_t*)malloc(sizeof(vec##flavor##_t)*size);        \
        pol->v_bar = (vec##flavor##_t*)malloc(sizeof(vec##flavor##_t)*size);    \
        for (int i = 0; i < size; ++i){                                         \
            _make_vec##flavor(&pol->v[i], flavor);                              \
            _make_vec##flavor(&pol->v_bar[i], flavor);                          \
        }                                                                       \
        return pol;                                                             \
    }                                                                           \
                                                                                \
    void free_pol##flavor(pol##flavor##_t *pol){                                \
        for (int i = 0; i < pol->size; ++i){                                    \
            _free_vec##flavor(&pol->v[i]);                                      \
            _free_vec##flavor(&pol->v_bar[i]);                                  \
        }                                                                       \
        free(pol->v);  free(pol->v_bar);                                        \
        free(pol);                                                              \
    }                                                                            


DECLARE_FUNC(2)
DECLARE_FUNC(3)
    
#undef DECLARE_FUNC

