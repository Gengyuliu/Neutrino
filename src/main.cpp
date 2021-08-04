#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "neu.hpp"

/*	driver program	*/
int main(int argc, char *argv[])
{
    phys_sys_t osc;
    num_sys_t sys; 
    bzero(&osc, sizeof(phys_sys_t));
    bzero(&sys, sizeof(num_sys_t));
    Neu::parsing(argc, argv, &osc, &sys); 
    Neu::Qke2 state(osc, sys);
    state.run(); 
    return 0;
}
