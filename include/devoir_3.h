#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct{
    double *u; // displacement (de taille 2*n)
    double *v; // velocity (de taille 2*n)
} State;

State* newMark(State *state_0, double T, double dt, double beta, double gamma, CSRMatrix *Ksp, CSRMatrix *Msp,int I);