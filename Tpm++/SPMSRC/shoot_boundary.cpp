/* $Id: shoot_boundary.c 4875 2012-08-30 20:04:30Z john $ */
/* (c) John Ashburner (2011) */

#include "StdAfx.h"
#include "shoot_boundary.h"

/* Neumann boundary condition */
static int neumann_boundary(int i, unsigned int m)
{
	int m2 = m*2;
    i = (i<0) ? m2-((-i-1)%m2)-1 : (i%m2);
    return((m<=i)? m2-i-1: i);
}

static int circulant_boundary(int i, unsigned int m)
{
    return((i>=0) ? i%((signed)m) : (((signed)m)+i%((signed)m))%(signed)m);
}

int (*bound)(int, unsigned int) = circulant_boundary;
static int bound_type = BOUND_CIRCULANT;

void set_bound(int t)
{
    bound_type = t;
    if (t==BOUND_CIRCULANT)
        bound = circulant_boundary;
    else if (t==BOUND_NEUMANN)
        bound = neumann_boundary;
    else
        ErrOutput("set_bound","Undefined boundary condition.");
}

int get_bound()
{
    return(bound_type);
}

