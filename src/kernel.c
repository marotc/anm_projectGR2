#include "kernel.h"

double eval_Cubic_kernel(double q, double h) {
    double alpha = 15.0/(7.0*M_PI*h*h);
    if(q >= 0. && q <= 1.) return alpha * (2./3. - q*q + q*q*q/2.);
    else if(q > 1. && q <= 2.) return alpha * ((2-q)*(2-q)*(2-q)/6.);
    else return 0.;
}

double eval_kernel(xy *p1, xy *p2, double kh, Kernel kernel) {
    double d = sqrt(squared(p1->x-p2->x) + squared(p1->y-p2->y));
    double h;
    if(kernel == Cubic) {
        h = kh/2.;
        return eval_Cubic_kernel(d/h, h);
    }
    else {
        fprintf(stderr, "Other kernel functions cannot be used.\n");
        exit(-1);
    }
}

double derivative_Cubic_kernel(double q,double h)
{
     double alpha = 15.0/(7.0*M_PI*pow(h, 2));
     double g;
     if (q >= 0 && q <= 1)
         g = -2 * q + 3./2.*q*q;
     else if (q > 1 && q <= 2)
         g = -0.5*(2 - q)*(2 - q);
     else
         g = 0;
     return alpha*g;
 }


xy* grad_kernel(xy* p1, xy* p2, double kh, Kernel kernel) {
    if(p1->x == p2->x && p1->y == p2->y) return xy_new(0,0);
    double d_x = p1->x-p2->x;
    double d_y = p1->y-p2->y;
	double d = sqrt(d_x*d_x + d_y*d_y);
	double g;
	double h;
	if (kernel == Cubic) {
		h = kh / 2;
		g = derivative_Cubic_kernel(d / h, h);
    }
    else {
        fprintf(stderr, "Other kernel functions cannot be used.\n");
        exit(-1);
    }
    
	double g_x = g*(d_x / (h*d));
	double g_y = g*(d_y / (h*d));

	return xy_new(g_x, g_y);
}

double deriv2_Cubic_kernel(xy* p1, xy* p2, double kh, Kernel kernel) {
    if(p1->x == p2->x && p1->y == p2->y) return 0.;
    double h = kh / 2;
    double d_x = p1->x-p2->x;
    double d_y = p1->y-p2->y;
    double d = sqrt(d_x*d_x + d_y*d_y);
    double q = d/h;
    double g;
    double alpha = 15.0/(7.0*M_PI*pow(h, 2));
    if (q >= 0 && q <= 1)
        g = -2. + 3.*q*q;
    else if (q > 1 && q <= 2)
        g = (2.-q);
    else
        g = 0;
    return alpha*g;
}
