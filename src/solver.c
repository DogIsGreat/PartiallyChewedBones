#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

double uE(double x, double t, double c){
    // Asin((pi*x)/L)cos((pi*c*t)/L)
    return M_PI;
}

double c1(double c, double dt, double dx){
    return c * (dt / dx);
}

double c2(double c){
    return gsl_pow_2(c);
}


int main()
{
// Given Mesh points as arrays x and t (x[i], t[n])
    int k;
    int n;
    int i;
    int c;
    int y;

    int* x = (int*)malloc(k * sizeof(int));
    double* t = (double*)malloc(n * sizeof(double));
    double* u = (double*)malloc(i * sizeof(double));
    double* u_n = (double*)malloc(i * sizeof(double));
    double* u_nm = (double*)malloc(i * sizeof(double));
    double* u_nm1 = (double*)malloc(i * sizeof(double));
    double* Iarray =(double*)malloc(i * sizeof(double));

    double dx;
    double dy;
    double dt;

    int nt = sizeof(t)/sizeof(double) -1;
    int nx = sizeof(x)/sizeof(double) -1;

    // Set initial condition u(x,0) = I(x)
    for ( int i = 0; i <= nx; ++i ){
        u_n[i] = Iarray[x[i]];
    } 

    // Apply special formula for first step, incorporating du/dt =0
    for ( int i = 0; i < nx; ++i ){
        u[i] = u_n[i] - 0.5*c1(c,dt,dx)*(u_n[i+1]-2*u_n[i]+u_n[i-1]);
    }

    // Enforce boundary conditions
    u[0] = 0;
    u[nx] = 0;

    // Switch Variables before next step
    for (int i = 0; i <= nx; ++i){
        u_nm1[i] = u_n[i];
    }
    for (int i = 0; i <= nx; ++i){
        u_n[i] = u[i];
    }

    for( int i = 0; i <= nt; ++i){
        // Update all inner mesh points at time t[n+1]
        for (int j = 0; j <= nx; ++j){
            u[j] = 2*u_n[j]-u_nm1[j] - c1(c,dt,dx)*(u_n[j+1] - 2*u_n[j]+ u_n[j-1]);
        }

        // Insert boundary conditions
        u[0] = 0;
        u[nx] = 0;

        // Switch Variables before next step
        for (int i = 0; i <= nx; ++i){
            u_nm1[i] = u_n[i];
        }
        for (int i = 0; i <= nx; ++i){
            u_n[i] = u[i];
        }

    }
/*
    if(!y){
        printf ("error: %s\n", gsl_strerror (y));
    }
*/
    printf ("J0(%g) = %.18e\n", x, y);

    free(x);
    free(y);
    free(u);
    free(u_n);
    free(u_nm);
    free(u_nm1);
    free(Iarray);
    return 0;

}

  
