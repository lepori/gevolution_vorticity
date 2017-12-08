#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

/** 
    differential equation for the growth rate D_1 
**/

double g1(double a, void *p)
{
  struct cosmology par=*(struct cosmology *) p;
  double g1=1./(pow(a,3)*pow(1 - par.Omega_m + par.Omega_m/pow(a,3),1.5));
  return g1;
}


int growth_rate_ode(
    double t, const double y[], 
    double f[], void *params
)
{
    (void)t;
    struct cosmology par = *(struct cosmology *) params;

    f[0] = y[1];
    f[1] = -3./2*(1+((-1+par.Omega_m)*((-1.)))/(1+(-1+pow(t,-3))*par.Omega_m))*y[1]/t
         +3./2*(pow(t,3*((-1.)))*par.Omega_m)/(1+(-1+pow(t,3*((-1.))))*par.Omega_m)*y[0]/pow(t,2);
    return GSL_SUCCESS;
}


/** 
    jacobian for the growth rate differential equation 
**/

int growth_rate_jac(
    double t, const double y[], double *dfdy,
    double dfdt[], void *params
)
{
    struct cosmology par = *(struct cosmology *) params;
    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, 3./2*(pow(t,3*(-1.))*par.Omega_m)/(1+(-1+pow(t,3*(-1.)))*par.Omega_m)/pow(t,2));
    gsl_matrix_set(m, 1, 1, -3./2*(1+((-1+par.Omega_m)*(-1.))/(1+(-1+pow(t,3*(-1.)))*par.Omega_m))/t);
    dfdt[0] = 0.0;
    dfdt[1] = (-3*t*pow(-1 + par.Omega_m,2)*(-1 + (-1.))*y[1] 
            + 3*pow(t,6*(-1.))*pow(par.Omega_m,2)*(-2*y[0] + t*y[1]) 
            + 3*pow(t,3*(-1.))*(-1 + par.Omega_m)*par.Omega_m*(-2 + 3*(-1.))
            *(-y[0] + t*(1 + (-1.))*y[1]))/
           (2.*pow(t,3)*pow(1 + (-1 + pow(t,3*(-1.)))*par.Omega_m,2));
    return GSL_SUCCESS;
}


/** 
    computes and stores all the background functions 
**/

double D1_prime(
    struct cosmology par, double a
)
{

  /*
    double temp_a = 0.005, temp_z;
    gsl_odeiv_system sys
        = {growth_rate_ode, growth_rate_jac, 2, par};

    const gsl_odeiv_step_type *step_type
        = gsl_odeiv_step_rk8pd;

    gsl_odeiv_step *step
        = gsl_odeiv_step_alloc(step_type, 2);
    gsl_odeiv_control *control
        = gsl_odeiv_control_y_new(1e-6, 0.0);
    gsl_odeiv_evolve *evolve
        = gsl_odeiv_evolve_alloc(2);

  
    double initial_values[2] = {0.005, 1.0};

    double h = 1E-6, prec = 1E-5;
#if 0
    for (size_t i = 0; i<MAX_BINS; ++i){
        temp_a = 0.05;
        initial_values[0] = 0.05;
        initial_values[1] = 1.0;

        temp_z = (par->z_mean + par->deltaz + par->deltaz/100.)*i/(double)(MAX_BINS);

        (temp_bg->z)[i] = temp_z;
        (temp_bg->a)[i] = 1./(1.+temp_z);
        (temp_bg->Hz)[i] = sqrt(
             par->Omega0_m*pow(1+temp_z, 3)
            +par->Omega0_gamma*pow(1+temp_z, 4)
            +par->Omega0_de*pow(1+temp_z, 3*(par->w+1))
        );
        (temp_bg->conformal_Hz)[i] = (temp_bg->a)[i]*(temp_bg->Hz)[i];
#endif
        while (temp_a < a){
            gsl_odeiv_evolve_apply(
                evolve, control, step,
                &sys, &temp_a, a,
                &h, initial_values
            );
        }

*/

    gsl_function G1 = {.function = &g1, .params = &par };
    double a_minus = a-0.005*a;
    double a_plus  = a+0.005*a;
    gsl_integration_workspace *w;
    double prec=1E-5;
    double result, error;
    w=gsl_integration_workspace_alloc(10000);
    gsl_integration_qag (&G1, 0, a_plus, 0, prec, 10000,
		         GSL_INTEG_GAUSS61, w, &result, &error);

    double result_plus=(5./2*par.Omega_m*sqrt(1 - par.Omega_m + par.Omega_m/pow(a_plus,3)))*result; 

    gsl_integration_qag (&G1, 0, a_minus, 0, prec, 10000,
                         GSL_INTEG_GAUSS61, w, &result, &error);

    double result_minus=(5./2*par.Omega_m*sqrt(1 - par.Omega_m + par.Omega_m/pow(a_minus,3)))*result;
 
    double result_der = (result_plus -  result_minus)/(0.01*a);

    result = result_der*Hconf(a, 1, par)*a;

    //gsl_odeiv_step_free(step);
    //gsl_odeiv_control_free(control);
    //gsl_odeiv_evolve_free(evolve);

    return result;
}
