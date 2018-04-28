/////////////////////////
// vorticity.hpp
//////////////////////////
// 
// Function use to extract the vorticity field from gevolution
// 
// Last modified: December 2017
//
//////////////////////////

#ifndef VORTICITY_HEADER
#define VORTICITY_HEADER

#include "prng_engine.hpp"
#include "d1_prime.hpp"

#include <gsl/gsl_spline.h>

using namespace std;
using namespace LATfield2;

#ifndef MAX_LINESIZE
#define MAX_LINESIZE 2048
#endif

#ifndef Cplx
#define Cplx Imag
#endif

using namespace std;
using namespace LATfield2;

// should be larger than maximum Ngrid                                                                                    
#ifndef HUGE_SKIP
#define HUGE_SKIP   65536
#endif
/////////////////////////////
// loadTransferFunctions_vel
/////////////////////////////
// Description:                                                                                                             
//   load the transfer function from CLASS for psi and theta
//  
//                                                                                                                           
// Arguments:                                                                                                              
//   tk_psi       reference to the transfer function tk_psi
//   tk_theta     reference to the transfer function tk_theta
//   qname        string of the species (e.g. 'cdm' or 'b') 
//   boxsize      boxsize in physical units (as in the setting files)
//   h            h from cosmology structure
// Returns:                                                                                                                       
//                                                                                                                                 //////////////////////////                                                                                                               
void loadTransferFunctionsVel(const char * filename, gsl_spline * & tk_psi, gsl_spline * & tk_theta, const char * qname, const double boxsize, const double h)
{
  int i = 0, numpoints = 0;
  double * k;
  double * tk_p;
  double * tk_t;

  if (parallel.grid_rank()[0] == 0) // read file
    {
      FILE * tkfile;
      char line[MAX_LINESIZE];
      char format[MAX_LINESIZE];
      char * ptr;
      double dummy[3];
      int kcol = -1, dcol = -1, tcol = -1, colmax;
      
      line[MAX_LINESIZE-1] = 0;
      
      tkfile = fopen(filename, "r");
      
      if (tkfile == NULL)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Unable to open file " << filename << "." << endl;
	  parallel.abortForce();
	}
      
      while (!feof(tkfile) && !ferror(tkfile))
	{
	  fgets(line, MAX_LINESIZE, tkfile);
	  if (line[MAX_LINESIZE-1] != 0)
	    {
	      cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Character limit (" << (MAX_LINESIZE-1) << "/line) exceeded in file " << filename << "." << endl;
	      fclose(tkfile);
	      parallel.abortForce();
	    }
	  
	  if (line[0] != '#' && !feof(tkfile) && !ferror(tkfile)) numpoints++;
	}
      
      if (numpoints < 2)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! No valid data found in file " << filename << "." << endl;
	  fclose(tkfile);
	  parallel.abortForce();
	}
      
      k = (double *) malloc(sizeof(double) * numpoints);
      tk_p = (double *) malloc(sizeof(double) * numpoints);
      tk_t = (double *) malloc(sizeof(double) * numpoints);
      
      if (k == NULL || tk_p == NULL || tk_t == NULL)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Memory error." << endl;
	  fclose(tkfile);
	  parallel.abortForce();
	}
      
      rewind(tkfile);
      
      while (!feof(tkfile) && !ferror(tkfile))
	{
	  fgets(line, MAX_LINESIZE, tkfile);
	  for (ptr = line, i = 0; (ptr = strchr(ptr, ':')) != NULL; i++)
	    {
	      ptr++;
	      if (*ptr == 'k') kcol = i;
	      else if (*ptr == 'p')
		{
		  if (strncmp(ptr, "psi", strlen("psi")) == 0) dcol = i;
		}
	      else if (*ptr == 't')
		{
		  if (strncmp(ptr+2, qname, strlen(qname)) == 0) tcol = i;
		}
	    }
	  
	  if (kcol >= 0 && dcol >= 0 && tcol >= 0) break;
	}
      
      if (kcol < 0 || dcol < 0 || tcol < 0)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Unable to identify requested columns!" << endl;
	  fclose(tkfile);
	  free(k);
	  free(tk_p);
	  free(tk_t);
	  parallel.abortForce();
	}
      
      colmax = i;
      for (i = 0, ptr=format; i < colmax; i++)
	{
	  if (i == kcol || i == dcol || i == tcol)
	    {
	      strncpy(ptr, " %lf", 4);
	      ptr += 4;
	    }
	  else
	    {
	      strncpy(ptr, " %*lf", 5);
	      ptr += 5;
	    }
	}
      *ptr = '\0';
      
      if (kcol < dcol && dcol < tcol)
	{
	  kcol = 0; dcol = 1; tcol = 2;
	}
      else if (kcol < tcol && tcol < dcol)
	{
	  kcol = 0; dcol = 2; tcol = 1;
	}
      else if (dcol < kcol && kcol < tcol)
	{
	  kcol = 1; dcol = 0; tcol = 2;
	}
      else if (dcol < tcol && tcol < kcol)
	{
	  kcol = 2; dcol = 0; tcol = 1;
	}
      else if (tcol < kcol && kcol < dcol)
	{
	  kcol = 1; dcol = 2; tcol = 0;
	}
      else if (tcol < dcol && dcol < kcol)
	{
	  kcol = 2; dcol = 1; tcol = 0;
	}
      else
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Inconsistent columns!" << endl;
	  fclose(tkfile);
	  free(k);
	  free(tk_p);
	  free(tk_t);
	  parallel.abortForce();
	}
      
      i = 0;
      while (!feof(tkfile) && !ferror(tkfile))
	{
	  fgets(line, MAX_LINESIZE, tkfile);
	  
	  if (sscanf(line, format, dummy, dummy+1, dummy+2) == 3 && !feof(tkfile) && !ferror(tkfile))
	    {
	      if (dummy[kcol] < 0.)
		{
		  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Negative k-value encountered." << endl;
		  free(k);
		  free(tk_p);
		  free(tk_t);
		  fclose(tkfile);
		  parallel.abortForce();
		}
	      
	      if (i > 0)
		{
		  if (k[i-1] >= dummy[kcol] * boxsize)
		    {
		      cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! k-values are not strictly ordered." << endl;
		      free(k);
		      free(tk_p);
		      free(tk_t);
		      fclose(tkfile);
		      parallel.abortForce();
		    }
		}
	      
	      k[i] = dummy[kcol] * boxsize;
	      tk_p[i] = dummy[dcol];
	      tk_t[i] = dummy[tcol] * boxsize / h;
	      i++;
	    }
	}
      
      fclose(tkfile);
      
      if (i != numpoints)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! File may have changed or file pointer corrupted." << endl;
	  free(k);
	  free(tk_p);
	  free(tk_t);
	  parallel.abortForce();
	}
      
      parallel.broadcast_dim0<int>(numpoints, 0);
    }
  else
    {
      parallel.broadcast_dim0<int>(numpoints, 0);
      
      if (numpoints < 2)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Communication error." << endl;
	  parallel.abortForce();
	}
      
      k = (double *) malloc(sizeof(double) * numpoints);
      tk_p = (double *) malloc(sizeof(double) * numpoints);
      tk_t = (double *) malloc(sizeof(double) * numpoints);
      
      if (k == NULL || tk_p == NULL || tk_t == NULL)
	{
	  cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Memory error." << endl;
	  parallel.abortForce();
	}
    }
  
  parallel.broadcast_dim0<double>(k, numpoints, 0);
  parallel.broadcast_dim0<double>(tk_p, numpoints, 0);
  parallel.broadcast_dim0<double>(tk_t, numpoints, 0);
  
  tk_psi = gsl_spline_alloc(gsl_interp_cspline, numpoints);
  tk_theta = gsl_spline_alloc(gsl_interp_cspline, numpoints);
  
  gsl_spline_init(tk_psi, k, tk_p, numpoints);
  gsl_spline_init(tk_theta, k, tk_t, numpoints);
  
  free(k);
  free(tk_p);
  free(tk_t);
}

#ifdef FFT3D

//////////////////////////
// generateDisplacementField for the velocity field (generateRealizationVel)
//////////////////////////
// Description:
//   generates realization from the linear velocity field theta
//   (copied from gevolution.hpp)
//
// Non-type template parameters:
//   ignorekernel  this is effectively an optimization flag defaulted to 0; instantiating with 1 instead will cause
//                 the function to ignore the convolution kernel, allowing the function to be used for generating
//                 realizations (generateRealization is simply an alias for generateDisplacementField<1>)
// 
// Arguments:
//   potFT         reference to allocated field that contains the convolution kernel relating the potential
//                 (generating the displacement field) with the bare density perturbation; will contain the
//                 Fourier image of the potential generating the displacement field
//   coeff         gauge correction coefficient "H_conformal^2"
//   pkspline      pointer to a gsl_spline which holds a tabulated power spectrum
//   seed          initial seed for random number generator
//   ksphere       flag to indicate that only a sphere in k-space should be initialized
//                 (default = 0: full k-space cube is initialized)
//   deconvolve_f  flag to indicate deconvolution function
//                 0: no deconvolution
//                 1: sinc (default)
//
// Returns:
// 
//////////////////////////

#ifndef generateReal
#define generateReal generateRealizationVel<1>
#endif

template<int ignorekernel = 1>
void generateRealizationVel(Field<Cplx> & potFT, const Real coeff, const gsl_spline * pkspline, const unsigned int seed, const int ksphere = 0, const int deconvolve_f = 1)
{
	const int linesize = potFT.lattice().size(1);
	const int kmax = (linesize / 2) - 1;
	rKSite k(potFT.lattice());
	int kx, ky, kz, i, j;
	int kymin, kymax, kzmin, kzmax;
	long jumpy, jumpz;
	float r1, r2, k2, s;
	float * sinc;
	sitmo::prng_engine prng;
	uint64_t huge_skip = HUGE_SKIP;
	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	
	sinc = (float *) malloc(linesize * sizeof(float));
	
	sinc[0] = 1.;
	if (deconvolve_f == 1)
	{
		for (i = 1; i < linesize; i++)
			sinc[i] = sin(M_PI * (float) i / (float) linesize) * (float) linesize / (M_PI * (float) i);
	}
	else
	{
		for (i = 1; i < linesize; i++)
			sinc[i] = 1.;
	}
	
	k.initialize(potFT.lattice(), potFT.lattice().siteLast());
	kymax = k.coord(1);
	kzmax = k.coord(2);
	k.initialize(potFT.lattice(), potFT.lattice().siteFirst());
	kymin = k.coord(1);
	kzmin = k.coord(2);
		
	if (kymin < (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
	{
		prng.seed(seed);
		   
		if (kymin == 0 && kzmin == 0)
		{
			k.setCoord(0, 0, 0);
			potFT(k) = Cplx(0.,0.);
			kx = 1;
		}
		else
		{
			kx = 0;
			prng.discard(((uint64_t) kzmin * huge_skip + (uint64_t) kymin) * huge_skip); 
		}
		
		for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				for (i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) (ky * ky) + (float) (kz * kz);
					
					if (kx >= kmax || ky >= kmax || kz >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[ky] * sinc[kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
				kx = 0;
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymax >= (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + (uint64_t) kzmin) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		
		for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				for (kx = 0, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
										
					k2 = (float) (kx * kx) + (float) ((linesize-ky) * (linesize-ky)) + (float) (kz * kz);
					
					if (kx >= kmax || (linesize-ky) >= kmax || kz >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[linesize-ky] * sinc[kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
						
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymin < (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) kymin) * huge_skip);
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) (ky * ky) + (float) ((linesize-kz) * (linesize-kz));
					
					if (kx >= kmax || ky >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[ky] * sinc[linesize-kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
		
		prng.seed(seed);
		prng.discard(((uint64_t) (linesize - kzmax) * huge_skip + (uint64_t) kymin) * huge_skip);
		kx = 0;
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				k.setCoord(kx, ky, kz);
					
				k2 = (float) (ky * ky) + (float) ((linesize-kz) * (linesize-kz));
				i = 0;
				
				if (ky >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
				{
					potFT(k) = Cplx(0., 0.);
				}
				else
				{
					s = sinc[ky] * sinc[linesize-kz];
					k2 *= 4. * M_PI * M_PI;
					do
					{
						r1 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					}
					while (r1 == 0.);
					r2 = (float) prng() / (float) sitmo::prng_engine::max();
					i++;
				
					potFT(k) = (ignorekernel? Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
				}
								
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymax >= (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) ((linesize-ky) * (linesize-ky)) + (float) ((linesize-kz) * (linesize-kz));
					
					if (kx >= kmax || (linesize-ky) >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[linesize-ky] * sinc[linesize-kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
		
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		kx = 0;
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				k.setCoord(kx, ky, kz);
					
				k2 = (float) ((linesize-ky) * (linesize-ky)) + (float) ((linesize-kz) * (linesize-kz));
				i = 0;
				
				if ((linesize-ky) >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
				{
					potFT(k) = Cplx(0., 0.);
				}
				else
				{
					s = sinc[linesize-ky] * sinc[linesize-kz];
					k2 *= 4. * M_PI * M_PI;
					do
					{
						r1 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					}
					while (r1 == 0.);
					r2 = (float) prng() / (float) sitmo::prng_engine::max();
					i++;
				
					potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
				}
				
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	gsl_interp_accel_free(acc);
	free(sinc);
}

#endif

//  Compute primordial Pk, used for the computation of the linear velocity field
inline double Pk_primordial_vort(const double k, const icsettings & ic)
{
  return ic.A_s * pow(k / ic.k_pivot, ic.n_s - 1.);  // note that k_pivot is in units of inverse Mpc!                         
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Description:                                                                                                                  
//   subtract the linear velocity field, computed from CLASS transfer function                                                        
//                                                                                                                                 
//                                                                                                                                   
// Arguments:                                                                                                                          
//   sim, ic, cosmo  sim, ic, cosmo structures   
//                 
//                 
//   viFTsub       reference to the Fourier image of the subtracted velocity
//   viFT          reference to the input Fourier image of the non-linear velocity field
//   thFT          reference to the Fourier image of the diverge of the velocity field.
//                 Here we store the Fourier image of the linear divergence field.
//   count         counter for the CLASS transfer function file
//
//   a             scale factor
//
// Returns:                                                                                                                        
//                                                                                                                                
////////////////////////// 

void subtract_velocity(metadata & sim, icsettings & ic, cosmology & cosmo,
                       Field<Cplx> & viFT_sub, Field<Cplx> & viFT, Field<Cplx> & thFT, 
                       int count = 1, Real a = 1.)
                       
{  

  const int linesize = viFT.lattice().size(1);
  rKSite k(viFT.lattice());  

  gsl_spline*thspline = NULL;
  gsl_spline*tk_psi    = NULL;
  gsl_spline*tk_th_cdm = NULL;
  gsl_spline*tk_th_b   = NULL;
  double * temp = NULL;
  size_t numpts3d = (size_t) sim.numpts * (size_t) sim.numpts * (size_t) sim.numpts;

  gsl_interp_accel*tk_psi_accel = gsl_interp_accel_alloc();
  gsl_interp_accel*tk_theta_accel = gsl_interp_accel_alloc();
  

  Real * gridk2;
  Cplx * kshift;
  Real k2;
  Real k_mod;

  char filename[100];  
  sprintf(filename, "/home/leporif7/Nbody/gevolution-1.1/output/Transfer_gevolution/test_z%d_tk.dat", count);
  loadTransferFunctionsVel(filename, tk_psi, tk_th_cdm, "cdm", sim.boxsize, cosmo.h);
  loadTransferFunctionsVel(filename, tk_psi, tk_th_b, "b", sim.boxsize, cosmo.h);
      
  
  if (tk_th_cdm->size != tk_th_b->size)
	{
	  COUT << " error: baryon transfer function line number mismatch!"\
	       << endl;
	  parallel.abortForce();
	}
   
  temp = (double *) malloc(tk_th_cdm->size * sizeof(double));
  thspline = gsl_spline_alloc(gsl_interp_cspline, tk_th_cdm->size);

  for (int i = 0; i < tk_th_cdm->size; i++)
   {
     cout << tk_th_cdm->y[i] << " " << tk_th_b->y[i] << "\n";      
     temp[i] = -((cosmo.Omega_cdm * tk_th_cdm->y[i] + cosmo.Omega_b* tk_th_b->y[i]) / 
              (cosmo.Omega_cdm + cosmo.Omega_b))*M_PI
              *sqrt(Pk_primordial_vort(tk_th_cdm->x[i]*cosmo.h/sim.boxsize, ic)/tk_th_cdm->x[i])/tk_th_cdm->x[i];
   }


  gridk2 = (Real *) malloc(linesize * sizeof(Real));
  kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

  for (int i = 0; i < linesize; i++)
    {
      gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
      kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
      gridk2[i] *= gridk2[i];
    }

  gsl_spline_init(thspline, tk_th_cdm->x, temp, tk_th_cdm->size);
  generateReal(thFT, 0., thspline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
 

  for (k.first(); k.test(); k.next())
    {
      k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
      k_mod = sqrt(k2);
      // Uncomment for computing the velocity from the CLASS linear transfer function
      //viFT_sub(k, 0) = -Cplx(0.0, 1.0)*kshift[k.coord(0)].conj()*a/k2*thFT(k)*numpts3d;
      //viFT_sub(k, 1) = -Cplx(0.0, 1.0)*kshift[k.coord(1)].conj()*a/k2*thFT(k)*numpts3d;
      //viFT_sub(k, 2) = -Cplx(0.0, 1.0)*kshift[k.coord(2)].conj()*a/k2*thFT(k)*numpts3d; 
      //      cout << "viTF, coorection: " << viFT(k, 0) <<    
      viFT_sub(k, 0) = viFT(k, 0) - Cplx(0.0, 1.0)*kshift[k.coord(0)].conj()*a/k2*thFT(k)*numpts3d;
      viFT_sub(k, 1) = viFT(k, 1) - Cplx(0.0, 1.0)*kshift[k.coord(1)].conj()*a/k2*thFT(k)*numpts3d;
      viFT_sub(k, 2) = viFT(k, 2) - Cplx(0.0, 1.0)*kshift[k.coord(2)].conj()*a/k2*thFT(k)*numpts3d;
        
    } 
 
  free(gridk2);
  free(kshift);
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Description:                                                                                                  
//   Compute the rotational part of the velocity field, in Fourier space
//                                                                                                            
//                                                                                                            
// Arguments:                                                                                                     
//                                                                                                           
//   vRFT          reference to the Fourier image of the rotational part of the velocity                     
//   viFT          reference to the input Fourier image of the velocity field                                  
//                                                                                                             
// Returns:                                                                                                     
//                                                                                                             
//////////////////////////

void projectFTvelocityVR(Field<Cplx> & vRFT, Field<Cplx> & viFT)
{

  const int linesize = vRFT.lattice().size(1);
  int i;
  Real * gridk2;
  Cplx * kshift;
  Real * gridk;
  rKSite k(vRFT.lattice());
  Real k2;
  Cplx tmp(0., 0.);

  gridk2 = (Real *) malloc(linesize * sizeof(Real));
  kshift = (Cplx *) malloc(linesize * sizeof(Cplx));
  gridk = (Real *) malloc(linesize * sizeof(Real));

  for (i = 0; i < linesize; i++)
    {
      gridk[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);
      kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
      gridk2[i] = gridk[i]*gridk[i];
    }


  k.first();
  if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
    {
      vRFT(k, 0) = Cplx(0.,0.);
      vRFT(k, 1) = Cplx(0.,0.);
      vRFT(k, 2) = Cplx(0.,0.);
      k.next();
    }
  for (; k.test(); k.next())
    {

      k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
      if ((k.coord(0) == 0 || k.coord(0) == linesize/2)&& 
          (k.coord(1) == 0 || k.coord(1) == linesize/2)&&
	  (k.coord(2) == 0 || k.coord(2) == linesize/2)  ) {
	vRFT(k, 0) = Cplx(0.,0.);
	vRFT(k, 1) = Cplx(0.,0.);
	vRFT(k, 2) = Cplx(0.,0.);
      }
      else {
      tmp = (gridk[k.coord(0)] * viFT(k, 0) + gridk[k.coord(1)] * viFT(k, 1) + gridk[k.coord(2)] * viFT(k, 2)) / k2;

      vRFT(k, 0) = (viFT(k, 0) - gridk[k.coord(0)] * tmp); 
      vRFT(k, 1) = (viFT(k, 1) - gridk[k.coord(1)] * tmp);
      vRFT(k, 2) = (viFT(k, 2) - gridk[k.coord(2)] * tmp);}
    }


  free(gridk2);
  free(kshift);

}

//////////////////////////                                                                                     
// projectFTvelocityTh                                                                                        
//////////////////////////                                                                                      
// Description:                                                                                                 
//   Compute the diverge of the velocity in Fourier space
//                                                                                                               
// Arguments:
//   thFT       reference to the Fourier image of the divergence of the velocity field
//   viFT       reference to the Fourier image of the velocity field                                             
//
// Returns:                                                                                                       
//                                                                                                                
//////////////////////////                                                                                                                

void projectFTvelocityTh(Field<Cplx> & thFT, Field<Cplx> & viFT)
{
  const int linesize = thFT.lattice().size(1);
  int i;
  Real * gridk2;
  Real * gridk;
  Cplx * kshift;
  rKSite k(thFT.lattice());
  Real k2;
  Cplx tmp(0., 0.);

  gridk2 = (Real *) malloc(linesize * sizeof(Real));
  kshift = (Cplx *) malloc(linesize * sizeof(Cplx));
  gridk = (Real *) malloc(linesize * sizeof(Real));
 
  for (i = 0; i < linesize; i++)
    {
      gridk[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);
      kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
      gridk2[i] = gridk[i]*gridk[i];
    }
  

  for (k.first(); k.test(); k.next())
    {
      thFT(k) = Cplx(0.,1.)*(gridk[k.coord(0)] * viFT(k, 0) + 
                 gridk[k.coord(1)] * viFT(k, 1) + 
                 gridk[k.coord(2)] * viFT(k, 2) );
    }

  free(gridk2);
  free(kshift);
  free(gridk);
}

//////////////////////////                                                                                                        
// compute_vi_zero
//////////////////////////                                                                                                         
// Description:                                                                                                                   
//   Compute the velocity field as v^i = T^i_0/T^0_0, if a = 1 then vi = a v^i
//   If T^0_0 = 0 the velocity is set to zero (velocity method = zero)                                                   
//                                                                                                                           
// Arguments:                                                                                                                 
//   viFT       reference to the velocity field                                                                
//   source     reference to the field source (a^3 T^0_0)
//   Bi         reference to the field Bi (a^4 T^0_i)
//   phi        reference to the field phi
//   chi        reference to the field chi 
// Returns:                                                                                                            
//                                                                                                                                 
////////////////////////// 

void compute_vi_zero(Field<Real> * vi, Field<Real> * source = NULL, Field<Real> * Ti0 = NULL)
{
  
  Site xvi(vi->lattice());
      
  for(xvi.first(); xvi.test(); xvi.next())
    {  


      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,0)= 0.0;}
      else {(*vi)(xvi,0) = (*Ti0)(xvi,0)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,1)= 0.0;}
      else {(*vi)(xvi,1) = (*Ti0)(xvi,1)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,2)= 0.0;}
      else {(*vi)(xvi,2) = (*Ti0)(xvi,2)/(*source)(xvi);}

    }
}

//////////////////////////                                                                                                      
// compute_vi_past                                                                                                                 
//////////////////////////                                                                                                        
// Description:                                                                                                               
//   Compute the velocity field as v^i = T^i_0/T^0_0, if a = 1 then vi = a v^i                                                 
//   If T^0_0 = 0 the velocity field is set to be the one at the previous time step (velocity method = past)             
//                                                                                                                              
// Arguments:                                                                                                                    
//   viFT       reference to the velocity field                                                                                    
//   source     reference to the field source (a^3 T^0_0)                                                                      
//   Bi         reference to the field Bi (a^4 T^0_i)                                                                           
//   phi        reference to the field phi                                                                                      
//   chi        reference to the field chi
//   vi_past    reference to the velocity field at the previous time step                                               
// Returns:                                                                                                                    
//                                                                                                                                 
//////////////////////////


void compute_vi_past(Field<Real> * vi, Field<Real> * source = NULL, Field<Real> * Ti0 = NULL, Field<Real> * vi_past = NULL)
{

  Real  localCubePhi[8];
  Real  localCubeChi[8];
  Real  localCubeT00[8];
  Real  localEdgeTi0[12];  

  Site xvi(vi->lattice());
      
  for(xvi.first(); xvi.test(); xvi.next())
    {  

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,0)= (*vi_past)(xvi,0);}
      else {(*vi)(xvi,0) = (*Ti0)(xvi,0)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,1)= (*vi_past)(xvi,1);}
      else {(*vi)(xvi,1) = (*Ti0)(xvi,1)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,2)= (*vi_past)(xvi,2);}
      else {(*vi)(xvi,2) = (*Ti0)(xvi,2)/(*source)(xvi);}

    }
}

//////////////////////////                                                                                         
// compute_vi_past_rescaled                                                                                                     
//////////////////////////                                                                                                    
// Description:                                                                                                                
//   Compute the velocity field as v^i = T^i_0/T^0_0, if a = 1 then vi = a v^i                                              
//   If T^0_0 = 0 the velocity field is set to be the one at the previous time step,
//   rescaled as v^i(a) = v^i(a_past) a*Hconf(a) dD1/da (velocity method = rescaled past)                  
//                                                                                                                            
// Arguments:                                                                                                                 
//   viFT       reference to the velocity field                                                                                    
//   source     reference to the field source (a^3 T^0_0)                                                                       
//   Bi         reference to the field Bi (a^4 T^0_i)                                                                           
//   phi        reference to the field phi                                                                                      
//   chi        reference to the field chi                                                                                    
//   vi_past    reference to the velocity field at the previous time step                                                     
// Returns:                                                                                                                    
//                                                                                                                                 
////////////////////////// 

void compute_vi_past_rescaled(cosmology & cosmo, Field<Real> * vi, Field<Real> * source = NULL, double a = 1., double a_past = 1., Field<Real> * Ti0 = NULL, Field<Real> * vi_past = NULL)
{

  Real  localCubePhi[8];
  Real  localCubeChi[8];
  Real  localCubeT00[8];
  Real  localEdgeTi0[12];  

  Site xvi(vi->lattice());

  Real rescale = D1_prime(cosmo, a)/D1_prime(cosmo, a_past)*a/a_past;
  
  for(xvi.first(); xvi.test(); xvi.next())
    {  

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,0)= (*vi_past)(xvi,0)*rescale;}
      else {(*vi)(xvi,0) = (*Ti0)(xvi,0)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,1)= (*vi_past)(xvi,1)*rescale;}
      else {(*vi)(xvi,1) = (*Ti0)(xvi,1)/(*source)(xvi);}

      if ( (*source)(xvi) < 1.E-300) {(*vi)(xvi,2)= (*vi_past)(xvi,2)*rescale;}
      else {(*vi)(xvi,2) = (*Ti0)(xvi,2)/(*source)(xvi);}
    }
}

// Store the velocity field at each time step in vi_past 
void store_vi(Field<Real> * vi_past, Field<Real> * vi = NULL)
{
  Site x(vi_past->lattice());

  if (vi != NULL)
  {
   for(x.first(); x.test(); x.next())
    {
      (*vi_past)(x,0) = (*vi)(x,0);
      (*vi_past)(x,1) = (*vi)(x,1);
      (*vi_past)(x,2) = (*vi)(x,2);
    }
  }
}

//////////////////////////                                                                                                             
// compute_count
//////////////////////////                                                                                                        
// Description:                                                                                                                        
//   Compute the number of particles in each cell, stored in part_in_cube.
//   Compute the number of empty cells.  
//                                                                                                                                     
// Arguments:                                                                                                                          
//   part          particles information
//   n_empty_size  number of empty cells
//   part_in_cube  reference to the field containing the number of particles in cell
//
// Returns:                                                                                                                          
//                                                                                                                              
////////////////////////// 

template<typename part, typename part_info, typename part_dataType>
void compute_count(Particles<part,part_info,part_dataType> * pcls, long* n_empty_size, Field<Real> *part_in_cube = NULL)
{
  typename std::list<part>::iterator it;
  Site xPart(pcls->lattice());
  Site xField(part_in_cube->lattice());

  for(xField.first(); xField.test(); xField.next())
    {
      (*part_in_cube)(xField) = 0.;
    }

  for(xPart.first(), xField.first(); xPart.test(); xPart.next(), xField.next())
    {
      (*part_in_cube)(xField) = (pcls->field())(xPart).size;
      if ((pcls->field())(xPart).size == 0) (*n_empty_size) += 1;
    }

  parallel.sum<long>(*n_empty_size); 

}

//////////////////////////                                                                                                         
// convolve_field
//////////////////////////                                                                                                          
// Description:                                                                                                                     
//   Compute the smoothed field in Fourier space through convolution.
//                                                                                                                                 
// Arguments:                                                                                                                    
//   fieldFT_conv          reference to the Fourier image of the convolved field
//   fieldFT               reference to the Fourier image of the field to be smoothed
//   dim_field             3 for vector fields, 1 for scalar fields
//   sigma                 size of the smoothing, in code units
//                                                                                                                                
// Returns:                                                                                                                      
//                                                                        
//////////////////////////

void convolve_field(Field<Cplx> * fieldFT_conv, Field<Cplx> * fieldFT, int dim_field, Real sigma) 
{

        rKSite k(fieldFT->lattice());

	const int linesize = fieldFT->lattice().size(1);
	int i;
	Real * gridk2;
	Real k2;
	
	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	
	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		gridk2[i] *= gridk2[i];
	}
	
	
	for (k.first(); k.test(); k.next())
	{
		k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
		for (i = 0; i < dim_field; i++){
		  (*fieldFT_conv)(k, i) = (*fieldFT)(k, i)*exp(-0.5*k2*sigma*sigma); 
		}
	}
	
	free(gridk2);

}

//////////////////////////                                                                                     
// compute_velocity_smooth
//////////////////////////                                                                                     
// Description:                                                                                              
//   Compute the smoothed velocity field from the smoothed fields Ti0 and T00
//                                                                                                                   
// Arguments:                                                                                                        
//   vi         reference to the velocity field vi 
//   Ti0        reference to the field Ti0 = a^4 T^i_0                                                         
//   T00        reference to the field T00 (a^3 T^0_0)                                        
//                                                                                                                 
// Returns:                                                                             
//                                                                                                               
//////////////////////////

void compute_velocity_smooth(Field<Real> * vi, Field<Real> * Ti0 = NULL, Field<Real> * T00 = NULL)
{
  Site xvi(vi->lattice());
      
  for(xvi.first(); xvi.test(); xvi.next())
    {  
      (*vi)(xvi,0) = (*Ti0)(xvi,0)/(*T00)(xvi);
      (*vi)(xvi,1) = (*Ti0)(xvi,1)/(*T00)(xvi);
      (*vi)(xvi,2) = (*Ti0)(xvi,2)/(*T00)(xvi);
    }
}


void compute_norm2_vR(
		      Field<Real> * vR,
		      Field<Real> * norm_vR2
		      )
{
  Site x(norm_vR2->lattice());
  for(x.first(); x.test(); x.next()){
    (*norm_vR2)(x) = pow((*vR)(x, 0),2)+ pow((*vR)(x, 1),2)+ pow((*vR)(x, 2),2);
  }
}

void compute_norm_w_old(
		    Field<Real> &vR,
		    Field<Real> &norm2_vR
		    )
{
  Site x(vR.lattice());
  for(x.first(); x.test(); x.next()){
    (norm2_vR)(x) = vR(x,0)* vR(x,0) + vR(x,1)* vR(x,1) + vR(x,2)* vR(x,2);
  }
}

void compute_norm_w(
                    Field<Real> * norm2_vR,
                    Field<Real> * norm_w
                    )
{
  Site x(norm_w->lattice());
  for(x.first(); x.test(); x.next()){
    (*norm_w)(x) = 8.0*(*norm2_vR)(x)-((*norm2_vR)(x-0-1-2) + (*norm2_vR)(x-0-1+2) + (*norm2_vR)(x-0+1-2) + (*norm2_vR)(x+0-1-2) + \
				       (*norm2_vR)(x-0+1+2) + (*norm2_vR)(x+0-1+2) + (*norm2_vR)(x+0+1-2) + (*norm2_vR)(x+0+1+2) );
    (*norm_w)(x) = sqrt(abs( (*norm_w)(x) ));

  }
}


void compute_v2(             
		    Field<Real> * v2,
                    Field<Real> * vi
                    )
{
  Site x(v2->lattice());
  for(x.first(); x.test(); x.next()){
    (*v2)(x) = (*vi)(x,0)*(*vi)(x,0)+(*vi)(x,1)*(*vi)(x,1)+(*vi)(x,2)*(*vi)(x,2) ;
  }
}



void compute_laplacianFT(
			 Field<Cplx> &source_scalar,
			 Field<Cplx> &dest_scalar
			 )
{
  const int linesize = source_scalar.lattice().size(1);
  int i;
  Real *gridk2;
  rKSite k(source_scalar.lattice());

  gridk2 = (Real *) malloc(linesize * sizeof(Real));

  double coeff = ((long) linesize * (long) linesize * (long) linesize);

  for (i = 0; i < linesize; i++)
    {
      gridk2[i] = (Real) linesize * sin(M_PI * 2.0 * (Real) i / (Real) linesize);
      gridk2[i] *= gridk2[i];
    }

  for (k.first(); k.test(); k.next())
    {
      dest_scalar(k) = -source_scalar(k) * coeff *(gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]);
    }

  free(gridk2);
}


//////////////////////////                                                                                                                    
// compute_sigma2_rescaled                                                                                                                    
//////////////////////////                                                                                                                    
// Description:                                                                                                                              
//   Compute the trace of the velocity dispertion tensor as sigma2 = -(T^1_1 + T^2_2 + T^3_3) / T^0_0 - <v>^2
//   If T^0_0 = 0 the velocity dispertion is set to be the one at the previous time step rescaled by linear velocity growth      
//                                                                                                                                            
// Arguments:                                                                                                                              
//   sigma2         reference to the velocity dispertion scalare
//   source         reference to the field source (a^3 T^0_0)                                                                
//   Sij            reference to the field Sij (a^3 T^i_j for the diagonal)
//   vi             reference to the velocity field                                                                     
//   sigma2_past    reference to the velocity dispertion at the previous time step                                                    
// Returns:                                                                                                                                   
//                                                                                                                                            
//////////////////////////                                                                                                                    

void compute_sigma2_rescaled(cosmology & cosmo, Field<Real> * sigma2, Field<Real> * source = NULL, Field<Real> * Sij = NULL, Field<Real> * vi = NULL, Field<Real> * sigma2_past = NULL, double a = 1., double a_past = 1.)
{

  Site xsigma(sigma2->lattice());

  Real rescale = D1_prime(cosmo, a)/D1_prime(cosmo, a_past)*a/a_past;

  for(xsigma.first(); xsigma.test(); xsigma.next())
    {

      if ( (*source)(xsigma) < 1.E-300) {(*sigma2)(xsigma)= (*sigma2_past)(xsigma)*rescale*rescale;}
      else {(*sigma2)(xsigma) = ((*Sij)(xsigma, 0, 0) + (*Sij)(xsigma, 1, 1) + (*Sij)(xsigma, 2, 2) )/(*source)(xsigma) 
	  -((*vi)(xsigma,0)*(*vi)(xsigma,0) + (*vi)(xsigma,1)*(*vi)(xsigma,1) + (*vi)(xsigma,2)*(*vi)(xsigma,2));
           }

    }
}


// Store the velocity field at each time step in vi_past                                                                                        
void store_sigma2(Field<Real> * sigma2_past, Field<Real> * sigma2 = NULL)
{
  Site x(sigma2_past->lattice());

  if (sigma2 != NULL)
    {
      for(x.first(); x.test(); x.next())
	{
	  (*sigma2_past)(x,0) = (*sigma2)(x,0);
	  (*sigma2_past)(x,1) = (*sigma2)(x,1);
	  (*sigma2_past)(x,2) = (*sigma2)(x,2);
	}
    }
}

#endif
