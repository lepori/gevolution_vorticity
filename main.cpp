//////////////////////////
// Copyright (c) 2015-2016 Julian Adamek (Université de Genève)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//  
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//  
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////

//////////////////////////
// main.cpp
//////////////////////////
// 
// main control sequence of Geneva N-body code with evolution of metric perturbations (gevolution)
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: February 2017
//
//////////////////////////

#include <stdlib.h>
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX            // due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "class_tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "tools.hpp"
#include "output.hpp"
#include "hibernation.hpp"

using namespace std;

using namespace LATfield2;

int main(int argc, char **argv)
{
    
#ifdef BENCHMARK
    //benchmarking variables
    
    double ref_time, ref2_time, cycle_start_time;
    
    double initialization_time;
    double run_time;
    double cycle_time=0;
    double projection_time = 0;
    double snapshot_output_time = 0;
    double spectra_output_time = 0;
    double gravity_solver_time = 0;
    double fft_time = 0;
    int fft_count = 0;   
    double update_q_time = 0;
    int update_q_count = 0;
    double moveParts_time = 0;
    int  moveParts_count =0;
    
#endif  //BENCHMARK
    
    int n = 0, m = 0;
    int io_size = 0;
    int io_group_size = 0;
    
    int i, j, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numspecies;
    int numsteps_ncdm[MAX_PCL_SPECIES-2];
    long numpts3d;
    int box[3];
    double dtau, dtau_old, dx, tau, a, fourpiG, tau_Lambda, tmp, start_time;
    double maxvel[MAX_PCL_SPECIES];
    FILE * outfile;
    char filename[2*PARAM_MAX_LENGTH+24];
    string h5filename;
    char * settingsfile = NULL;
    parameter * params = NULL;
    metadata sim;
    cosmology cosmo;
    icsettings ic;
    gadget2_header hdr;
    Real T00hom;

#ifndef H5_DEBUG
    H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif
    
    for (i=1 ; i < argc ; i++ ){
        if ( argv[i][0] != '-' )
            continue;
        switch(argv[i][1]) {
            case 's':
                settingsfile = argv[++i]; //settings file name
                break;
            case 'n':
                n = atoi(argv[++i]); //size of the dim 1 of the processor grid
                break;
            case 'm':
                m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
                break;
            case 'i':
#ifndef EXTERNAL_IO
                cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
                exit(-1000);
#endif
                io_size =  atoi(argv[++i]);
                break;
            case 'g':
#ifndef EXTERNAL_IO
                cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
                exit(-1000);
#endif
                io_group_size = atoi(argv[++i]);
        }
    }

#ifndef EXTERNAL_IO
    parallel.initialize(n,m);
#else
    parallel.initialize(n,m,io_size,io_group_size);
    if(parallel.isIO()) ioserver.start();
    else
    {
#endif
    
    COUT << COLORTEXT_WHITE << endl;    
    COUT << "  _   _      _         __ ,  _" << endl;
    COUT << " (_| (-' \\/ (_) (_ (_| (  ( (_) /\\/    version 1.1         running on " << n*m << " cores." << endl;
    COUT << "  -'" << endl << COLORTEXT_RESET << endl;
    
    if (settingsfile == NULL)
    {
        COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
        parallel.abortForce();
    }
    
    COUT << " initializing..." << endl;
    
    start_time = MPI_Wtime();
    
    numparam = loadParameterFile(settingsfile, params);
    
    usedparams = parseMetadata(params, numparam, sim, cosmo, ic);
    
    COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;
    
    sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
    saveParameterFile(filename, params, numparam);
    
    free(params);

#ifdef HAVE_CLASS
    background class_background;
      perturbs class_perturbs;
      spectra class_spectra;
#endif
    
    h5filename.reserve(2*PARAM_MAX_LENGTH);
    h5filename.assign(sim.output_path);
    h5filename += sim.basename_snapshot;
    
    box[0] = sim.numpts;
    box[1] = sim.numpts;
    box[2] = sim.numpts;
    
    Lattice lat(3,box,1);
    Lattice latFT;
    latFT.initializeRealFFT(lat,0);
    
    Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
    Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_b;
    Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES-2];
    Field<Real> * update_cdm_fields[3];
    Field<Real> * update_b_fields[3];
    Field<Real> * update_ncdm_fields[3];
    double f_params[5];

    Field<Real> phi;
    Field<Real> source;
    Field<Real> chi;
    Field<Real> Sij;
    Field<Real> Bi;
        Field<Real> vi;
    Field<Cplx> scalarFT;
    Field<Cplx> SijFT;
    Field<Cplx> BiFT;
        Field<Cplx> viFT;
     
    source.initialize(lat,1);
    phi.initialize(lat,1);
    chi.initialize(lat,1);
    scalarFT.initialize(latFT,1);
    PlanFFT<Cplx> plan_source(&source, &scalarFT);
    PlanFFT<Cplx> plan_phi(&phi, &scalarFT);
    PlanFFT<Cplx> plan_chi(&chi, &scalarFT);
    Sij.initialize(lat,3,3,symmetric);
    SijFT.initialize(latFT,3,3,symmetric);
    PlanFFT<Cplx> plan_Sij(&Sij, &SijFT);
        vi.initialize(lat,3);
    Bi.initialize(lat,3);
        viFT.initialize(latFT,3);
    BiFT.initialize(latFT,3);
    PlanFFT<Cplx> plan_Bi(&Bi, &BiFT);
        PlanFFT<Cplx> plan_vi(&vi, &viFT);
        //PlanFFT<Cplx> plan_wi(&wiFT);

#ifdef CHECK_B
    Field<Real> Bi_check;
        Field<Real> vi_check;
    Field<Cplx> BiFT_check;
        Field<Cplx> viFT_check;
    Bi_check.initialize(lat,3);
        vi_check.initialize(lat,3);
    BiFT_check.initialize(latFT,3);
        viFT_check.initialize(latFT,3);
    PlanFFT<Cplx> plan_Bi_check(&Bi_check, &BiFT_check);
        PlanFFT<Cplx> plan_vi_check(&vi_check, &viFT_check);
#endif

    update_cdm_fields[0] = &phi;
    update_cdm_fields[1] = &chi;
    update_cdm_fields[2] = &Bi;
        
    update_b_fields[0] = &phi;
    update_b_fields[1] = &chi;
    update_b_fields[2] = &Bi;
            
    update_ncdm_fields[0] = &phi;
    update_ncdm_fields[1] = &chi;
    update_ncdm_fields[2] = &Bi;
        
    Site x(lat);
    rKSite kFT(latFT);
    
        //projection_init(&NORM_smooth);
    dx = 1.0 / (double) sim.numpts;
    numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

        //sigma smoothing                                                                                                                 
    /*
        Real NORM_smooth[sim.numpts][sim.numpts][sim.numpts]= {0.};
        Real sigma_smooth = 0.1*dx;
        Site y(lat);
        for (x.first(); x.test(); x.next())
          { //cout << "x0 " << x.coord(0) << "y0 " << x.coord(1) << "z0 " << x.coord(2)<< "\n";                                  
            for (y.first(); y.test(); y.next()){
              cout << "Computing some shit \n";
              NORM_smooth[x.coord(0)][x.coord(1)][x.coord(2)] += exp(-(
                                       pow(x.coord(0)-y.coord(0),2) +
                                       pow(x.coord(1)-y.coord(1),2) +
                                       pow(x.coord(2)-y.coord(2),2)
                                       ) /(2.*pow(sigma_smooth,2))   );
            }
          }    
     */

    for (i = 0; i < 3; i++) // particles may never move farther than to the adjacent domain
    {
        if (lat.sizeLocal(i)-1 < sim.movelimit)
            sim.movelimit = lat.sizeLocal(i)-1;
    }
    parallel.min(sim.movelimit);

    fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
    a = 1. / (1. + sim.z_in);
    tau = particleHorizon(a, fourpiG, cosmo);
    tau_Lambda = -1.0;
    
    if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
        dtau = sim.Cf * dx;
    else
        dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
        
    dtau_old = 0.;
    
    if (ic.generator == ICGEN_BASIC)
      generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi,&vi, &source, &Sij, &scalarFT, &BiFT, &viFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi,&plan_vi, &plan_source, &plan_Sij); // generates ICs on the fly
    else if (ic.generator == ICGEN_READ_FROM_DISK)
        readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
#ifdef ICGEN_PREVOLUTION
    else if (ic.generator == ICGEN_PREVOLUTION)
        generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
#ifdef ICGEN_FALCONIC
    else if (ic.generator == ICGEN_FALCONIC)
        maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
    else
    {
        COUT << " error: IC generator not implemented!" << endl;
        parallel.abortForce();
    }
    
    if (sim.baryon_flag > 1)
    {
        COUT << " error: baryon_flag > 1 after IC generation, something went wrong in IC generator!" << endl;
        parallel.abortForce();
    }
    
    numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;    
    parallel.max<double>(maxvel, numspecies);
    
    if (sim.gr_flag > 0)
    {
        for (i = 0; i < numspecies; i++)
            maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
    }

#ifdef CHECK_B
    if (sim.vector_flag == VECTOR_ELLIPTIC)
    {
        for (kFT.first(); kFT.test(); kFT.next())
        {
            BiFT_check(kFT, 0) = BiFT(kFT, 0);
            BiFT_check(kFT, 1) = BiFT(kFT, 1);
            BiFT_check(kFT, 2) = BiFT(kFT, 2);
        }
    }
#endif

#ifdef CHECK_v
        if (sim.vector_flag == VECTOR_ELLIPTIC)
      {
        for (kFT.first(); kFT.test(); kFT.next())
          {
        viFT_check(kFT, 0) = viFT(kFT, 0);
        viFT_check(kFT, 1) = viFT(kFT, 1);
        viFT_check(kFT, 2) = viFT(kFT, 2);
          }
      }
#endif
    
    for (i = 0; i < 6; i++)
    {
        hdr.npart[i] = 0;
        hdr.npartTotal[i] = 0;
        hdr.mass[i] = 0.;
    }
    hdr.num_files = 1;
    hdr.Omega0 = cosmo.Omega_m;
    hdr.OmegaLambda = cosmo.Omega_Lambda;
    hdr.HubbleParam = cosmo.h;
    hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
    hdr.flag_sfr = 0;
    hdr.flag_cooling = 0;
    hdr.flag_feedback = 0;
    for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; i++)
        hdr.fill[i] = 0;
    
#ifdef BENCHMARK
    initialization_time = MPI_Wtime() - start_time;
    parallel.sum(initialization_time);
    COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
    COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

#ifdef HAVE_CLASS
    if (sim.radiation_flag > 0)
    {
        initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra);
        if (sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
        {
            prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
            plan_source.execute(FFT_BACKWARD);
            for (x.first(); x.test(); x.next())
                chi(x) += source(x);
            chi.updateHalo();
        }
    }
#endif

    while (true)    // main loop
    {
#ifdef BENCHMARK        
        cycle_start_time = MPI_Wtime();
#endif
        // construct stress-energy tensor
        projection_init(&source);
#ifdef HAVE_CLASS
        if (sim.radiation_flag > 0)
            projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
#endif
        if (sim.gr_flag > 0)
        {
            projection_T00_project(&pcls_cdm, &source, a, &phi);
            if (sim.baryon_flag)
                projection_T00_project(&pcls_b, &source, a, &phi);
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
                    projection_T00_project(pcls_ncdm+i, &source, a, &phi);
                else if (sim.radiation_flag == 0)
                {
                    tmp = bg_ncdm(a, cosmo, i);
                    for(x.first(); x.test(); x.next())
                        source(x) += tmp;
                }
            }
        }
        else
        {
            scalarProjectionCIC_project(&pcls_cdm, &source);
            if (sim.baryon_flag)
                scalarProjectionCIC_project(&pcls_b, &source);
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
                    scalarProjectionCIC_project(pcls_ncdm+i, &source);
            }
        }
        projection_T00_comm(&source);

        if (sim.vector_flag == VECTOR_ELLIPTIC)
          {     
            projection_init(&Bi);
            projection_T0i_project(&pcls_cdm, &Bi, &phi);
            if (sim.baryon_flag)
                projection_T0i_project(&pcls_b, &Bi, &phi);
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.))
                    projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
            }
            projection_T0i_comm(&Bi);
        }

                if (sim.vector_flag == VECTOR_ELLIPTIC)
 
          {                     
                    compute_vi_project(&vi, &source, a, &Bi,  &phi);
                  }  

        
        projection_init(&Sij);
        projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
        if (sim.baryon_flag)
            projection_Tij_project(&pcls_b, &Sij, a, &phi);
        if (a >= 1. / (sim.z_switch_linearchi + 1.))
        {
            for (i = 0; i < cosmo.num_ncdm; i++)
                projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
        }
        projection_Tij_comm(&Sij);
        //Here we add something to compute the velocity field
        //compute_vi_project

        
#ifdef BENCHMARK 
        projection_time += MPI_Wtime() - cycle_start_time;
        ref_time = MPI_Wtime();
#endif
        
        if (sim.gr_flag > 0)
        {    
            T00hom = 0.;
            for (x.first(); x.test(); x.next())
                T00hom += source(x);
            parallel.sum<Real>(T00hom);
            T00hom /= (Real) numpts3d;
            
            if (cycle % CYCLE_INFO_INTERVAL == 0)
            {
                COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
            }
            
            if (dtau_old > 0.)
            {
                prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update

#ifdef BENCHMARK
                ref2_time= MPI_Wtime();
#endif
                plan_source.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
                fft_time += MPI_Wtime() - ref2_time;
                fft_count++;
#endif
        
                solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);  // phi update (k-space)

#ifdef BENCHMARK
                ref2_time= MPI_Wtime();
#endif        
                plan_phi.execute(FFT_BACKWARD);     // go back to position space
#ifdef BENCHMARK
                fft_time += MPI_Wtime() - ref2_time;
                fft_count++;
#endif    
            }
        }
        else
        {
#ifdef BENCHMARK
            ref2_time= MPI_Wtime();
#endif
            plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space
#ifdef BENCHMARK
            fft_time += MPI_Wtime() - ref2_time;
            fft_count++;
#endif
        
            solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)

#ifdef BENCHMARK
            ref2_time= MPI_Wtime();
#endif        
            plan_phi.execute(FFT_BACKWARD);     // go back to position space
#ifdef BENCHMARK
            fft_time += MPI_Wtime() - ref2_time;
            fft_count++;
#endif    
        }

        phi.updateHalo();  // communicate halo values

        // record some background data
        if (kFT.setCoord(0, 0, 0))
        {
            sprintf(filename, "%s%s_background.dat", sim.output_path, sim.basename_generic);
            outfile = fopen(filename, "a");
            if (outfile == NULL)
            {
                cout << " error opening file for background output!" << endl;
            }
            else
            {
                if (cycle == 0)
                    fprintf(outfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
                fprintf(outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau, a, Hconf(a, fourpiG, cosmo) / Hconf(1., fourpiG, cosmo), scalarFT(kFT).real(), T00hom);
                fclose(outfile);
            }
        }
        // done recording background data
        
        prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations

#ifdef BENCHMARK
        ref2_time= MPI_Wtime();
#endif        
        plan_Sij.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
        fft_time += MPI_Wtime() - ref2_time;
        fft_count += 6;
#endif

#ifdef HAVE_CLASS
        if (sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
        {
            prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
            projectFTscalar(SijFT, scalarFT, 1);
        }
        else
#endif        
        projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)

#ifdef BENCHMARK
        ref2_time= MPI_Wtime();
#endif        
        plan_chi.execute(FFT_BACKWARD);     // go back to position space
#ifdef BENCHMARK
        fft_time += MPI_Wtime() - ref2_time;
        fft_count++;
#endif    
        chi.updateHalo();  // communicate halo values

        if (sim.vector_flag == VECTOR_ELLIPTIC)
        {
#ifdef BENCHMARK
            ref2_time= MPI_Wtime();
#endif
            plan_Bi.execute(FFT_FORWARD);
                        plan_vi.execute(FFT_FORWARD);
#ifdef BENCHMARK
            fft_time += MPI_Wtime() - ref2_time;
            fft_count++;
#endif
            projectFTvector(BiFT, BiFT, fourpiG * dx * dx); // solve B using elliptic constraint (k-space)
                projectFTvelocity(viFT, viFT, dx*dx);           // compute the vorticity field 
#ifdef CHECK_B
            evolveFTvector(SijFT, BiFT_check, a * a * dtau_old); 
#endif
        }
        else
            evolveFTvector(SijFT, BiFT, a * a * dtau_old);  // evolve B using vector projection (k-space)

        if (sim.gr_flag > 0)
        {
#ifdef BENCHMARK
            ref2_time= MPI_Wtime();
#endif                
            plan_Bi.execute(FFT_BACKWARD);  // go back to position space
                        plan_vi.execute(FFT_BACKWARD); 
#ifdef BENCHMARK
            fft_time += MPI_Wtime() - ref2_time;
            fft_count += 3;
#endif
             Bi.updateHalo();  // communicate halo values
                        vi.updateHalo(); 
        }

#ifdef BENCHMARK 
        gravity_solver_time += MPI_Wtime() - ref_time;
        ref_time = MPI_Wtime();
#endif

        // snapshot output
        if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
        {
            COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
            writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
            writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif

            snapcount++;
        }
        
#ifdef BENCHMARK
        snapshot_output_time += MPI_Wtime() - ref_time;
        ref_time = MPI_Wtime();
#endif
        
        // power spectra
        if (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
        {
            COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B            
            writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &vi,&viFT, &plan_vi, &Bi_check, &BiFT_check, &plan_Bi_check, &vi_check, &viFT_check, &plan_vi_check );
#else
            writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &vi,&viFT, &plan_vi);
#endif

            pkcount++;
        }
        
#ifdef BENCHMARK
        spectra_output_time += MPI_Wtime() - ref_time;
#endif 

        if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot) break; // simulation complete
        
        // compute number of step subdivisions for particle updates
        numsteps = 1;
        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
                numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
            else numsteps_ncdm[i] = 1;
            
            if (numsteps < numsteps_ncdm[i]) numsteps = numsteps_ncdm[i];
        }
        if (numsteps > 1 && numsteps % 2 > 0) numsteps++;   // if >1, make it an even number
        
        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (numsteps / numsteps_ncdm[i] <= 1) numsteps_ncdm[i] = numsteps;
            else if (numsteps_ncdm[i] > 1) numsteps_ncdm[i] = numsteps / 2;
        }
        
        if (cycle % CYCLE_INFO_INTERVAL == 0)
        {
            COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
            if (sim.baryon_flag)
            {
                COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
            }
            
            COUT << "), time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau;
            
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (i == 0)
                {
                    COUT << endl << " time step subdivision for ncdm species: ";
                }
                COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i+1+sim.baryon_flag] << ")";
                if (i < cosmo.num_ncdm-1)
                {
                    COUT << ", ";
                }
            }
            
            COUT << endl;
        }

        for (j = 0; j < numsteps; j++) // particle update
        {
#ifdef BENCHMARK
            ref2_time = MPI_Wtime();
#endif
            f_params[0] = a;
            f_params[1] = a * a * sim.numpts;
            if (j == 0)
            {
                if (sim.gr_flag > 0)
                {
                    maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
                    if (sim.baryon_flag)
                        maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
                }
                else
                {
                    maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
                    if (sim.baryon_flag)
                        maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
                }

#ifdef BENCHMARK
                update_q_count++;
#endif
            }
                
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (j % (numsteps / numsteps_ncdm[i]) == 0)
                {
                    if (sim.gr_flag > 0)
                        maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
                    else
                        maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

#ifdef BENCHMARK
                    update_q_count++;
#endif
                }
            }
#ifdef BENCHMARK
            update_q_time += MPI_Wtime() - ref2_time;
            ref2_time = MPI_Wtime();
#endif
        
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (numsteps > 1 && ((numsteps_ncdm[i] == 1 && j == numsteps / 2) || (numsteps_ncdm[i] == numsteps / 2 && j % 2 > 0)))
                {
                    if (sim.gr_flag > 0)
                        pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
                    else
                        pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
                        moveParts_count++;
                        moveParts_time += MPI_Wtime() - ref2_time;
                        ref2_time = MPI_Wtime();
#endif
                }      
            }

            if (numsteps == 1)
                rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step

            f_params[0] = a;
            f_params[1] = a * a * sim.numpts;
            if (numsteps == 1 || j == numsteps / 2)
            {
                if (sim.gr_flag > 0)
                {
                    pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
                    if (sim.baryon_flag)
                        pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
                }
                else
                {
                    pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
                    if (sim.baryon_flag)
                        pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
                }

#ifdef BENCHMARK
                moveParts_count++;
                moveParts_time += MPI_Wtime() - ref2_time;
                ref2_time = MPI_Wtime();
#endif
            }

            if (numsteps != 1)
                rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step
            
            f_params[0] = a;
            f_params[1] = a * a * sim.numpts;    
            for (i = 0; i < cosmo.num_ncdm; i++)
            {
                if (numsteps_ncdm[i] == numsteps)
                {
                    if (sim.gr_flag > 0)
                        pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
                    else
                        pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
                        moveParts_count++;
                        moveParts_time += MPI_Wtime() - ref2_time;
                        ref2_time = MPI_Wtime();
#endif
                }
            }

            rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step
        }   // particle update done
        
        parallel.max<double>(maxvel, numspecies);
        
        if (sim.gr_flag > 0)
        {
            for (i = 0; i < numspecies; i++)
                maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
        }
        
        tau += dtau;
        
        if (tau_Lambda < 0. && (cosmo.Omega_m / a / a / a) < cosmo.Omega_Lambda)
        {
            tau_Lambda = tau;
            COUT << "matter-dark energy equality at z=" << ((1./a) - 1.) << endl;
        }
        
        if (sim.wallclocklimit > 0.)   // check for wallclock time limit
        {
            tmp = MPI_Wtime() - start_time;
            parallel.max(tmp);
            if (tmp > sim.wallclocklimit)   // hibernate
            {
                COUT << COLORTEXT_YELLOW << " reaching hibernation wallclock limit, hibernating..." << COLORTEXT_RESET << endl;
                COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
                if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
                    plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
                if (sim.vector_flag == VECTOR_ELLIPTIC)
                {
                    plan_Bi_check.execute(FFT_BACKWARD);
                    hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle);
                }
                else
#endif
                hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle);
                break;
            }
        }
        
        if (restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
        {
            COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
            if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
                plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
            if (sim.vector_flag == VECTOR_ELLIPTIC)
            {
                plan_Bi_check.execute(FFT_BACKWARD);
                hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle, restartcount);
            }
            else
#endif
            hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle, restartcount);
            restartcount++;
        }
        
        dtau_old = dtau;
        
        if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
            dtau = sim.Cf * dx;
        else
            dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
           
        cycle++;
        
#ifdef BENCHMARK
        cycle_time += MPI_Wtime()-cycle_start_time;
#endif
    }
    
    COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef HAVE_CLASS
    if (sim.radiation_flag > 0)
        freeCLASSstructures(class_background, class_perturbs, class_spectra);
#endif

#ifdef BENCHMARK
    run_time = MPI_Wtime() - start_time;

    parallel.sum(run_time);
    parallel.sum(cycle_time);
    parallel.sum(projection_time);
    parallel.sum(snapshot_output_time);
    parallel.sum(spectra_output_time);
    parallel.sum(gravity_solver_time);
    parallel.sum(fft_time);
    parallel.sum(update_q_time);
    parallel.sum(moveParts_time);
    
    COUT << endl << "BENCHMARK" << endl;   
    COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
    COUT << "total number of cycles: "<< cycle << endl;
    COUT << "time consumption breakdown:" << endl;
    COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
    COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;
    
    COUT << "----------- main loop: components -----------"<<endl;
    COUT << "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
    COUT << "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
    COUT << "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
    COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time <<"%."<<endl;
    COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time <<"%."<<endl;
    COUT << "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
    COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
#endif

#ifdef EXTERNAL_IO    
        ioserver.stop();
    }
#endif

    return 0;
}

