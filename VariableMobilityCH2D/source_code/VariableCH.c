/*

    This code solves the variable mobility Cahn-Hilliard equation which results in 
    morphological evolution during Spinodal Decomposition of a binary alloy and 
    subsequent coarsening due to bulk-diffusion controlled dynamics and interface-diffusion 
    controlled dynamics. For more information regarding the equations used in the algorithm,
    please refer to the paper by: J Zhu, LQ Chen, J Shen, V Tikare - Physical Review E, 1999 - APS.
    DOI: 10.1103/PhysRevE.60.3564

    Copyright (C) 2020  Abhinav Roy, M.P. Gururajan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


*/

//defining the required input header files

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "fftw3.h"
#include<gsl/gsl_rng.h>
#include "../header_file/headers.h"

//#define INTERFACE 

/*The above comment is to be used for selective compilation. #define INTERFACE :- to implement interface-diffusion controlled dynamics.*/

int main(void)
{
  FILE *fr, *fw;
  int Nx, Ny;
  double dx,dy;
  double A,a,kappa;
  int time_step;
  double dt;
  double c_zero, c_noise;
  double average_comp;
 
  fftw_complex *comp;  
  gsl_rng * ran_num;
  const gsl_rng_type * Taus;
  int i1,i2;

  //Clearing the output directory in order to print fresh result from every new simulation run.

  (void) system("rm -rf output/*");

  //Opening the file simulation_data where all the physical constants and system information related to the simulation is stored.

  if ((fw = fopen("output/simulation_data","w")) == NULL)
  {
    printf("Unable to open output/simulation_data. Exiting.");
    exit(0);
  }
  else
  {
    fw = fopen("output/simulation_data","w");
  }
  #ifdef INTERFACE /*If INTERFACE macro is defined, then this part of the code will be selectively compiled in order
                    to produce results for interface-diffusion controlled dynamics.*/
  {
    if ((fr = fopen("input/constants_interface","r")) == NULL)
    {
      printf("Unable to open input/constants_interface. Exiting.");
      exit(0);
    }
    else
    {
      fr = fopen("input/constants_interface","r");
    }
    (void)fscanf(fr,"%le%le%le",&kappa,&A,&a);
    fclose(fr);
    (void)fprintf(fw,"kappa = %le\nA = %le\na = %le\n",kappa,A,a);
  }
  #else //if INTERFACE macro not defined, then by default the bulk-diffusion-controlled dynamics will be the product of the simulation.
  {
    if ((fr = fopen("input/constants_bulk","r")) == NULL)
    {
      printf("Unable to open input/constants_bulk. Exiting.");
      exit(0);
    }
    else
    {
      fr = fopen("input/constants_bulk","r");
    }
    (void)fscanf(fr,"%le%le%le",&kappa,&A,&a);
    (void)fclose(fr);
    fprintf(fw,"kappa = %le\nA = %le\na = %le\n",kappa,A,a);
  }
  #endif

  //Taking the input for time information.

  if ((fr = fopen("input/time_info","r")) == NULL)
  {
    printf("Unable to open input/time_info. Exiting.");
    exit(0);
  }
  else
  {
    fr = fopen("input/time_info","r");
  }
  (void)fscanf(fr,"%d%le",&time_step,&dt);
  (void)fclose(fr);
  fprintf(fw,"time_step = %d\ndt = %le\n",time_step,dt);

  //Taking the input for system information, i.e.,grid points and grid spacing.

  if ((fr = fopen("input/system_info","r")) == NULL)
  {
    printf("Unable to open input/system_info. Exiting.");
    exit(0);
  }
  else
  {
    fr = fopen("input/system_info","r");
  }
  (void)fscanf(fr,"%d%d%le%le",&Nx,&Ny,&dx,&dy);
  fclose(fr);
  (void)fprintf(fw,"Nx = %d\nNy = %d\n", Nx, Ny);
  (void)fprintf(fw,"dx = %le\ndy = %le\n", dx, dy);

  //Taking the input of compistion profile, i.e., nominal alloy composition and composition noise strength.

  if ((fr = fopen("input/composition_profile","r")) == NULL)
  {
    printf("Unable to open input/composition_profile. Exiting.");
    exit(0);
  }
  else
  {
    fr = fopen("input/composition_profile","r");
  }
  (void)fscanf(fr,"%le%le", &c_zero, &c_noise);
  (void)fclose(fr);
  fprintf(fw,"c_zero = %le\nc_noise = %le\n", c_zero, c_noise);

  //Closing the file pointer for writing the simulation data.

  (void)fclose(fw);
  fflush(fw);
  //Performing the input data test to make sure all the simulation data provided are within acceptable range.

  input_data_test(kappa, A, a, time_step, dt, Nx, Ny, dx, dy, c_zero, c_noise);

  /*Defining the complex composition variable of fftw_complex type.
  For more information, kindly refer to FFTW3 manual.*/

  comp = fftw_malloc(Nx*Ny* sizeof(fftw_complex));

  /*Defining variable for the Tausworthe uniform random number generator.
  For more information, kindly refer to the GNU Scientific Library (GSL) documentation.*/

  (void) gsl_rng_env_setup();
  Taus = gsl_rng_taus;
  ran_num = gsl_rng_alloc (Taus);

  //Setting the initial composition profile.

  average_comp = 0.0;
  for(i1=0; i1 < Nx; i1++)
  {
    for(i2=0; i2 < Ny; i2++)
    {
      __real__(comp[i2+Ny*i1]) = c_zero + c_noise*(0.5 - gsl_rng_uniform_pos(ran_num));
      __imag__(comp[i2+Ny*i1]) = 0.0;
      average_comp = average_comp + __real__(comp[i2+Ny*i1]);
    }
  }
  average_comp = average_comp/(Nx*Ny);

  //Writing the average composition as part of simulation data. Average composition is the actual alloy composition.

  if ((fw = fopen("output/simulation_data","a")) == NULL)
  {
    printf("Unable to open output/simulation_data. Exiting.");
    exit(0);
  }
  else
  {
    fw = fopen("output/simulation_data","a");
  }
  fprintf(fw,"average_comp = %le\n", average_comp);
  (void)fclose(fw);
  //Performing the evolution using Variable Mobility Cahn-Hilliard equation.

  cahn_hilliard_evolution(kappa, A, a, time_step, dt, Nx, Ny, dx, dy, comp);

  //free the composition and the random number variable.
  fftw_free(comp);
  gsl_rng_free(ran_num);

  return 0;
}
/*------------------------------------------------End of CODE---------------------------------------------------------*/
