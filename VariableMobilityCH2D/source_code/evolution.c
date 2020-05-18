/*

    This code performs evolution of the composition profile guided by 
    the variable mobility Cahn-Hilliard equaiton.

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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "fftw3.h"
#include<gsl/gsl_math.h>
#include "../header_file/headers.h"

void cahn_hilliard_evolution(double kappa, double A, double a, double time_step, double dt, int Nx, int Ny, double dx, double dy, fftw_complex *comp)
{
  FILE *fw;
  char file_name[50];
  char ps_file_name[50];
  double *c;
  int halfNx, halfNy;
  float kx, ky, delkx, delky;
  double k2, k4, k;
  double denominator;
  int i1, i2;
  int temp = 0;
  fftw_complex *g;
  fftw_complex *r;
  fftw_complex *q;
  fftw_complex *comp2;

  //double pointer variable for writing the composition profile to a file.
  c = (double *)malloc((size_t) Nx*Ny* sizeof(double));
  //allocating memory for the fftw_complex type variables.
  g = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  r = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  q = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  comp2 = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  //Writing the initial composition profile into an output file and generating the PS file as well.
  for (i1 = 0; i1<Nx; ++i1)
  {
    for (i2 = 0; i2<Ny; ++i2)
    {
      c[i2 + i1*Ny] = __real__(comp[i2 + i1*Ny]);
    }
  }
  sprintf(file_name, "output/time%d.dat", temp);
  fw = fopen(file_name, "wb");
  fwrite(&c[0], sizeof(double),(size_t) Nx*Ny,fw);
  (void) fclose (fw);
  fflush(fw);
  sprintf(ps_file_name,"output/time%d.ps",temp);
  ps_file(file_name, ps_file_name, Nx, Ny, Nx);

  //Defining the plans for fourier transforms.
  fftw_plan plan1, plan2, plan3, plan4, plan5;

  plan1 = fftw_plan_dft_2d(Nx, Ny, comp, comp, FFTW_FORWARD, FFTW_ESTIMATE);
  plan2 = fftw_plan_dft_2d(Nx, Ny, g, g, FFTW_FORWARD, FFTW_ESTIMATE);
  plan3 = fftw_plan_dft_2d(Nx, Ny, r, r, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan4 = fftw_plan_dft_2d(Nx, Ny, q, q, FFTW_FORWARD, FFTW_ESTIMATE);
  plan5 = fftw_plan_dft_2d(Nx, Ny, comp, comp, FFTW_BACKWARD, FFTW_ESTIMATE);


 
  //half of the grid points required for Periodic Boundary Condition (PBC) to get rid of the surface effects during simulation.
  halfNx = (int) Nx/2;
  halfNy = (int) Ny/2;

  //delta kx and delta ky - for defining the Fourier space vectors.
  delkx = (2*M_PI)/(Nx*dx);
  delky = (2*M_PI)/(Ny*dy);

  //Loop for temporal evolution.
  for (temp = 1; temp <time_step+1; ++temp)
  {
  for (i1 = 0; i1<Nx; ++i1)
  {
    for (i2 = 0; i2<Ny; ++i2)
    {
      __real__(comp2[i2 + i1*Ny]) = __real(comp[i2 + i1*Ny]);
      __imag__(comp2[i2 + i1*Ny]) = 0.0;

    }
  }
    for(i1=0; i1<Nx; ++i1)
    {
      for(i2=0; i2<Ny; ++i2)
      {
        g[i2+Ny*i1] = 2.0*A*(comp2[i2 + i1*Ny])*(1.0-(comp2[i2 + i1*Ny]))*(1.0-2.0*(comp2[i2 + i1*Ny])); //the derivative of the bulk free energy density.
      }
    }

    fftw_execute(plan1);
    fftw_execute(plan2);
    //Loop for implementing the PBC and subsequent evolution of the composition profile.
    for(i1=0; i1 < Nx; ++i1)
    {
	     if(i1 < halfNx)
       { kx = i1*delkx; }
	      else
       { kx = (i1-Nx)*delkx; }
       for(i2=0; i2 < Ny; ++i2)
       {
	        if(i2 < halfNy)
          { ky = i2*delky; }
	        else
          { ky = (i2-Ny)*delky; }
	         k2 = kx*kx + ky*ky;
	         k4 = k2*k2;
                 k = sqrt(k2);
          r[i2 + Ny*i1] = (0.0+1.0*I)*k*(g[i2 + Ny*i1] + kappa*k2*comp[i2 + Ny*i1]);
        }
    }
    fftw_execute(plan3); //inverse fourier transform of r.
    //Normalizing the r complex array after the inverse fourier transform (as required by the fftw algorithm).
    for(i1=0; i1<Nx; ++i1)
    {
          for(i2=0; i2<Ny; ++i2)
          {
           r[i2+Ny*i1] = r[i2+Ny*i1]/(Nx*Ny);
          }
    }
    for (i1=0; i1 < Nx; ++i1)
    {
      for (i2=0; i2 < Ny; ++i2)
      {
        q[i2 + Ny*i1] = (1.0 - a*(comp2[i2 + i1*Ny])*(comp2[i2 + i1*Ny]))*r[i2 + Ny*i1];
      }
    }
    fftw_execute(plan4);//Fourier transform of q.
    //evolving the composition profile.
    for(i1=0; i1 < Nx; ++i1)
    {
	     if(i1 < halfNx)
       { kx = i1*delkx; }
	      else
       { kx = (i1-Nx)*delkx; }
       for(i2=0; i2 < Ny; ++i2)
       {
	        if(i2 < halfNy)
          { ky = i2*delky; }
	        else
          { ky = (i2-Ny)*delky; }
	         k2 = kx*kx + ky*ky;
	         k4 = k2*k2;
          k = sqrt(k2);
          denominator = (1.0 + 2.0*0.5*dt*kappa*k4);
          comp[i2 + Ny*i1] = comp[i2 + Ny*i1] + ((0.0+1.0*I)*dt*k*q[i2 + Ny*i1])/denominator;
        }
    }
    fftw_execute(plan5);
    //Normalizing the composition.
    for(i1=0; i1<Nx; ++i1)
    {
          for(i2=0; i2<Ny; ++i2)
          {
           comp[i2+Ny*i1] = comp[i2+Ny*i1]/(Nx*Ny);
          }
    }

    //Taking the output every 1000 units interval of time.
      if (temp%1000 == 0)
      {
        for (i1 = 0; i1<Nx; ++i1)
        {
          for (i2 = 0; i2<Ny; ++i2)
          {
            c[i2 + i1*Ny] = __real__(comp[i2 + i1*Ny]);
          }
        }
        sprintf(file_name,"output/time%d.dat",temp);
        fw = fopen(file_name,"wb");
        fwrite(&c[0], sizeof(double),(size_t)Nx*Ny, fw);
	sprintf(ps_file_name,"output/time%d.ps",temp);
        ps_file(file_name, ps_file_name, Nx, Ny, Nx);
        (void) fclose (fw);
	fflush(fw);
      }
  }
  //free the memory allocated for all the variables.
  fftw_free(g);
  fftw_free(r);
  fftw_free(q);
  fftw_free(comp2);
  free(c);
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  fftw_destroy_plan(plan3);
  fftw_destroy_plan(plan4);
  fftw_destroy_plan(plan5);
  
}
/*-------------------------------------------End of CODE------------------------------------------------------------------*/
