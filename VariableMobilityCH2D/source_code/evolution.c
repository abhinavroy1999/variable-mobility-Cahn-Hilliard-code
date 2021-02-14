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
	double k2, k4;
	double denominator;
	int i1, i2;
	int temp = 0;
	fftw_complex *g;
	fftw_complex *r;
	fftw_complex *y1, *y2;
	fftw_complex *qx, *qy;
	fftw_complex *comp2, *xcomp;

	// Stability factor
	double alpha = 0.5;
	
	//double pointer variable for writing the composition profile to a file.
	c = (double *)malloc((size_t) Nx*Ny* sizeof(double));
	
	//allocating memory for the fftw_complex type variables.
	g = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	y1 = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	y2 = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	r = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	qx = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	qy = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	comp2 = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	xcomp = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
	
	//Writing the initial composition profile into an output file and generating the PS file as well.
	for (i1 = 0; i1<Nx; ++i1)
	{
		for (i2 = 0; i2<Ny; ++i2)
		{
			c[i2 + i1*Ny] = __real__(comp[i2 + i1*Ny]);
		}
	}
	for (i1 = 0; i1<Nx; ++i1)
	{
		for (i2 = 0; i2<Ny; ++i2)
		{
			xcomp[i2 + i1*Ny] = __real__(comp[i2 + i1*Ny]);
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
	fftw_plan plan1, plan2, plan3, plan4, plan5, plan6, plan7;

	plan1 = fftw_plan_dft_2d(Nx, Ny, comp, comp, FFTW_FORWARD, FFTW_ESTIMATE);
	plan2 = fftw_plan_dft_2d(Nx, Ny, g, g, FFTW_FORWARD, FFTW_ESTIMATE);
	plan3 = fftw_plan_dft_2d(Nx, Ny, y1, y1, FFTW_BACKWARD, FFTW_ESTIMATE);
	plan4 = fftw_plan_dft_2d(Nx, Ny, y2, y2, FFTW_BACKWARD, FFTW_ESTIMATE);
	plan5 = fftw_plan_dft_2d(Nx, Ny, qx, qx, FFTW_FORWARD, FFTW_ESTIMATE);
	plan6 = fftw_plan_dft_2d(Nx, Ny, qy, qy, FFTW_FORWARD, FFTW_ESTIMATE);
	plan7 = fftw_plan_dft_2d(Nx, Ny, comp, comp, FFTW_BACKWARD, FFTW_ESTIMATE);

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
				__real__(comp2[i2 + i1*Ny]) = __real__(comp[i2 + i1*Ny]);
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
		
		for(i1=0; i1 < Nx; ++i1)
		{
			if(i1 <= halfNx)
			{ kx = i1*delkx; }
			else
			{ kx = (i1-Nx)*delkx; }
			for(i2=0; i2 < Ny; ++i2)
			{
				if(i2 <= halfNy)
				{ ky = i2*delky; }
				else
				{ ky = (i2-Ny)*delky; }
				k2 = kx*kx + ky*ky;
				r[i2 + Ny*i1] = g[i2 + Ny*i1] + 2.0*kappa*k2*comp[i2 + Ny*i1];
			}
		}
		for(i1=0; i1 < Nx; ++i1)
		{
			if(i1 <= halfNx)
			{ kx = i1*delkx; }
			else
			{ kx = (i1-Nx)*delkx; }
			for(i2=0; i2 < Ny; ++i2)
			{
				if(i2 <= halfNy)
				{ ky = i2*delky; }
				else
				{ ky = (i2-Ny)*delky; }
				y1[i2 + Ny*i1] = (_Complex_I)*kx*r[i2 + Ny*i1];
				y2[i2 + Ny*i1] = (_Complex_I)*ky*r[i2 + Ny*i1];
			}
		}

		fftw_execute(plan3);	// inverse fourier transform of y1
		fftw_execute(plan4);	// inverse fourier transform of y2
		
		//Normalizing the r complex array after the inverse fourier transform (as required by the fftw algorithm).
		for(i1=0; i1<Nx; ++i1)
		{
			for(i2=0; i2<Ny; ++i2)
			{
				y1[i2+Ny*i1] = y1[i2+Ny*i1]/(Nx*Ny);
			}
		}
		for(i1=0; i1<Nx; ++i1)
		{
			for(i2=0; i2<Ny; ++i2)
			{
				y2[i2+Ny*i1] = y2[i2+Ny*i1]/(Nx*Ny);
			}
		}
		for (i1=0; i1 < Nx; ++i1)
		{
			for (i2=0; i2 < Ny; ++i2)
			{
				qx[i2 + Ny*i1] = (1.0 - a*xcomp[i2 + Ny*i1]*xcomp[i2 + Ny*i1])*y1[i2 + Ny*i1];
				qy[i2 + Ny*i1] = (1.0 - a*xcomp[i2 + Ny*i1]*xcomp[i2 + Ny*i1])*y2[i2 + Ny*i1];
			}
		}
		
		fftw_execute(plan5);//Fourier transform of q.
		fftw_execute(plan6);

		//evolving the composition profile.
		for(i1=0; i1 < Nx; ++i1)
		{
			if(i1 <= halfNx)
			{ kx = i1*delkx; }
			else
			{ kx = (i1-Nx)*delkx; }
			for(i2=0; i2 < Ny; ++i2)
			{
				if(i2 <= halfNy)
				{ ky = i2*delky; }
				else
				{ ky = (i2-Ny)*delky; }
				k2 = kx*kx + ky*ky;
				k4 = k2*k2;
				denominator = (1.0 + 2.0*alpha*dt*kappa*k4);
				comp[i2+Ny*i1] = comp[i2+Ny*i1] + (dt*((_Complex_I)*kx*qx[i2 + Ny*i1] + (_Complex_I)*ky*qy[i2 + Ny*i1]))/denominator;
			}
		}
		fftw_execute(plan7);
		
		for(i1=0; i1<Nx; ++i1)
		{
			for(i2=0; i2<Ny; ++i2)
			{
				comp[i2+Ny*i1] = comp[i2+Ny*i1]/(Nx*Ny);
			}
		}

		//Taking the output every 100 units interval of time.
		if (temp%100 == 0)
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
	fftw_free(y1);
	fftw_free(y2);
	fftw_free(qx);
	fftw_free(qy);
	fftw_free(comp2);
	fftw_free(xcomp);
	free(c);
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
	fftw_destroy_plan(plan3);
	fftw_destroy_plan(plan4);
	fftw_destroy_plan(plan5);
	fftw_destroy_plan(plan6);
	fftw_destroy_plan(plan7);
  
}
/*-------------------------------------------End of CODE------------------------------------------------------------------*/
