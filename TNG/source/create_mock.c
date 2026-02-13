#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include "hdf5.h"

#include "constant.h"
#include "para.h"
#include "proto.h"

void save_data_cube(char *fout, float *flist, float *intensity, int Nx, int Ny, int Nz);
void check_image(char *fname, float *data_cube, int Nx, int Ny, int Nz);
int binary_search(float arr[], int size, float value);
void read_filter_table(char *fname, float fmin, float fmax, float *flist, float *dflist, int *Nz);
void make_2d_gaussian_kernel(float w[], float sigma, int nw);

int main(int argc, char *argv[])
{
	int i;
	int Nx, Nz;
	if(!(argc == 11))
	{
		printf("Please input path, model id, side length [arcsec], fmin [GHz], fmax [GHz], angular resolution [arcsec], spectral resolution (R = f/df or the name of filter table if USE_FILTER_TABLE), alpha, output file id, and threshold (If threshold > 0, create galaxy number density with flux > threshold; otherwise create intensity map)\n");
		exit(1);
	}

	char path[256];
	sprintf(path, "%s", argv[1]);
	int model_id = atoi(argv[2]);
	double theta[2] = {atof(argv[3]), atof(argv[3])}; // arcsec
	double fmin = atof(argv[4]) * 1e9; // Hz
	double fmax = atof(argv[5]) * 1e9; // Hz
	double dtheta[2] = {atof(argv[6]), atof(argv[6])}; // arcsec

	double threshold = atof(argv[10]);
	if( threshold < 0 )
	{
		printf("# create intensity map\n");
	}
	else
	{
		printf("# create galaxy number density with flux > %e [erg/s/cm2]\n", threshold);
	}

	Nx = (int)( theta[0] / dtheta[0] );

	printf("# frequency: %.1f - %.1f [GHz]\n", fmin*1e-9, fmax*1e-9); 
	printf("# area: %.1f arcsec x %.1f arcsec\n", theta[0], theta[1]);
	printf("# angular resolution: %.1f arcsec x %.1f arcsec \n", dtheta[0], dtheta[1]);

	// spectral resolution //
	float *flist, *dflist;
#if defined USE_FILTER_TABLE
	char fname_filter[256];
	sprintf(fname_filter, "%s", argv[7]);
	read_filter_table(fname_filter, fmin, fmax, flist, dflist, &Nz);
	printf("# spectral resolution (GHz): %.1f \n", dflist[0]*1e-9);
#else
	double R = atof(argv[7]); // f/df	
							  
	Nz = 0;
	float fnow = fmin;
	while( fnow < fmax ) fnow += fnow / R, Nz += 1;

	flist = (float *)malloc((Nz+1) * sizeof(float));
	dflist = (float *)malloc(Nz * sizeof(float));
	for(fnow = fmin, i = 0; i < Nz; i++)
	{
		flist[i] = fnow; 
		dflist[i] = fnow / R;
		fnow += fnow / R;
	}
	flist[Nz] = fnow;
	printf("# spectral resolution (R): %.1f \n", R);
#endif

	printf("# Nx: %d, Nz: %d\n", Nx, Nz);

	// Dust-rich universe model parameters //
	double z_middle = 8; // Redshift where the sfr_factor becomes half of the maximum
	double alpha = atof(argv[8]); // Maximum factor for sfr. At highest redshift, SFR' = alpha * SFR
	printf("# z_middle: %.1f, alpha: %.1f\n", z_middle, alpha);

	// output file id //
	char fout_id[256];
	sprintf(fout_id, "%s", argv[9]);
	printf("# output file id: %s\n", fout_id);
		
	// allocate memory //
	float *intensity, *intensity_line[Nline]; 
	intensity = (float *)malloc(Nx*Nx*Nz * sizeof(float));
	for(i = 0; i < Nline; i++) intensity_line[i] = (float *)malloc(Nx*Nx*Nz * sizeof(float));

	for (i = 0; i < Nx*Nx*Nz; i++) 
	{
		intensity[i] = 0.0;
		for(int j = 0; j < Nline; j++) intensity_line[j][i] = 0.0;
	}	
	
	// read lambda and nu of each line
	double lambda[Nline], nu[Nline];
	read_line_lambda(lambda, nu);

	// read data //
	printf("### reading lightcone data %s ...\n", path);
	char fname[256];
	FILE *fp;
	sprintf(fname, "%s", path);
	fp = fopen(fname, "r");
	if(fp == NULL)
	{
		fprintf(stderr, "can not open %s\n", fname);
		exit(1);
	}

	// read header data
	double zmin, zmax, theta_lc[2];
	fscanf(fp, "%*s %*s %lf %*s %lf", &zmin, &zmax);
	fscanf(fp, "%*s %*s %*s %*s");
	fscanf(fp, "%*s %*s %lf %*s %*s %lf %*s", &theta_lc[0], &theta_lc[1]);
	for(i = 0; i < 10; i++) fscanf(fp, "%*s");
	
	printf("# redshift: %.1f - %.1f\n", zmin, zmax);
	printf("# area of lightcone: %.1f arcsec x %.1f arcsec\n", theta_lc[0], theta_lc[1]);
	double x[2], redshift, redshift_obs, log_lumi_dis, vz, rmass_half, log_sfr, log_mstar, log_Z;
	int count = 0;

	while( fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &x[0], &x[1], &redshift, &log_lumi_dis, &vz, &rmass_half, &log_sfr, &log_mstar, &log_Z) != EOF )
	{
		int ix = (int)( x[0] / dtheta[0] ); 
		int iy = (int)( x[1] / dtheta[1] );
		redshift_obs = redshift + vz * 1e5 / cspeed * (1 + redshift); // redshift space; the factor 1e5 is for converting vz from [km/s] to [cm/s]
		
		// parameters for angular size of the galaxy 
		double irx = rmass_half / dtheta[0];
		int nw = (int)(irx * 3 + 1);
		float w[nw*nw];
		make_2d_gaussian_kernel(w, irx / 0.674, nw); // 50 % is within 0.674 sigma 
		
		// compute luminosity //
		double log_lumi[Nline];
		double sfr_factor = (alpha - 1) * (0.5 * ( tanh( 0.5 * (redshift - z_middle) ) + 1 )) + 1;
		log_sfr += log10(sfr_factor); // [Msun/yr]
		if( log_mstar < -10 ) log_mstar = -10;
		if( log_Z < -10 ) log_Z = -10;
		line_luminosity_CO(log_lumi, nu, pow(10.0, log_sfr), pow(10.0, log_mstar), pow(10.0, log_Z), model_id);

		// assign line fluxes to pixels //
		for(i = 0; i < Nline; i++)
		{
			double flux = pow(10.0, log_lumi[i] - 2 * log_lumi_dis) / ( 4. * M_PI ); // erg/s/cm2

			float nu_obs = nu[i] / (1.0 + redshift_obs); // [Hz]
			int iz = binary_search(flist, Nz, nu_obs);

			if( threshold < 0 ) // create intensity map
			{
				if( ix >= -nw && ix < Nx + nw && iy >= -nw && iy < Nx + nw && iz != -1)
				{
					int di[2];
					for(di[1] = -nw+1; di[1] < nw; di[1]++)
					{
						for(di[0] = -nw+1; di[0] < nw; di[0]++) 
						{
							int iw = abs(di[0]) + abs(di[1])*nw;

							int itmp = (ix+di[0])*Nx*Nz + (iy+di[1])*Nz + iz;
							float tmp = flux / dflist[iz] / Jy / (dtheta[0]*arcsec) / (dtheta[1]*arcsec) * w[iw]; // [Jy/sr]
							
							if( tmp > 0 )
							{
								intensity[itmp] += tmp; // [Jy/sr]
								intensity_line[i][itmp] += tmp; // [Jy/sr]
							}
						}
					}
				}
			}
			else // create number density map
			{
				if( flux > threshold && ix >= 0 && ix < Nx && iy >= 0 && iy < Nx && iz != -1)
				{
					//printf("%d %.4f %.4f %.4f \n", i, redshift, log_lumi[i], log10(flux));
				
					int itmp = ix*Nx*Nz + iy*Nz + iz;
					intensity[itmp] ++;
					intensity_line[i][itmp] ++;						
				}
			}
		}
		count ++;
#if defined DEBUG
	if(count > 10000) break;
#endif
	}
	printf("# number of galaxies: %d\n", count);
	fclose(fp);

	check_image("check_image.txt", intensity, Nx, Nx, Nz);

	char fout[256];
	sprintf(fout, "%s.h5", fout_id);
	save_data_cube(fout, flist, intensity, Nx, Nx, Nz);

	for(i = 0; i < Nline; i++)
	{
		sprintf(fout, "%s_%d.h5", fout_id, i);
		save_data_cube(fout, flist, intensity_line[i], Nx, Nx, Nz);
	}

	free(intensity);
	for(i = 0; i < Nline; i++) free(intensity_line[i]);
	free(flist);
	free(dflist);
	
	return 0;
}

void save_data_cube(char *fout, float *flist, float *intensity, int Nx, int Ny, int Nz)
{
	hid_t file_id = H5Fcreate(fout, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if( file_id < 0 )
	{
		fprintf(stderr, "can not create %s\n", fout);
		exit(1);
	}

	hid_t dataspace_id;
	hid_t dataset_id;
	
	dataspace_id = H5Screate_simple(1, (hsize_t[]){Nz+1}, NULL);
	if (dataspace_id < 0) {
        printf("Error creating dataspace\n");
        H5Fclose(file_id);
        exit(1);
    }
	dataset_id = H5Dcreate(file_id, "/frequency", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dataset_id < 0) {
        printf("Error creating dataset\n");
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        exit(1);
    }

	H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flist);

	H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

	dataspace_id = H5Screate_simple(3, (hsize_t[]){Nx, Ny, Nz}, NULL);
	if (dataspace_id < 0) {
        printf("Error creating dataspace\n");
        H5Fclose(file_id);
        exit(1);
    }
	dataset_id = H5Dcreate(file_id, "/intensity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dataset_id < 0) {
        printf("Error creating dataset\n");
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        exit(1);
    }

	H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intensity);

	H5Dclose(dataset_id);
    H5Sclose(dataspace_id);


    H5Fclose(file_id);

	/*
	long num[3] = {Nx, Ny, Nz};
	int status = 0;
	long fpixel[] = {1,1,1};
	int ntot = num[0] * num[1] * num[2];
	
	fitsfile *fptr;
	fits_create_file(&fptr, fout, &status);
	fits_create_img(fptr, -32, 3, num, &status);
	
	fits_write_pix(fptr, TFLOAT, fpixel, ntot, intensity, &status);
	
	fits_close_file(fptr, &status);
	*/

	printf("# output: %s\n", fout);
}

void check_image(char *fname, float *data_cube, int Nx, int Ny, int Nz)
{
	FILE *fp;
	fp = fopen(fname, "w");

	for(int i = 0; i < Nx; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			float tot = 0;
			for(int k = 0; k < Nz; k++)
			{
				tot += data_cube[i*Ny*Nz + j*Nz + k];
			}
			fprintf(fp, "%e ", tot);
		}
		fprintf(fp, "\n");
	}
	
	printf("# output: %s\n", fname);
}

// binary search
// return the index of the array where the value is located
int binary_search(float arr[], int size, float value)
{
	int left = 0;
	int right = size;

	if( value < arr[left] || value >= arr[right] )
	{
		return -1;
	}
	
	while (left <= right)
	{
		int mid = left + (right - left) / 2;
		if( value >= arr[mid] && value < arr[mid+1] )
		{
			return mid;
		}
		else if (value < arr[mid])
		{
			right = mid - 1;
		}
		else
		{
			left = mid + 1;
		}
	}
	
	return -1;

}

void read_filter_table(char *fname, float fmin, float fmax, float *flist, float *dflist, int *Nz)
{
	FILE *fp;
	float temp;

	fp = fopen(fname, "r");
	if( fp == NULL )
	{
		printf("cannot open %s\n", fname);
		exit(1);
	}

	while( fscanf(fp, "%f", &temp) != EOF ) 
	{
		if( temp > fmin && temp < fmax ) (*Nz) ++;
	}

	(*Nz) -= 1; // flist[Nz] < fmax

	fclose(fp);

	fp = fopen(fname, "r");
	
	flist = (float *)malloc((*Nz+1) * sizeof(float));
	dflist = (float *)malloc(*Nz * sizeof(float));
	int count = 0;
	while( fscanf(fp, "%f", &temp) != EOF ) 
	{
		if( temp > fmin && temp < fmax )
		{
			flist[count] = temp;
			count ++;
		}
	}

	for(int i = 0; i < *Nz; i++)
	{
		dflist[i] = flist[i+1] - flist[i];
	}

	fclose(fp);
	printf("# read filter table %s\n", fname);
}

double w_coef_2d(int i, int j)
{
	if( i*i+j*j == 0 ){ return 1.0; } // i=j=0
	else if( i*j == 0 ){ return 2.0; } // i=0 or i=j
	else{ return 4.0; } // i!=0 and j!=0
}

void make_2d_gaussian_kernel(float w[], float sigma, int nw)
{
	int i, j;
	float wtot = 0.0;

	for(i = 0; i < nw*nw; i++) w[i] = 0.0;

	for(j = 0; j < nw; j++)
	{
		for(i = 0; i < nw; i++)
		{
			double r2 = (double)( i*i + j*j );
			int itmp = i + j*nw;
			if( r2 <= (nw-1)*(nw-1) )
			{
				w[itmp] = exp( - r2 / (2.0 * sigma * sigma) );
				wtot += w_coef_2d(i,j) * w[itmp];
			}
			else
			{ 
				w[itmp] = 0.0; 
			}
		}
	}

	for(i = 0; i < nw*nw; i++) w[i] /= wtot;
}
