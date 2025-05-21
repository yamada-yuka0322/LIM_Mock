#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include <string.h>
#include "hdf5.h"

#include "constant.h"
#include "para.h"
#include "proto.h"

#define NSNAP_MAX 100

void read_snapshot_ids(int snapshot_ids[], double r_box[], int *Nsnap);
void read_comdis_table(double redshift[], double *comdis_min, double *comdis_max);

int main(int argc, char *argv[])
{
	int i, j;

	if(!(argc == 9))
	{
		fprintf(stderr, "Please input path (including a part of the file name),side length of the observation area [arcsec], minimum and maximum redshifts, x0, y0, tan_LoS_x, and tan_LoS_y \n");

		exit(1);
	}

	char path[256];
	sprintf(path, "%s", argv[1]);
	double dtheta[2] = {atof(argv[2]), atof(argv[2])}; // arcsec
	double zmin = atof(argv[3]), zmax = atof(argv[4]);
	
	double x0[2] = {atof(argv[5]), atof(argv[6])}; // [cMpc/h]
	double tan_LoS[2] = {atof(argv[7]), atof(argv[8])}; //line-of-sight direction of the lightcone with respect to simulation coordinate

	fprintf(stderr, "# read data from %s.{snapshot_id}.txt\n", path);
	fprintf(stderr, "# area: %.1f arcsec x %.1f arcsec\n", dtheta[0], dtheta[1]);

	// initialization //
	double h, lmd, omb, om;
	Cosmology cosm;
	h = 0.6774;
	lmd = 0.6911;
	om = 0.3089;
	omb = 0.0486;
	init_cosm(&cosm, om, omb, lmd, h);

	int Nsnap, snapshot_ids[NSNAP_MAX];
    double redshifts[NSNAP_MAX];
	double r_box[NSNAP_MAX]; // This fraction of the boxsize is used.
	read_redshifts(redshifts, snapshot_ids, r_box, &Nsnap);
	
	// read comoving distance table
	double redshift[Ntable], comdis_min, comdis_max;
	read_comdis_table(redshift, &comdis_min, &comdis_max);
	double dcdis = (comdis_max - comdis_min) / (double)Ntable;

	// read data //
	fprintf(stderr, "# snapshot: redshift r_start r_end com_dis_start com_dis_end com_dis_ref ... ngal_in_FoV / ngal\n");
	printf("# redshift: %.1f - %.1f\n", zmin, zmax);
	printf("# tan_LoS: %.2f %.2f\n", tan_LoS[0], tan_LoS[1]);
	printf("# area: %.1f arcsec x %.1f arcsec\n", dtheta[0], dtheta[1]);
	printf("# theta_x[arcsec] theta_y[arcsec] redshift log_lumi_dis vz[km/s] rmass_half log_sfr log_mstar log_Z\n");

	double com_dis_start, com_dis_end, r_start;
    char name[NSNAP_MAX];

	for(i = 0; i < Nsnap; i++)
	{
		char fname[256];
		FILE *fp;

        sprintf(name, "%.2f", redshifts[i]);
        for (j = 0; j < strlen(name); j++)
        {
            if (name[j] == '.')
            {
                name[j] = 'p';
            }
        }
        sprintf(fname, "%s.MiniUchuu_halolist_%s.h5", path, name);
		//sprintf(fname, "%s.MiniUchuu_halolist_%d.h5", path, snapshot_ids[i]);
		
        hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

        hid_t dset = H5Dopen(file_id, "x", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_x = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _x);

        dset = H5Dopen(file_id, "y", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_y = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _y);

        dset = H5Dopen(file_id, "z", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_z = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _z);

        dset = H5Dopen(file_id, "vx", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_vx = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vx);

        dset = H5Dopen(file_id, "vy", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_vy = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vy);

        dset = H5Dopen(file_id, "vz", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_vz = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vz);

        dset = H5Dopen(file_id, "Mvir", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *_mass = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _mass);

        dset = H5Dopen(file_id, "Halfmass_Radius", H5P_DEFAULT);
        space = H5Dget_space(dset);
        H5Sget_simple_extent_dims(space, dims, NULL);
        double *r_half = malloc(dims[0] * sizeof(double));
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, r_half);

        H5Dclose(dset);
        H5Sclose(space);
        H5Fclose(file_id);

        z0 = redshifts[i];
        boxsize = 400.0;
        ngal = dims[0];


		if(i == 0)
		{
			com_dis_ref = comoving_distance(cosm, z0) / ( mpc / h );
			com_dis_start = com_dis_ref - 0.5 * r_box[0] * boxsize; //cMpc/h
			com_dis_end = com_dis_start + r_box[0] * boxsize; //cMpc/h
			r_start = 0.0;
		}
		else 
		{
			com_dis_ref = comoving_distance(cosm, z0) / ( mpc / h );
			com_dis_start = com_dis_end;
			com_dis_end += r_box[i] * boxsize; //cMpc/h
			r_start += r_box[i-1];
			r_start = r_start - floor(r_start);
		}

		fprintf(stderr, "%d: %.1f %.2f %.2f %.2f %.2f %.2f ... ", snapshot_ids[i], z0, r_start, r_start + r_box[i], com_dis_start, com_dis_end, com_dis_ref);

		// read galaxy data
		int count = 0;
		double x[3], v[3], rmass_half, log_mass, log_sfr, log_mstar, log_Z;
        for (int k = 0; k < ngal; j++)
        {
            x[0] = _x[k];
            x[1] = _y[k];
            x[2] = _z[k];
            v[0] = _vx[k];
            v[1] = _vy[k];
            v[2] = _vz[k];
            rmass_half = r_half[k];
            log_mass = log10(_mass[k]);

            x[0] -= x0[0]; 
			x[1] -= x0[1]; 
			x[2] -= r_start * boxsize;
			if( x[2] < 0 ) x[2] += boxsize; // periodic boundary condition
			if( r_box[i] < 1.0)
			{
				// Galaxies in [r_start, r_start + r_box] of the boxsize are used.
				if( x[2] > r_box[i] * boxsize ) continue;
			}
			x[2] += com_dis_start; 

			// rotate coordinate //
			// first in z-x plane, then in z-y plane 
			for( j = 0; j < 2; j++)
			{
				double s = sin( atan(tan_LoS[j]) ), c = cos( atan(tan_LoS[j]) );

				// replicate boxes //
				while( x[j] < x[2] * tan_LoS[j] ) x[j] += boxsize;

				// rotate //
				double x_tmp = x[j] * c - x[2] * s; 
				x[2] = x[j] * s + x[2] * c; 
				x[j] = x_tmp;
				v[2] = v[j] * s + v[2] * c; // km/s
			}
			
			// compute comoving distance //
			double com_dis = 0.0;
			for(j = 0; j < 3; j++) 
			{
				com_dis += x[j]*x[j]; // cMpc/h
			}
			com_dis = sqrt( com_dis );
			
			// convert x, y [cMpc] to x, y [arcsec] //
			double x_arcsec[2];
			double com_dis_pc = com_dis * 1e6 / h; //pc
			for(j = 0; j < 2; j++)
			{
				double x_au = x[j] * mpc / h / au; //au
				x_arcsec[j] = x_au / com_dis_pc; //arcsec
			}
			if( x_arcsec[0] > dtheta[0] ) continue;
			if( x_arcsec[1] > dtheta[1] ) continue;
			rmass_half *= mpc / h / au / com_dis_pc; //arcsec

			// convert z [cMpc] to redshift //
			double ctmp = ( com_dis - comdis_min ) / dcdis;
			int itmp = (int)( ctmp );
			if(itmp + 1 > Ntable)
			{
				itmp = Ntable - 2;
			}
			else if(itmp < 0)
			{
				itmp = 0;
			}
			ctmp -= itmp;
			double znow = redshift[itmp+1] * ctmp + redshift[itmp] * (1. - ctmp);
			if( znow < zmin || znow >= zmax ) continue;
			
			// print //
			double log_lumi_dis = log10( com_dis * ( 1.0 + znow ) * mpc / h) ; //cm
		
			fprintf(stdout, "%e %e %.4f %.4f %e %e %.4f %.4f %.4f\n", 
				x_arcsec[0], 
				x_arcsec[1], 
				znow, 
				log_lumi_dis, 
				v[2], 
				rmass_half, 
				log_mass
			);
			
			count ++;
        }
		
		fprintf(stderr, "%d / %d\n", count, ngal);

		fclose(fp);
	}

	return 0;
}

void read_snapshot_ids(int snapshot_ids[], double r_box[], int *Nsnap)
{
	FILE *fp;
	int i;

	fp = fopen("./snapshot_ids.txt", "r");
	if(fp == NULL)
	{
		fprintf(stderr, "can not open snapshot_ids.txt\n");
		exit(1);
	}

	fscanf(fp, "%*s %d", Nsnap);
	if( *Nsnap > NSNAP_MAX )
	{
		fprintf(stderr, "Error: Nsnap > NSNAP_MAX\n");
		exit(1);
	}

	for(i = 0; i < *Nsnap; i++)
	{
		fscanf(fp, "%d %lf", &snapshot_ids[i], &r_box[i]);
	}

	fclose(fp);
	fprintf(stderr, "# Loaded %d snapshot ids from ./snapshot_ids.txt\n", *Nsnap);
}

void read_redshifts(int snapshot_ids[], float redshifts[], double r_box[], int *Nsnap)
{
	FILE *fp;
	int i;

	fp = fopen("./snapshot_ids_UCHUU.txt", "r");
	if(fp == NULL)
	{
		fprintf(stderr, "can not open snapshot_ids_UCHUU.txt\n");
		exit(1);
	}

	fscanf(fp, "%*s %d", Nsnap);
	if( *Nsnap > NSNAP_MAX )
	{
		fprintf(stderr, "Error: Nsnap > NSNAP_MAX\n");
		exit(1);
	}

	for(i = 0; i < *Nsnap; i++)
	{
		fscanf(fp, "%d %lf %lf", &snapshot_ids[i], &r_box[i], &redshifts[i]);
	}

	fclose(fp);
	fprintf(stderr, "# Loaded %d snapshot ids from ./snapshot_ids_UCHUU.txt\n", *Nsnap);
}

void read_comdis_table(double redshift[], double *comdis_min, double *comdis_max)
{
	FILE *fp;
	int Ntemp, i;
	double temp;

	fp = fopen("cdis_table.txt", "r");
	if(fp == NULL)
	{
		fprintf(stderr, "can not open comdis_table.txt\n");
		exit(1);
	}

	fscanf(fp, "%*s %d", &Ntemp);
	if( Ntemp != Ntable )
	{
		fprintf(stderr, "Ntemp != Ntable\n");
		exit(1);
	}
	
	fscanf(fp, "%*s %lf", &temp);
	*comdis_min = temp;
	fscanf(fp, "%*s %lf", &temp);
	*comdis_max = temp;

	for(i = 0; i < Ntable; i++)
	{
		fscanf(fp, "%*s %lf", &redshift[i]);
	}
	
	fclose(fp);
}