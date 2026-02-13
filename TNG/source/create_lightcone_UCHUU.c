#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include <string.h>
#include <time.h>
#include "hdf5.h"

#include "constant.h"
#include "para.h"
#include "proto.h"

#define NSNAP_MAX 100
#define MAX_PARAMS 100
#define MAX_KEY_LEN 64

typedef struct {
    char key[MAX_KEY_LEN];
    double value;
} Param;

typedef struct {
    Param params[MAX_PARAMS];
    int count;
} ParamDict;

void read_redshifts(int snapshot_ids[], float redshifts[], double r_box[], int *Nsnap);
void read_comdis_table(double redshift[], double *comdis_min, double *comdis_max);
double get_mstar(double log_mass, double redshift, double smass_params[]);
int get_smass_params(const ParamDict *dict, double smass_params[]);
double get_SFR(double log_mass, double log_smass, double redshift, double SFR_params[]);
int get_SFR_params(const ParamDict *dict, double SFR_params[]);
double fmin(double a, double b);
double fmax(double a, double b);

int main(int argc, char *argv[])
{
	srand(time(NULL));
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

	ParamDict smass_dict;
    load_params("../params/Girelli.par", &smass_dict);
	double smass_params[8];
	get_smass_params(smass_dict, smass_params);

	ParamDict SFR_dict;
	load_params("../params/SIDES.par", &SFR_dict);
	double SFR_params[24];
	get_SFR_params(SFR_dict, SFR_params);




	for(i = 0; i < Nsnap; i++)
	{
		char fname[256];
		FILE *fp;
		double z0, boxsize, com_dis_ref;
		hsize_t dims[1];
		int ngal;

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
		
        hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

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
        for (int k = 0; k < ngal; k++)
        {
            x[0] = _x[k];
            x[1] = _y[k];
            x[2] = _z[k];
            v[0] = _vx[k];
            v[1] = _vy[k];
            v[2] = _vz[k];
            rmass_half = r_half[k];
            log_mass = log10(_mass[k]/h);

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
			log_Z = -0.50;
			log_mstar = get_mstar(log_mass, znow, smass_params[]);
			log_sfr = get_SFR(log_mass, log_mstar, znow, SFR_params);
		
			fprintf(stdout, "%e %e %.4f %.4f %e %e %.4f %.4f %.4f\n", 
				x_arcsec[0], 
				x_arcsec[1], 
				znow, 
				log_lumi_dis, 
				v[2], 
				rmass_half, 
				log_sfr, 
				log_mstar, 
				log_Z
			);
			
			count ++;
        }
		
		fprintf(stderr, "%d / %d\n", count, ngal);

		fclose(fp);
	}

	return 0;
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

void load_params(const char *filename, ParamDict *dict) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open file");
        exit(1);
    }

    char line[256];
    dict->count = 0;

    while (fgets(line, sizeof(line), fp)) {
        // Skip comments and empty lines
        if (line[0] == '#' || strlen(line) < 3) continue;

        char key[MAX_KEY_LEN];
        double value;
        if (sscanf(line, "%s = %lf", key, &value) == 2) {
            strncpy(dict->params[dict->count].key, key, MAX_KEY_LEN);
            dict->params[dict->count].value = value;
            dict->count++;
        }
    }

    fclose(fp);
}

int get_param(const ParamDict *dict, const char *key, double *out_value) {
    for (int i = 0; i < dict->count; i++) {
        if (strcmp(dict->params[i].key, key) == 0) {
            *out_value = dict->params[i].value;
            return 1;
        }
    }
    return 0;
}

int get_smass_params(const ParamDict *dict, double smass_params[]){
	const char *keys[] = { "B", "mu", "C", "nu", "D", "eta", "F", "E" };
    const int nkeys = sizeof(keys) / sizeof(keys[0]);

    for (int i = 0; i < nkeys; i++) {
        if (get_param(dict, keys[i], &smass_params[i])) {
            printf("%s = %f\n", keys[i], smass_params[i]);
        } else {
            printf("%s not found\n", keys[i]);
            return 1;
        }
    }
    return 0;
}

int get_SFR_params(const ParamDict *dict, double SFR_params[]){
	const char *keys[] = { "alpha1",
		"alpha2",
    	"beta1",
    	"beta2",
    	"gamma",
    	"Mt0",
    	"qfrac0",
    	"sigma0",
    	"m0",
    	"a0",
    	"a1",
    	"m1",
    	"a2",
    	"sigma_MS",
    	"logBsb",
    	"logx0",
    	"Psb_hz",
    	"slope_Psb",
    	"z_Psb_knee",
    	"Chab2Salp",
    	"zmean_lowzcorr",
    	"corr_zmean_lowzcorr",
    	"zmax_lowzcorr",
		"SFR_max" };
    const int nkeys = sizeof(keys) / sizeof(keys[0]);

    for (int i = 0; i < nkeys; i++) {
        if (get_param(dict, keys[i], &SFR_params[i])) {
            printf("%s = %f\n", keys[i], SFR_params[i]);
        } else {
            printf("%s not found\n", keys[i]);
            return 1;
        }
    }
    return 0;
}

double get_mstar(double log_mass, double redshift, double smass_params[]){
	double B, mu, C, nu, D, eta, F, E;
	double log_MA, A, gamma, beta;

	B = smass_params[0];
	mu = smass_params[0];
	C = smass_params[0];
	nu = smass_params[0];
	D = smass_params[0];
	eta = smass_params[0];
	F = smass_params[0];
	E = smass_params[0];

	log_MA = B + redshift*mu;
	A = C * pow((1+redshift),nu);
	gamma = D * pow((1+redshift), eta);
	beta = F * redshift + E;
	double log_hA = log_mass - log_MA;
	double hA = pow(10.0, log_hA);

	double log_ratio = log10(2.0*A) -log10(pow(hA, -beta) + pow(hA, gamma));
	double log_smass = log_mass + log_ratio;
	return log_smass;
}

double get_SFR(double log_mass, double log_smass, double redshift, double SFR_params[]){
	//Parameters to derive the quenched fraction
	double alpha1 = SFR_params[0];
	double alpha2 = SFR_params[1];
	double beta1 = SFR_params[2];
	double beta2 = SFR_params[3];
	double gamma = SFR_params[4];
	double Mt0 = SFR_params[5];
	double qfrac0 = SFR_params[6];
	double sigma0 = SFR_params[7];

	//Parameters for the evolution of the main sequence fom Schreiber et al.
	double m0 = SFR_params[8];
	double a0 = SFR_params[9];
	double a1 = SFR_params[10];
	double m1 = SFR_params[11];
	double a2 = SFR_params[12];

	//Main sequence scatter parameters
	double sigma_MS = SFR_params[13];
	double logBsb = SFR_params[14];
	double logx0 = SFR_params[15];

	//Evolution of the fraction of starburst (based on Sargent+12)
	double Psb_hz = SFR_params[16];
	double slope_Psb = SFR_params[17];
	double z_Psb_knee = SFR_params[18];

	//conversion from Chabrier to Salpeter IMF, used because Schreiber is in Salpeter
	double Chab2Salp = SFR_params[19];

	//Correction at low-z of the Schreiber relation (see Bethermin+17)
	double zmean_lowzcorr = SFR_params[20];
	double corr_zmean_lowzcorr = SFR_params[21];
	double zmax_lowzcorr = SFR_params[22];

	//Maximum SFR allowed
	double SFR_max = SFR_params[23];

	//Quenched flag
	int qflag;
	double Mtz = Mt0 + alpha1 * redshift + alpha2 * pow(redshift, 2.0);
	double sigmaz = sigma0 + beta1 * redshift + beta2 * pow(redshift, 2.0);
	double qfrac0z = qfrac0 * pow((1.0 + redshift), gamma);

	double erf_arg = (log_smass - Mtz) / sigmaz;
    double erf_val = erf(erf_arg);
    double Prob_SF = (1. - qfrac0z) * 0.5 * (1. - erf_val);

    double Xuni = (double)rand() / (double)RAND_MAX;
    qflag = (Xuni > Prob_SF);
	if (qflag){
		return 0.0;
	}else{
		double m_SF = log_smass + log10(Chab2Salp) - 9.0;
		double r = log10(1 + redshift);
		double expr = fmax(m_SF - m1 - a2 * r, 0.0);
		double logSFRms_SF = m_SF - m0 + a0 * r - a1 * pow(expr, 2.0) - log10(Chab2Salp);
		logSFRms_SF += corr_zmean_lowzcorr * (zmax_lowzcorr - fmin(redshift, zmax_lowzcorr) / (zmax_lowzcorr - zmean_lowzcorr));

    	double Psb = Psb_hz + slope_Psb * (z_Psb_knee - fmin(redshift, z_Psb_knee));
		Xuni = (double)rand() / (double)RAND_MAX;
		int issb = (Xuni < Psb );

		Xuni = (double)rand() / (double)RAND_MAX;
		double SFR_SF = pow(10.0,  ( logSFRms_SF + sigma_MS * Xuni + logx0 + (double)issb * (logBsb - logx0) ));

		int too_high_SFRs = (SFR_SF > SFR_max);
		while(too_high_SFRs)
		{
			Xuni = (double)rand() / (double)RAND_MAX;
			SFR_SF = pow(10.0,  ( logSFRms_SF + sigma_MS * Xuni + logx0 + (double)issb * (logBsb - logx0) ));
			too_high_SFRs = (SFR_SF > SFR_max);
		}

		return SFR_SF;
	}

}

double fmax(double a, double b){
	if (a>b)
	{
		return a;
	}else{
		return b;
	}
}

double fmin(double a, double b){
	if (a<b)
	{
		return a;
	}else{
		return b;
	}
}
