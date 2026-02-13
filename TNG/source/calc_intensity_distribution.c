#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

#include "constant.h"
#include "para.h"
#include "proto.h"

#define nsnap 17
#define boxsize 205.0 //cMpc/h

#define nbin 50
#define fmin 0.001 //[mJy]
#define fmax 100 //[mJy]

#define xmax 60.0 //[arcmin]
#define xmin 0.0 //[arcmin]

void read_comdis_table(double redshift[], double *comdis_min, double *comdis_max);

int main(int argc, char *argv[])
{
	int i, j, k, ids[nsnap], ngal, itemp, imax;
	char path[256], fname[256], fbase[256];
	double redshift[Ntable], x[3], xtemp, ytemp, ztemp, temp, lambda[Nline], nu[Nline], nu_obs, numin, numax;
	double com_len, lumi_dis, com_dis, com_dis0, z0, comdis_min, comdis_max, dcdis, flux, dlogf, nbeam, num[Nline+1][nbin], num_cum;
	double h, lmd, omb, om;
	Cosmology cosm;
	FILE *fp;

	if(!(argc == 4))
	{
		fprintf(stderr, "please input numax, numin and beam size\n");
		fprintf(stderr, "band3: 84.176 GHz - 114.928 GHz\n");
		fprintf(stderr, "band6: 212.032 GHz - 272.011 GHz\n");
		exit(1);
	}
	numax = atof(argv[1]);
	numin = atof(argv[2]);
	nbeam = ( xmax - xmin ) * ( xmax - xmin ) / atof(argv[3]);

	fprintf(stderr, "observed frequency: %.2f - %.2f (GHz)\n", numin, numax);
	fprintf(stderr, "beam size: %.2f (arcmin^2)\n", atof(argv[3]));

	dlogf = log10( fmax / fmin ) / (double)nbin;

	h = 0.703;
	lmd = 0.729;
	om = 0.271;
	omb = 0.0451;
	init_cosm(&cosm, om, omb, lmd, h);

	sprintf(path, "/Users/kanamoriwaki/TNG/subgroup_data");
	sprintf(fbase, "subgroup_data.all.mJy_DL_DL_1+z");

	ids[0] = 33;
	ids[1] = 31;
	ids[2] = 29;
	ids[3] = 27;
	ids[4] = 25;
	ids[5] = 24;
	ids[6] = 23;
	ids[7] = 21;
	ids[8] = 19;
	ids[9] = 17;
	ids[10] = 15;
	ids[11] = 13;
	ids[12] = 11;
	ids[13] = 9;
	ids[14] = 7;
	ids[15] = 5;
	ids[16] = 3;
	
	read_line_lambda(lambda, nu);


	for(i = 0; i < nbin; i++)
	{
		for(j = 0; j < Nline+1; j++)
		{
			num[j][i] = 0.0;
		}
	}

	z0 = 2.0;
	com_dis0 = comoving_distance(cosm, z0) * h / mpc - 0.5 * boxsize; //cMpc/h
	read_comdis_table(redshift, &comdis_min, &comdis_max);
	dcdis = (comdis_max - comdis_min) / (double)Ntable;

	for(i = 0; i < nsnap; i++)
	{
		sprintf(fname, "%s/%s.%d.txt", path, fbase, ids[i]);
		fp = fopen(fname, "r");
		if(fp == NULL)
		{
			fprintf(stderr, "can not open %s\n", fname);
			exit(1);
		}

		fscanf(fp, "%*s %lf", &z0);
		fscanf(fp, "%*s %*f");
		fscanf(fp, "%*s %d", &ngal);

		fprintf(stderr, "%d: %.1f %d ", ids[i], z0, ngal);

		com_dis = boxsize * (double)(i+0.5) + com_dis0;//cMpc/h
		while( fscanf(fp, "%lf %lf %lf %*s", &x[0], &x[1], &x[2]) != EOF )//cMpc/h
		{
			com_dis = x[2] + boxsize*(double)i + com_dis0; //cMpc/h
			itemp = (int)( (com_dis - comdis_min) / dcdis );
			if(itemp < 0 || itemp > Ntable)
			{
				fprintf(stderr, "not enough Ntable\n");
				exit(1);
			}
			ztemp = redshift[ itemp ];
			lumi_dis = com_dis * ( 1.0 + ztemp ) * mpc / h; //cm

			com_len = com_dis * 1.0e6; //pc/h
			xtemp = x[0] * mpc / au / com_len / 60.0;//arcmin
			ytemp = x[1] * mpc / au / com_len / 60.0;//arcmin

			for(j = 0, imax = 0; j < Nline; j++)
			{
				nu_obs = nu[j] / (1.0 + ztemp) / 1.0e9; //GHz
				fscanf(fp, "%lf ", &temp);

				if( nu_obs > numin && nu_obs < numax 
						&& xtemp < xmax && xtemp > xmin 
						&& ytemp < xmax && ytemp > xmin )
				{
					flux = pow(10.0, temp) * (1.0+ztemp) / lumi_dis / lumi_dis;//mJy
					if( flux > fmin && flux < fmax )
					{
						itemp = (int)( log10( flux / fmin ) / dlogf );
						num[j][itemp] += 1.0;
						if( itemp > imax ) imax = itemp;	
					}	
				}
			}		

			num[Nline][imax] += 1.0;
		}
	
		fprintf(stderr, "%.0f\n", num[Nline][30]);
		fclose(fp);
	}

	for(i = 0; i < nbin; i++)
	{
		flux = 	fmin * pow(10.0, (i+0.5) * dlogf);
		printf("%e ", flux);

		for(j = 0; j < Nline+1; j++)
		{
			for(k = i, num_cum = 0.0; k < nbin; k++)
			{
				num_cum += num[j][k] / nbeam;
			}
			printf("%e ", num_cum);
		}
		printf("\n");
	}
	
	return 0;
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

