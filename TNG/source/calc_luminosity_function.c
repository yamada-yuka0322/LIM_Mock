#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Lsun 3.826e33
#define h 0.6711
#define LF_BIN 40
#define Nline 16
#define log_lumi_min 7

int main(int argc, char *argv[])
{
	int i, j, k, num[Nline][LF_BIN]={}, cumnum[Nline][LF_BIN], cum_temp, ngalaxy;
	float *log_lumi[Nline], dlogl, logl_start, box3, lf[LF_BIN][1+Nline], logmax, logmin, z, boxsize;

	FILE *fp;
	char fname[256];
	double temp;

// set up
	


// read data

	if(argc < 2)
	{
		printf("please input file name\n");
		exit(1);
	}
	sprintf(fname, "%s", argv[1]);

	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		fprintf(stderr, "can not open input file\n");
		exit(1);
	}
		
	fscanf(fp, "%*s %f", &z);
	fscanf(fp, "%*s %f", &boxsize);
	fscanf(fp, "%*s %d", &ngalaxy);

//	box3 = boxsize * boxsize * boxsize / h / h / h;
    box3 = boxsize * boxsize * boxsize;

	fprintf(stderr, "ngalaxy = %d\n", ngalaxy);
	fprintf(stderr, "boxsize = %f\n", boxsize);
	fprintf(stderr, "redshfit = %f\n", z);
	for(i=0;i<Nline;i++)
	{
		log_lumi[i] = malloc(sizeof(float) * ngalaxy);
	}
	if( log_lumi[0] == NULL )
	{
		printf("cannot allocate log_lumi\n");
		exit(1);
	}


	for(i=0;i<ngalaxy;i++)
	{
		for(j = 0; j < 4; j++) fscanf(fp, "%*s ");
		for(j = 0; j < Nline; j++)
		{
			fscanf(fp, "%lf ", &temp);
			log_lumi[j][i] = temp;
//			log_lumi[j][i] = log10( temp );
		}
	}

	fclose(fp);

// lumi_max
	
	logmax = -100.0;
	logmin = 100.0;
	
	/*
	for(j=0;j<Nline;j++){
	for(i=0;i<ngalaxy;i++)
	{
		if( log_lumi[j][i] >= log_lumi_min )
		{
			if( logmax < log_lumi[j][i] )
			{
				logmax = log_lumi[j][i];
			}

			if( logmin > log_lumi[j][i] )
			{
				logmin = log_lumi[j][i];
			}
		}
	}}*/


	logmax = 12.0;
	logmin = 7.0;
	dlogl = (logmax - logmin) / (double) LF_BIN;
	logl_start = logmin + 0.5 * dlogl;

// calculate LF
	
	for(j=0;j<Nline;j++)
	{
		for(i=0;i<ngalaxy;i++)
		{
			if( log_lumi[j][i] >= log_lumi_min)
			{
				k = (int) ((log_lumi[j][i] - logmin) / dlogl);
				if( k<LF_BIN && k>=0 ) num[j][k]++;
			}
		}

		for(i=0, cum_temp = 0 ; i < LF_BIN ; i++)
		{
			cum_temp += num[j][LF_BIN - 1 - i];
			cumnum[j][LF_BIN - 1 - i] = cum_temp;
		}
	}
	
	for(k=0;k<LF_BIN;k++)
	{
		lf[k][0] = logl_start + dlogl * (double) k;
		printf("%f ", lf[k][0]);

		for(j=0;j<Nline;j++)
		{
			//lf[k][1+j] = (double) num[j][k] / box3 / dlogl;
			lf[k][1+j] = (double) cumnum[j][k] / box3;
			printf("%f ", log10( lf[k][1+j] ));
		}
		printf("\n");
	}

	for(j=0;j<Nline;j++)
	{
		free(log_lumi[j]);
	}

	return 0;

}

