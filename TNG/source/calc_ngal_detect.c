#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

#define Nline 17
#define lmax 60.0

int main(int argc, char *argv[])
{
	int i, j, ngal_tot, itemp, *num, Nmesh, ix[2];
	double flux[Nline], Vsurv, x[3], dx;

	FILE *fp;
	char fname[256];

	if(!(argc==3))
	{
		fprintf(stderr, "please input file name and size of observation area\n");
		exit(1);
	}
	sprintf(fname, "%s", argv[1]);
	dx = sqrt(atof(argv[2]));

	fprintf(stderr, "dx = %.2f\n", dx);
	Nmesh = (int)( lmax / dx );

	num = malloc(sizeof(int) * Nmesh * Nmesh);
	if(num == NULL)
	{
		fprintf(stderr, "can not allocate\n");
		exit(1);
	}
	for(i = 0; i < Nmesh*Nmesh; i++)
	{
		num[i] = 0;
	}

	fprintf(stderr, "%s\n", fname);
	fp = fopen(fname, "r");
	if(fp == NULL)
	{
		fprintf(stderr, "can not open input file\n");
		exit(1);
	}

	ngal_tot = 0;
	while( fscanf(fp, "%lf %lf %lf %*f ", &x[0], &x[1], &x[2]) != EOF )
	{
		for(j = 0, itemp = 0; j < Nline; j++)
		{
			fscanf(fp, "%lf ", &flux[j]);
			if( flux[j] > 0 )
			{
				itemp ++;
			}
		}
		
		if( itemp > 0 && x[0] < lmax && x[1] < lmax )//&& x[2] > 5 )
		{
			ix[0] = (int)( x[0] / dx );
			ix[1] = (int)( x[1] / dx );
			num[ix[0] + ix[1]*Nmesh] ++;
			ngal_tot ++;
		}
	}

	fclose(fp);

	for(i = 0; i < Nmesh*Nmesh; i++)
	{
		printf("%d %d\n", i, num[i]);
	}
	fprintf(stderr, "ngal_tot = %d\n", ngal_tot);
	fprintf(stderr, "ngal_mean = %.3f\n", (double)ngal_tot / (double)(Nmesh*Nmesh));


	return 0;
}
