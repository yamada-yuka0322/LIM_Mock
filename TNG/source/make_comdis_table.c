#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "cosmology.c"

#define mpc              3.086e24        // megaparsec (cm) //

#define comdis_max 7200.0
#define comdis_min 100.0//mpc/h

int main(int argc, char **argv)
{
	int i;
	double z, dz, cdis;
	double om, omb, lmd, h, time, H, a;
	Cosmology cosm;

	int itemp;
	double dd = (comdis_max - comdis_min) / (double)Ntable;
	double comdis_temp, ztable[Ntable];

	h = 0.6774;
	lmd = 0.6911;
	om = 0.3089;
	omb = 0.0486;

	init_cosm(&cosm, om, omb, lmd, h);

	for(i = 0; i < Ntable; i++)
	{
		ztable[i] = -1.0;
	}

	printf("#N %d\n", Ntable);
	printf("#min %f\n", comdis_min);
	printf("#max %f\n", comdis_max);
 
	dz = 0.001;
	for(i = 0; i < Ntable/dz; i++)
	  {
	    z = 0.0 + (double)i * dz;

	    cdis = comoving_distance(cosm, z);
		itemp = (int)( (cdis/mpc*h - comdis_min) / dd );
		if( itemp > -1 && itemp < Ntable )
		{
			ztable[itemp] = z;
		}
	}

	for(i = 0; i < Ntable; i++)
	{
		comdis_temp = comdis_min + dd * (double)i;

		if(ztable[i] > 0)
		{	
	    	printf("%f %f\n"
				, comdis_temp
				, ztable[i]
				);
		}
		else if(i > 0 && i < Ntable-1)
		{
			fprintf(stderr, "ztable[%d] = %f\n", i, ztable[i]);
			printf("%f %f\n"
					, comdis_temp
					, ( ztable[i+1] + ztable[i-1] ) / 2.0
				  );
		}
    }
}

