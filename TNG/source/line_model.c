#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

#include "constant.h"
#include "para.h"

#include "mt19937ar.c"
#include "mt19937ar.h"


void read_line_lambda(double lambda[], double nu[])
{
	int i;

	lambda[ID_CO_1_0] = 2601.7 * micron;
	lambda[ID_CO_2_1] = 1300.9 * micron;
	lambda[ID_CO_3_2] = 867.3 * micron;
	lambda[ID_CO_4_3] = 650.5 * micron;
	lambda[ID_CO_5_4] = 521.0 * micron;
	lambda[ID_CO_6_5] = 433.7 * micron;
	lambda[ID_CO_7_6] = 371.8 * micron;
	lambda[ID_CO_8_7] = 325.0 * micron;
	lambda[ID_CO_9_8] = 289.0 * micron;
	lambda[ID_CO_10_9] = 260.0 * micron;
	lambda[ID_CO_11_10] = 237.0 * micron;
	lambda[ID_CO_12_11] = 217.0 * micron;
	lambda[ID_CO_13_12] = 200.0 * micron;

	lambda[ID_CII] = 158.0 * micron;
	lambda[ID_OIII88] = 88.0 * micron;
	lambda[ID_NII205] = 205.0 * micron;
	lambda[ID_NII122] = 122.0 * micron;

	lambda[ID_CI_1_0] = 609.14 * micron;
	lambda[ID_CI_2_1] = 370.42 * micron;

	lambda[ID_OIII5007] = 5007.0 * angstrom;

	for(i = 0; i < Nline; i++)
	{
		nu[i] = cspeed / lambda[i];
	}
}

void line_luminosity_CO(double log_lumi[], double nu[], double sfr, double mstar, double metallicity, int model_id)
{
	double Z = metallicity;
	double ssfr = sfr / mstar; // [Gyr^-1]
	double logL_FIR, temp;
	
	//logL_FIR = log10( sfr / SFR_FIR / Lsun ); // [Lsun]
	logL_FIR = log10( sfr * 2.22e43 / Lsun ); //[Lsun]: See eq. 23 of Fonseca+2017
	
	// CO(1-0)
	log_lumi[ID_CO_1_0] = 0.81 * logL_FIR + 0.54; // Sargent+ 14 [K km s-1 pc2]
	if( ssfr > 0.5e-8 )  // starbursting galaxies. Criteria is from Fig. 1 of Sargent+14
	{
		log_lumi[ID_CO_1_0] -= 0.46; // See discussion after Eq. 5 of Bethermin+22 and Sargent+14
	}
			
	// CO(2-1)
	log_lumi[ID_CO_2_1] = log10( 0.76 ) + log_lumi[ID_CO_1_0]; // CO ratio from Daddi+15 [K km s-1 pc2]
			
	// CO(3-2)
	//log_lumi[ID_CO_3_2] = log10( 1.0e8 * sfr );//Papadopoulos 12
	//log_lumi[ID_CO_3_2] = log10( 3.2e8 * sfr ); //Popping+ 18 average [K km s-1 pc2]
	log_lumi[ID_CO_3_2] = log10( 0.42 ) + log_lumi[ID_CO_1_0]; // CO ratio from Daddi+15 [K km s-1 pc2]
			
	// CO(4-3)
	log_lumi[ID_CO_4_3] = ( logL_FIR - 1.49 ) / 1.06; // Liu+15 [K km s-1 pc2] 
			
	// CO(5-4)
	log_lumi[ID_CO_5_4] = ( logL_FIR - 1.71 ) / 1.07; // Liu+15 [K km s-1 pc2]
			
	// CO(6-5)
	log_lumi[ID_CO_6_5] = ( logL_FIR - 1.79 ) / 1.10; // Liu+15 [K km s-1 pc2]		

	// CO(7-6)
	log_lumi[ID_CO_7_6] = ( logL_FIR - 2.62 ) / 1.03; // Liu+15 [K km s-1 pc2]		

	// CO(8-7)
	log_lumi[ID_CO_8_7] = ( logL_FIR - 2.82 ) / 1.02; // Liu+15 [K km s-1 pc2]		

	// CO(9-8)
	log_lumi[ID_CO_9_8] = ( logL_FIR - 3.10 ) / 1.01; // Liu+15 [K km s-1 pc2]		

	// CO(10-9)
	log_lumi[ID_CO_10_9] = ( logL_FIR - 3.67 ) / 0.96; // Liu+15 [K km s-1 pc2]		

	// CO(11-10)
	log_lumi[ID_CO_11_10] = ( logL_FIR - 3.51 ) / 1.00; // Liu+15 [K km s-1 pc2]		

	// CO(12-11)
	log_lumi[ID_CO_12_11] = ( logL_FIR - 3.83 ) / 0.99; // Liu+15 [K km s-1 pc2]		

	// CO(13-12)
	log_lumi[ID_CO_13_12] = logL_FIR - 5.;

	for(int i = 0; i <= ID_CO_13_12; i++)
	{
		log_lumi[i] += log10( 1.227e-4 ) + 3.0 * log10( nu[i] ); //see eq. 28 of Fonseca+2017 or Carilli+2013 for this conversion.
		
		//double lambda_km = cspeed / nu[i] / km;
		//log_lumi[i] += log10( ( 2. * k_B * nu[i] * nu[i] / cspeed / cspeed ) / lambda_km * pc * pc ); 
	}// [erg/s]: 

	// [CII]158
	if( (model_id % 10) == 0 )
	{
		log_lumi[ID_CII] = (6.99 + log10(sfr) ) / 1.01 + log10( Lsun ); //DeLooze+14 [erg/s]
	}
	else if( (model_id % 10) == 1 )
	{
		if( Z < 1e-10 ) Z = 1e-10;

		temp = 7.0 + 1.2 * log10(sfr) + 0.021 * log10(Z/Zsun) + 0.012 * log10(sfr) * log10(Z/Zsun) - 0.74 * log10(Z/Zsun) * log10(Z/Zsun);
		log_lumi[ID_CII] = temp + log10( Lsun );// Vallini+ 2015 [erg/s]
	}

	// [OIII]88
	log_lumi[ID_OIII88] = ( 7.48 + log10(sfr) ) / 1.12 + log10( Lsun ); //DeLooze+14 [erg/s]
			
	// [NII]205
	log_lumi[ID_NII205] = log10( 2.5 * 1.0e5 * sfr * Lsun ); //Visbal Loeb [erg/s]

	// [NII]122
	log_lumi[ID_NII122] = log10( 7.9 * 1.0e5 * sfr * Lsun ); //Visbal Loeb [erg/s]

	// [CI](1-0)
	log_lumi[ID_CI_1_0] = 1.07 * ( log_lumi[ID_CO_4_3] - log10( Lsun ) - logL_FIR ) + 0.14 + logL_FIR + log10( Lsun ); // Bethermin+22 Eq. 9 [erg/s]
			
	// [CI](2-1)
	log_lumi[ID_CI_2_1] = log_lumi[ID_CI_1_0] + 0.63 * ( log_lumi[ID_CO_7_6] - log_lumi[ID_CO_4_3] ) + 0.17; // Bethermin+22 Eq. 10 [erg/s]
	double r1 = genrand_real3();
	double r2 = genrand_real3();
	log_lumi[ID_CI_2_1] += 0.19 * sqrt( -2.0 * log(r1) ) * sin( 2.0 * M_PI * r2 );

	// [OIII]5007
	log_lumi[ID_OIII5007] = log10( 1.32e41 ) + log10( sfr ); // [erg/s]

	// add scatter
	double sigma = 0.2;
	if( (int)(model_id / 10) == 1)
	{
		sigma = 0.4;
	}

	for(int i = 0; i < Nline; i++)
	{
		double r1 = genrand_real3();
		double r2 = genrand_real3();
		log_lumi[i] += sigma * sqrt( -2.0 * log(r1) ) * sin( 2.0 * M_PI * r2 );
	}
}
