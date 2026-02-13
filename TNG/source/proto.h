#ifndef __PROTO__
#define __PROTO__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
//#include <hdf5.h>
#include "constant.h"
#include "para.h"

extern void read_line_lambda(double lambda[], double nu[]);
extern void line_luminosity_CO(double log_lumi[], double nu[], double sfr, double mstar, double metallicity, int model_id);

// cosmological calculation : cosmology.c //
extern void init_cosm(Cosmology *cosm, double om, double omb, double lmd, double h);
extern double timetoz(double tnow, Cosmology cosm);
extern double timetoa(double tnow, Cosmology cosm);
extern double ztotime(double znow, Cosmology cosm);
extern double atotime(double anow, Cosmology cosm);
extern double luminosity_distance(Cosmology cosm,double z);
extern double angular_distance(Cosmology cosm,double z);
extern double comoving_distance(Cosmology cosm,double z);

#endif

