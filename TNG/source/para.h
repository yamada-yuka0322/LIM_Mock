#ifndef __PARAM__
#define __PARAM__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#define Ntable 800

#define Nline 20

#define ID_CO_1_0 0
#define ID_CO_2_1 1
#define ID_CO_3_2 2
#define ID_CO_4_3 3
#define ID_CO_5_4 4
#define ID_CO_6_5 5
#define ID_CO_7_6 6
#define ID_CO_8_7 7
#define ID_CO_9_8 8
#define ID_CO_10_9 9
#define ID_CO_11_10 10
#define ID_CO_12_11 11
#define ID_CO_13_12 12
#define ID_CII 13
#define ID_OIII88 14
#define ID_NII205 15
#define ID_NII122 16
#define ID_CI_1_0 17
#define ID_CI_2_1 18

#define ID_OIII5007 19

#define GAS 0
#define DM 1
#define TRACER 3
#define STAR 4
#define BH 5

// new parameter 
#define SFR_FIR 3.5e-44

// parameter for cloudy table//
#define N_METAL_CLOUDY             6
#define N_IONIZING_CLOUDY          93

// parameter for simulation code //
#define RUNIT                  1.0e-3
#define MUNIT                  1.0e10
#define SOFTNING_LENGTH        0.01

// parameter for subgroup data //
#define NMETAL                 10
#define NCONTINUUM             1221
#define NLINE                  124
#define NLINE_PEAGSE2          61
#define NCONTINUUM_DUST        200
#define NPHOTO 8

struct Subgroup_temp
{
	int id, parent_flag, subhalo_flag, group_id, npart[6];
	float pos[3], vel[3], pos_star[3];
	float mass, mass_part[6], vsigma, rmass_half, rmass_half_part[6], global_density;
	float mass_metal_gas, mass_dust, elementabundance_gas[NMETAL], sfr_gas;
	float mass_metal_star, elementabundance_star[NMETAL];
	float sfr[3], log_NLy, line_int[NLINE];
	float star_mass_weighted_age, mdot_bh, photometrics[NPHOTO];
};
//element abundance: H, He, C, N, O, Ne, Mg, Si, Fe, total
//photometrics: U, B, V, K, g, r, i, z

struct Continuum_line_data
{
  float continuum_int[NCONTINUUM], continuum_nebular_int[NCONTINUUM];
#ifdef LINE_DUST_DEPLETION
  float line_int[NLINE], line_int_new[NLINE];
#else
  float line_int[NLINE];
#endif
};

struct Hash_table
{
  int temp_id, ptype;
};

struct my_header
{
  int ngroup_total, nsubgroup_total;
  double redshift, boxsize;
  double omega, omegab, lambda, hubble;
};

struct Header
{
	int ngroup_total, nsubgroup_total;
	double redshift, boxsize, lumi_dis, ang_dis;
	double omega, omegab, lambda, hubble;
};

typedef struct
{
  double omega, lambda, omegab, hubble;
}Cosmology;

struct Dark_matter
{
  int id;
  float x[3], v[3];
};

struct Gas //PartType0
{
	int id;
	float x[3], v[3];
	float mass, density, sfr, temperature, internalenergy, hsml;
	float metallicity, elementabundance[NMETAL], e_abundance;
#ifdef FULL_SNAP
	float cooling_rate;
#endif
};

// element type : 0 = hydrogen, 1 = helium, 2 = carbon, 3 = nitrogen, 4 = oxygen, //
// 5 = neon, 6 = magnesium, 7 = silicon, 8 = iron //

struct Star
{
	int id;
	float x[3], v[3];
	float init_mass, mass, mass_component[3];
	float metallicity, elementabundance[NMETAL];
	float formation_time, age;
};

struct my_Star
{
	int id;
	float x[3], v[3];
	float init_mass, mass, metallicity, age;
	float log_Hb, log_NLy, log_habing_field;
	double line_int[NLINE], continuum_int[NCONTINUUM];
};

struct Black_Hole
{
	int id;
	float x[3], v[3];
	float mass, accretion_rate;
};

struct Unit
{
  double length, mass, velocity, density, energy, pressure, time;
};

struct Compare
{
  float aa, bb;
};

#endif
