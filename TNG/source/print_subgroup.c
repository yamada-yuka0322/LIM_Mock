#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

#include "constant.h"
#include "para.h"
#include "proto.h"

#include "mt19937ar.c"
#include "mt19937ar.h"

void read_header_data(char *path, struct Header *header);
void read_subgroup_temp_data(char *path, struct Subgroup_temp **subgroup_temp, struct Header header);


int main(int argc, char *argv[])
{
	if(!(argc == 2))
	{
		printf("Please input and input_dir\n");
		exit(1);
	}
	char path[256];
	sprintf(path, "%s", argv[1]);

	int seed = 12345;
	srand(seed);
	init_genrand(seed);

	struct Header header;
	struct Subgroup_temp *subgroup_temp;
	fprintf(stderr, "#\n");
	read_header_data(path, &header);
	read_subgroup_temp_data(path, &subgroup_temp, header);
		
	double total_sfr = 0.0;

	printf("#z: %.3f\n", header.redshift);
	printf("#boxsize: %.3f\n", header.boxsize/1e3);
	printf("#ngal: %d\n", header.nsubgroup_total);
	printf("#x y z vx vy vz rmass_half sfr mstar metallicity\n");
	for(int i = 0; i < header.nsubgroup_total; i++)
	{
		if( subgroup_temp[i].mass_part[DM] * MUNIT > 1e9 )
		{
			if( subgroup_temp[i].sfr_gas <= 0 ) continue;

			printf("%f %f %f %f %f %f %f %f %f %f\n"
				, subgroup_temp[i].pos[0] * RUNIT
				, subgroup_temp[i].pos[1] * RUNIT
				, subgroup_temp[i].pos[2] * RUNIT
				, subgroup_temp[i].vel[0] 
				, subgroup_temp[i].vel[1] 
				, subgroup_temp[i].vel[2] 
				, subgroup_temp[i].rmass_half_part[STAR] * RUNIT
				, log10( subgroup_temp[i].sfr_gas )
				, log10( subgroup_temp[i].mass_part[STAR] * MUNIT / header.hubble ) // Msun
				, log10( subgroup_temp[i].mass_metal_gas / subgroup_temp[i].mass_part[GAS] )
			  );
		}
		
		total_sfr += subgroup_temp[i].sfr_gas;
	}
	
	FILE *fp = fopen("./sfrd.txt", "a");
	double volume = pow(header.boxsize / 1e3 / header.hubble, 3); // [Mpc3]
	fprintf(fp, "%f %e\n", header.redshift, total_sfr / volume);
	fclose(fp);

	free(subgroup_temp);
	
	return 0;
}


void read_header_data(char *path, struct Header *header)
{
  char fin[256];
  FILE *fp;
  struct my_header header_s;
  Cosmology lcdm;
  
  sprintf(fin, "%s/header.txt", path);
  fprintf(stderr, "# read %s\n", fin);
  fp = fopen(fin, "rb");
  
  if(fp == NULL)
    {
      printf("can not open inputfile : %s\n", fin);
      exit(1);
    }
    
  fread(&header_s, sizeof(struct my_header), 1, fp);
 
  fclose(fp);

  header->redshift = header_s.redshift;
  header->boxsize = header_s.boxsize;
  header->ngroup_total = header_s.ngroup_total;
  header->nsubgroup_total = header_s.nsubgroup_total;
  header->omega = header_s.omega;
  header->omegab = header_s.omegab;
  header->lambda = header_s.lambda;
  header->hubble = header_s.hubble;

  init_cosm(&lcdm, header->omega, header->omegab, header->lambda, header->hubble);
  header->lumi_dis = luminosity_distance(lcdm, header->redshift);
  header->ang_dis = angular_distance(lcdm, header->redshift);

}

void read_subgroup_temp_data(char *path, struct Subgroup_temp **subgroup_temp, struct Header header)
{
	int i;
	FILE *fp;
	char input[256];
	struct Subgroup_temp temp;

	*subgroup_temp = malloc(sizeof(struct Subgroup_temp) * header.nsubgroup_total);
	if(*subgroup_temp==NULL)
	{
		printf("can not allocate subgroup\n");
		exit(1);
	}
	
	sprintf(input, "%s/subgroup_temp_data.txt", path);
	fprintf(stderr, "# read %s\n", input);

	fp = fopen(input, "rb");
	if(fp==NULL)
	{
		printf("can not open file : %s\n", input);
		exit(1);
	}

	for(i=0;i<header.nsubgroup_total;i++)
	{
		fread(&temp, sizeof(struct Subgroup_temp), 1, fp);
		(*subgroup_temp)[i] = temp;
	}

	fclose(fp);
}
