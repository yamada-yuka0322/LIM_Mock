#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <float.h>
#include <mpi.h>
#include "constant.h"
#include "para.h"
#include "proto.h"

void output_header_data(char *path, int ngroup_total, int nsubgroup_total, struct io_header header, Cosmology cosm)
{
  char fout[256];
  struct my_header temp;
  FILE *fp;
  
  sprintf(fout, "%s/header.txt", path);
  fp = fopen(fout, "wb");
  
  if(fp == NULL)
    {
      printf("can not open outputfile : %s\n", fout);
      exit(1);
    }
    
  temp.ngroup_total = ngroup_total;
  temp.nsubgroup_total = nsubgroup_total;
    
  temp.redshift = header.redshift;
  temp.boxsize = header.BoxSize;
  temp.omega = header.Omega0;
  temp.omegab = cosm.omegab;
  temp.lambda = header.OmegaLambda;
  temp.hubble = header.HubbleParam;
    
  fwrite(&temp, sizeof(struct my_header), 1, fp);
  
  fclose(fp);
}

void read_header_data(char *path, struct my_header *header)
{
  char fin[256];
  FILE *fp;
  
  sprintf(fin, "%s/header.txt", path);
  fp = fopen(fin, "rb");
  
  if(fp == NULL)
    {
      printf("can not open inputfile : %s\n", fin);
      exit(1);
    }
    
  fread(header, sizeof(struct my_header), 1, fp);
  
  fclose(fp);
}
