# Codes for creating mocks for SUBLIME

Create mocks with the following steps. See job.sh to see what to do in more
detail.

If you have lightcone data, skip steps 0-2 and go to step 3.

## 0. Preparation

### 0-0. Create fof catalogs

Use job_fof.sh to run output_fof_data.c in the upper directory. This creates 
header.txt and subgroup_temp_data.txt (note these are not ascii but binary)

### 0-1. Create comving distance table
If a file "cdis_table.txt" does not exist, create it using make_comdis_table.c.
This should be a text file with comovig distance and redshift in its first and second column. 
This text file is required when creating lightcone. 
This table is used to get redshift from comoving distance, so it is better to have fine resolution for redshift. 

### 0-2. Create snapshot id list
If a file "snapshot_ids.txt" does not exsit, create it. 
This should be a text file with the snapshot number and the fraction of volume to be used in its first and second column.
This text file is required when creating lightcone.

### 0-3. Define lines
Define the line names and the function in line_model.py

### 0-3'. Define Lines (using c files; old version)
Define the followings
- line ids in para.h
- rest-frame wavelengths and L-SFR relations in create_lightcone.c

### 0-3. Makefile

Makefile option:
- USE_FILTER_TABLE: use pre-defined list of observed frequencies

### 0-3.

## 1. Print subgroup properties

Run print_subgroup.c. This reads header.txt and subgroup_temp_data.txt and output a galaxy catalog

Output is an ascii file containing the following: 
0: x [Mpc/h]
1: y [Mpc/h]
2: z [Mpc/h]
3: vx [km/s]
4: vy [km/s]
5: vz [km/s]
6: log sfr 
7: log mstar
8: log Z

## 2. Create lightcone data 

Run  
```
./create_lightcone output_path side_length zmin zmax x0 y0 tan_LoS_x tan_LoS_y
```
where side_length, x0, and y0 are in units of [arcsec], [Mpc/h], and [Mpc/h] to generate a lightcone with side length of side_length.

Change the function `line_luminosity_CO` to use 
different line luminosity calculation method.  


Output:
0: theta_x [arcsec] 
1: theta_y [arcsec] 
2: redshift -- not corrected with vz
3: log luminosity distance [cm]
4: vz [km /s ] -- LoS velocity, galaxies going away from us have vz > 0 
5: log sfr 
6: log mstar
7: log Z

### 3. Create mock line intensity map with python files

Run 
```
python create_mock.py 
```

Check parameters that can be changed by running
```
python create_mock.py --h
```

The flexible filter method is not yet implemented. 

The intensity maps generated can be plotted with plot.ipynb

Set z_dependenc_alpha = 1 to use the original SFR in the input file. If it is not unity, larger/smaller SFR is used for high redshift.


## 3'. Create mock line intensity map (using c files; old version)

Run 
```
./create_mock input_fname side_length fmin fmax angular_resolution R z_dependence_alpha fout_id threshold
```
where side_length, fmin, fmax, and angular_resolution are in unit of [arcsec], [Hz], [Hz], and [arcsec].

If you want to use a filter table, add USE_FILTER_TABLE in
the Makefile option and specify the name of table. 

--- (Example of the filter table ) ---
300 #<- minimum frequeny
310
320
330
340
350
360
370
380
390 #< maximum frequency
--- (End of example) ---
Frequencies are edge values (i.e., you need Nfilter+1 values)

Set z_dependenc_alpha = 1 to use the original SFR in the input file. If it is not unity, larger/smaller SFR is used for high redshift.

Set threshold < 0 to create intensity map. If threshold >= 0, then the number 
density of galaxies with flux > threshold [erg/s/cm2] will be created.

The intensity maps generated can be plotted with xxx/plot.ipynb



## LEGACY:

calc_ngal:
