
sim_name=TNG300-1
odir=../output_${sim_name}

mkdir -p $odir

### output subgroup data
for i in 15 # 86 80 76 72 67 63 59 55 52 49 46 43 40 38 35 33 31 29 27 25 24 23 21 19 17 15 14 13 12 11 10 9 8 7 6 5 4 3 2
do
	idir=../subgroup_data_${sim_name}/${i}
	#./print_subgroup $idir > ${odir}/subgroup.$i.txt
	#echo "output ${odir}/subgroup.$i.txt"
done

### create lightcone
side_length_lc=10800 #[arcsec]
zmin=0.17 # minimum redshift of lightcone
zmax=3.0 # maximum redshift of lightcone
x0=0.0 # x position of the observer [cMpc/h]: x0 and y0 should be changed when generating multiple realizations.
y0=0.0 # y position of the observer [cMpc/h]
tan_LoS_x=0.06 # LoS direction: tan_LoS is needed to avoid having similar structures in a line of site
tan_LoS_y=0.1 # LoS direction
label= "" # label of the lightcone. Define only when needed
./create_lightcone_UCHUU ${odir}/ $side_length_lc $zmin $zmax $x0 $y0 $tan_LoS_x $tan_LoS_y > ${odir}/lightcone_${side_length_lc}sec${label}.txt 2> ${odir}/log_${side_length_lc}sec${label}.txt


### create mock data from lightcone
fname_lightcone=$odir/lightcone_${side_length_lc}sec${label}.txt # lightcone name to use

model_id=0 # 0: default -- DeLooze+14, sigma=0.2
#model_id=1 # 1: Vallini+15, sigma=0.2
#model_id=10 # 10: DeLooze+14, sigma=0.4
#model_id=11 # 11: Vallini+15, sigma=0.4

side_length=10800 #[arcsec]
fmin=6e4 #[GHz]
fmax=4e5 #[GHz]
angular_resolution=4.2 #[arcsec]
spectral_resolution=100 # R = f/df = lambda / dlambda
z_dependence_alpha=15 # factor for sfr at high redshift

#fmin=78000
#fmax=95100

threshold=-1 # intensity map if threshold < 0
fout_id=${odir}/map_${side_length}sec_fmin${fmin}_fmax${fmax}_angular${angular_resolution}_spectral${spectral_resolution}${label}_model${model_id}_alpha${z_dependence_alpha}
#./create_mock $fname_lightcone $model_id $side_length $fmin $fmax $angular_resolution $spectral_resolution $z_dependence_alpha $fout_id $threshold

threshold=1e-18 # [erg/s/cm2] galaxy number density map if threshold >= 0
fout_id=${odir}/nden_map_${side_length}sec_fmin${fmin}_fmax${fmax}_angular${angular_resolution}_spectral${spectral_resolution}${label}_model${model_id}_alpha${z_dependence_alpha}
#./create_mock $fname_lightcone $model_id $side_length $fmin $fmax $angular_resolution $spectral_resolution $z_dependence_alpha $fout_id $threshold

python3 create_mock.py --path $fname_lightcone --output_dir $odir --side_length $side_length --fmin $fmin --fmax $fmax --angular_resolution $angular_resolution --R $spectral_resolution --z_dependence_alpha $z_dependence_alpha --threshold -1
