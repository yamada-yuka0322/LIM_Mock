SYSTYPE ?= "MyLapTop"
ifeq ($(SYSTYPE), Docker)
CC       =  cc
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib -lgsl -lgslcblas
FITS_INCL = -I/usr/include
FITS_LIBS = -L/usr/lib -lcfitsio
HDF5INCL = -I/usr/include/hdf5/serial
HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lz
endif

ifeq ($(SYSTYPE), "MyLapTop")
CC       =  cc
OPTIMIZE = -O3 -Wall
#GSL_INCL = -I/${HOME}/local/include
#GSL_LIBS = -L/${HOME}/local/lib -lgsl -lgslcblas
GSL_INCL = -I/opt/homebrew/include
GSL_LIBS = -L/opt/homebrew/lib -lgsl -lgslcblas
FITS_INCL = -I/opt/homebrew/include
FITS_LIBS = -L/opt/homebrew/lib -lcfitsio
#HDF5INCL = -I/opt/local/include
#HDF5LIB  = -L/opt/local/lib -lhdf5 -lz
HDF5INCL = -I/opt/homebrew/include
HDF5LIB  = -L/opt/homebrew/lib -lhdf5 -lz
endif

ifeq ($(SYSTYPE), "MYiMac")
CC       =  cc
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/opt/homebrew/opt/gsl/include
GSL_LIBS = -L/opt/homebrew/opt/gsl/lib -lgsl -lgslcblas -lcfitsio
#HDF5INCL = -I/${HOME}/local/include
#HDF5LIB  = -L/${HOME}/local/lib -lhdf5 -lz
HDF5INCL = -I/opt/homebrew/opt/hdf5/include
HDF5LIB  = -L/opt/homebrew/opt/hdf5/lib -lhdf5 -lz
endif

CFLAGS = $(HDF5INCL) $(FFTW_INCL) $(FITS_INCL) $(GSL_INCL) $(OPTIMIZE) $(DFLAGS) 
LIBS   = -lm $(HDF5LIB) $(FFTW_LIB) $(FITS_LIBS) $(GSL_LIBS)

# etc
SHELL	= /bin/sh
RM	= rm -f
TAR	= tar -cvzf
ZIP	= zip
 
SRCH = constant.h para.h proto.h

PRINT_SUBGROUPOBJ = cosmology.o print_subgroup.o
CREATE_LIGHTCONEOBJ = cosmology.o line_model.o create_lightcone.o
CREATE_LIGHTCONE_UCHUUOBJ = cosmology.o line_model.o create_lightcone_UCHUU.o
CREATE_MOCKOBJ = cosmology.o line_model.o create_mock.o 
CALC_INTENSITY_DISTRIBUTIONOBJ = cosmology.o line_model.o calc_intensity_distribution.o

ALL =  print_subgroup create_lightcone create_lightcone_UCHUU create_mock calc_intensity_distribution 

all: $(ALL)


print_subgroup: $(PRINT_SUBGROUPOBJ)
	$(CC) $(CFLAGS) $(PRINT_SUBGROUPOBJ) $(LIBS) -o print_subgroup

create_lightcone: $(CREATE_LIGHTCONEOBJ)
	$(CC) $(CFLAGS) $(CREATE_LIGHTCONEOBJ) $(LIBS) -o create_lightcone

create_lightcone_UCHUU: $(CREATE_LIGHTCONE_UCHUUOBJ)
	$(CC) $(CFLAGS) $(CREATE_LIGHTCONE_UCHUUOBJ) $(LIBS) -o create_lightcone_UCHUU
	
create_mock: $(CREATE_MOCKOBJ)
	$(CC) $(CFLAGS) $(CREATE_MOCKOBJ) $(LIBS) -o create_mock

calc_intensity_distribution: $(CALC_INTENSITY_DISTRIBUTIONOBJ)
	$(CC) $(CFLAGS) $(CALC_INTENSITY_DISTRIBUTIONOBJ) $(LIBS) -o calc_intensity_distribution

$(PRINT_SUBGROUPOBJ): $(SRCH)
$(CREATE_LIGHTCONEOBJ): $(SRCH)
$(CREATE_MOCKOBJ): $(SRCH)
$(CALC_INTENSITY_DISTRIBUTIONOBJ): $(SRCH)

clean:
	$(RM) $(ALL) 
	$(RM) *.o 	

.c.o:
	$(CC) $(CFLAGS) -c $<
