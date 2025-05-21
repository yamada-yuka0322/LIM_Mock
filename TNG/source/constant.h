
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

// Astronomical parameter //
#define au               1.496e13        // astronomical unit (cm) //
#define Msun             1.989e33        // solor mass //
#define Lsun             3.826e33        // solor luminosty (erg/s) //
#define Zsun             0.0127            // solor metallicity (for TNG)//
#define cspeed           2.998e10        // light speed(cm/s) //
#define H0               3.226e-18       // Hubble constant(/s) //
#define yr               3.156e7         // 1 year(s) //
#define Gyr              3.156e16         // 1 year(s) //
#define pc               3.086e18        // parsec(cm) //
#define kpc              3.086e21        // kiroparsec(cm) //
#define mpc              3.086e24        // megaparsec(cm) //
#define km               1.0e5           // kirometer(cm) //
#define Angstrom         1.0e-8          // angstrom(cm) //
#define ANGSTROM2MICRON  0.0001
#define MICRON2ANGSTROM  1.0e4
#define micron           1.0e-4          // micron(cm) //
#define angstrom         1.0e-8          // angstrom(cm) //
#define Jy               1.0e-23         // jansky (erg/s/cm2/Hz)//
#define mJy              1.0e-26         
#define GHz              1.0e9           // Hz //
#define G0               1.6e-3          // habing field (erg/cm2/s)

// physics parameters //
#define hplanck          6.626e-27       // Planck constant (erg s = g cm2 / s)//
#define G                6.673e-8        // gravitational constant (cm3 g-1 s-2)//
#define k_B              1.381e-16       // Boltzmann constant (erg/K = g cm2 / s2 / K)//
#define m_p              1.6726219e-24   //proton[g]//
#define alpha_B          2.6e-13         //recombination coefficient [cm3/s]
#define eV               1.602176e-12    // electron volt (erg)

// mathematics parameters //
#define PI               3.141592653589    // 2.0 * acos(0.0) //
#define arcmin           2.908882086656e-4 // [rad] ... PI / 180 / 60 //
#define arcsec           4.848136811094e-6 // [rad] ... arcmin / 60 //
