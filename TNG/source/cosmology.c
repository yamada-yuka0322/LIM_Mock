#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <float.h>
#include "constant.h"
#include "para.h"
#include "proto.h"

double kick_factor(double time1, double time2, Cosmology cosm);
double drift_factor(double time1, double time2, Cosmology cosm);
double drift(double time, Cosmology cosm);
double dtda(double anow, Cosmology cosm);
void cosm_polint(double *xa, double *ya, int n, double x, double *y, double *dy);
double cosm_trapzd(double (*func)(double, Cosmology), Cosmology cosm, double a, double b, int n);
double kai_distance(Cosmology cosm,double z);
double d_kai(double z, Cosmology cosm);

void init_cosm(Cosmology *cosm, double om, double omb, double lmd, double h)
{
  cosm->omega = om;
  cosm->omegab = omb;
  cosm->lambda = lmd;
  cosm->hubble = h;
}

#define FUNC(x,y) ((*func)(x,y))
double cosm_trapzd(double (*func)(double, Cosmology), Cosmology cosm, double a, double b, int n)
{
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if(n == 1)
    {
      return (s = 0.5 * (b - a) * (FUNC(a, cosm) + FUNC(b, cosm)));
    }
  else
    {
      for(it = 1, j = 1; j < n - 1; j++)
        {
          it <<= 1;
        }
        
      tnm = it;
      del = (b - a) / tnm;
      x = a + 0.5 * del;
      
      for(sum = 0.0, j = 1; j <= it; j++, x += del)
        {
          sum += FUNC(x, cosm);
        }
      
      s = 0.5 * (s + (b - a) * sum / tnm);
      
    return s;
  }
}
#undef FUNC

void cosm_polint(double *xa, double *ya, int n, double x, double *y, double *dy)
{
  int i, m, ns = 1;
  double den, dif, dift, ho, hp, w;
  double *e, *d;
   
  double *etmp = (double *)malloc(n * sizeof(double));
  double *dtmp = (double *)malloc(n * sizeof(double));
  
  e = etmp - 1;
  d = dtmp - 1;
   
  dif = fabs(x - xa[1]);

  for(i = 1; i <= n; i++)
    {
      if((dift = fabs(x - xa[i])) < dif)
        {
          ns = i;
          dif = dift;
        }
    
      e[i] = ya[i];
      d[i] = ya[i];
    }

  *y = ya[(ns--)];
  
  for(m = 1; m < n; m++)
    {
      for(i=1; i <= n - m; i++)
        {
          ho = xa[i] - x;
          hp = xa[i + m] - x;
          w = e[i + 1] - d[i];
          den = ho - hp;
          
          if(den == 0.0)
            {
              printf("Error in routine cosm_polint\n");
              exit(1);
            }
        
          den = w / den;
          d[i] = hp * den;
          e[i] = ho * den;
        }
        
      *y += (*dy = (2 * ns < (n - m) ? e[ns + 1] : d[ns--]));
    }
   
  free((void *)etmp);
  free((void *)dtmp);
}

double dtda(double anow, Cosmology cosm)
{
  double eta;
  double om, ov;

  om = cosm.omega;
  ov = cosm.lambda;

  eta = sqrt(anow / (om + (1.0 - om - ov) * anow + ov * anow * anow * anow));
  
  return  eta;
}

double atotime(double anow, Cosmology cosm)
{
  double ss, dss;
  static double h[21], s[21];
  int i, k = 5;
  double eps = 1.0e-6;
  
  h[1] = 1.e0;
  
  for(i = 1; i <= 20; i++)
    {
      s[i] = cosm_trapzd(dtda, cosm, 0.0, anow, i);
      
      if(i >= k)
        {
          cosm_polint(&h[i-k], &s[i-k], k, 0.e0, &ss, &dss);
          
          if(fabs(dss) <= eps * fabs(ss))
            {
              return ss;
            }
        }
        
      s[i + 1] = s[i];
      h[i + 1] = 0.25 * h[i];
    }

  fprintf(stderr,"too many steps in atotime..\n");
  exit(1);
}

double ztotime(double znow, Cosmology cosm)
{
  double ss, dss, anow;
  static double h[21], s[21];
  int i, k = 5;
  double eps = 1.0e-6;
  
  anow = 1.0 / (1.0 + znow);

  h[1]=1.e0;
  
  for(i = 1; i <= 20; i++)
    {
      s[i] = cosm_trapzd(dtda, cosm, 0.0, anow, i);
      
      if(i >= k)
        {
          cosm_polint(&h[i-k], &s[i-k], k, 0.e0, &ss, &dss);
          
          if(fabs(dss) <= eps * fabs(ss))
            {
              return ss;
            }
        }
        
      s[i + 1] = s[i];
      h[i + 1] = 0.25 * h[i];
    }

  printf("too many steps in ztotime..\n");
  exit(1);
}

#define MAXIT 100
double timetoa(double tnow, Cosmology cosm)
{
  int j;
  double anow;
  double xacc;

  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;
  double x1, x2;

  xacc = 1.0e-7;

  x1 = 1.e-4;
  x2 = 2.0;

  fl = atotime(x1, cosm) - tnow;
  fh = atotime(x2, cosm) - tnow;

  if((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    {
      printf("Root must be bracketed in timetoa\n");
      exit(1);
    }
  
  if(fl == 0.0)
    {
      return x1;
    }
    
  if(fh == 0.0)
    {
      return x2;
    }

  if(fl < 0.0)
    {
      xl = x1;
      xh = x2;
    }
  else
    {
      xh = x1;
      xl = x2;
    }

  rts = 0.5 * (x1 + x2);
  dxold = fabs(x2 - x1);
  dx = dxold;
  f = atotime(rts, cosm) - tnow;
  anow = 1.0 / (1.0 + rts);
  df = -dtda(anow, cosm) / (1.0 + rts) / (1.0 + rts);

  for (j = 1; j <= MAXIT; j++)
    {
      if((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) || (fabs(2.0 * f) > fabs(dxold * df)))
        {
          dxold = dx;
          dx = 0.5 * (xh - xl);
          rts = xl + dx;
          
          if(xl == rts)
            {
              return rts;
            }
        }
      else
        {
          dxold = dx;
          dx = f / df;
          temp = rts;
          rts -= dx;
          
          if(temp == rts)
            {
              return rts;
            }
        }
        
      if(fabs(dx) < xacc)
        {
          return rts;
        }
        
      f = atotime(rts, cosm) - tnow;
      anow = 1.0 / (1.0 + rts);
      df = -dtda(anow, cosm) / (1.0 + rts) / (1.0 + rts);
      
      if(f < 0.0)
        {
          xl = rts;
        }
      else
        {
          xh = rts;
        }
    }
    
  printf("Maximum number of iterations exceeded in timetoa\n");
  exit(1);
}

double timetoz(double tnow, Cosmology cosm)
{
  int j;
  double anow;
  double xacc;

  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;
  double x1, x2;

  xacc = 1.0e-5;

  x1 = -0.1;
  x2 = 500.0;

  fl = ztotime(x1, cosm) - tnow;
  fh = ztotime(x2, cosm) - tnow;

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    {
      printf("Root must be bracketed in timetoz. tnow = %12.4e\n", tnow);
      printf("F(x1) = %14.6e at x1 = %14.6e\n", fl, x1);
      printf("F(x2) = %14.6e at x2 = %14.6e\n", fh, x2);
      exit(1);
    }
  
  if(fl == 0.0)
    {
      return x1;
    }
    
  if(fh == 0.0)
    {
      return x2;
    }

  if(fl < 0.0)
    {
      xl = x1;
      xh = x2;
    }
  else
    {
      xh = x1;
      xl = x2;
    }
    
  rts = 0.5 * (x1 + x2);
  dxold = fabs(x2 - x1);
  dx = dxold;
  f = ztotime(rts, cosm) - tnow;
  anow = 1.0 / (1.0 + rts);
  df = -dtda(anow, cosm) / (1.0 + rts) / (1.0 + rts);

  for(j = 1; j <= MAXIT; j++)
    {
      if((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) || (fabs(2.0 * f) > fabs(dxold * df)))
        {
          dxold = dx;
          dx = 0.5 * (xh - xl);
          rts = xl + dx;
          
          if(xl == rts)
            {
              return rts;
            }
        }
      else
        {
          dxold = dx;
          dx = f / df;
          temp = rts;
          rts -= dx;

          if(temp == rts)
            {
              return rts;
            }
        }

      if(fabs(dx) < xacc)
        {
          return rts;
        }
        
      f = ztotime(rts, cosm) - tnow;
      anow = 1.0 / (1.0 + rts);
      df = -dtda(anow, cosm) / (1.0 + rts) / (1.0 + rts);
      
      if(f < 0.0)
        {
          xl = rts;
        }
      else
        {
          xh = rts;
        }
    }

  printf("Maximum number of iterations exceeded in timetoz\n");
  exit(1);
}
#undef MAXIT

double drift(double time, Cosmology cosm)
{
  double anow;
  anow = timetoa(time, cosm);

  return (1.0 / anow / anow);
}

double kick(double time, Cosmology cosm)
{
  double anow;
  anow = timetoa(time, cosm);

  return (1.0 / anow);
}


double drift_factor(double time1, double time2, Cosmology cosm)
{

  double ss, dss;
  static double h[21], s[21];
  int i, k = 5;
  double eps = 1.0e-6;

  h[1] = 1.e0;

  for(i = 1; i <= 20; i++)
    {
      s[i] = cosm_trapzd(drift, cosm, time1, time2, i);
      
      if(i >= k)
        {
          cosm_polint(&h[i-k], &s[i-k], k, 0.e0, &ss, &dss);
          
          if(fabs(dss) <= eps * fabs(ss))
            {
              return ss;
            }
        }
        
      s[i + 1] = s[i];
      h[i + 1] = 0.25 * h[i];
    }

  printf("too many steps in drift_factor..\n");
  exit(1);
}

double kick_factor(double time1, double time2, Cosmology cosm)
{

  double ss, dss;
  static double h[21], s[21];
  int i, k = 5;
  double eps = 1.0e-6;

  h[1] = 1.e0;
  
  for(i = 1; i <= 20; i++)
    {
      s[i] = cosm_trapzd(kick, cosm, time1, time2, i);
      
      if(i >= k)
        {
          cosm_polint(&h[i - k], &s[i - k], k, 0.e0, &ss, &dss);
          
          if(fabs(dss) <= eps * fabs(ss))
            {
              return ss;
            }
        }

      s[i + 1] = s[i];
      h[i + 1] = 0.25 * h[i];
    }

  printf("too many steps in kick_factor..\n");
  exit(1);
}

double d_kai(double z, Cosmology cosm)
{
   double delta_kai;
   double om, ov, zp1;

   om  = cosm.omega;
   ov  = cosm.lambda;
   zp1 = 1.0 + z;

   delta_kai = 1.0 / sqrt(om * zp1 * zp1 * zp1 + (1.0 - om - ov) * zp1 * zp1 + ov);

   return delta_kai;
}

double kai_distance(Cosmology cosm, double z)
{
  double kai;
  double ss, dss, eps;
  static double h[21], s[21];
  int j , k = 5;
  
  if(z < 1.0e-5)
    {
      kai = z;
      return kai;
    }
  
  eps = 1.0e-6;
   
  h[1] = 1.0;

  for(j = 1; j <= 20; j++)
    {
      s[j] = cosm_trapzd(d_kai, cosm, 0.0, z, j);
      
      if(j >= k)
        {
          cosm_polint(&h[j-k], &s[j-k], k, 0.0, &ss, &dss);
          
          if(fabs(dss) < eps * fabs(ss))
            {
	      kai = ss;
	      return kai;
            }
        }
        
    s[j + 1] = s[j];
    h[j + 1] = 0.25 * h[j];
  }

  fprintf(stderr,"too many steps in kai_distance \n");
  exit(1);
}

double comoving_distance(Cosmology cosm, double z)
{
  double ss, k0, sqk0;
  double r_z;
  
  r_z = 0.0;

  ss = kai_distance(cosm, z);

  k0 = cosm.omega + cosm.lambda - 1.0;
  sqk0 = sqrt(fabs(k0));

  if(k0 < 0.0 && fabs(k0) > 0.01)
    {
      r_z = sinh(sqk0 * ss) / sqk0 * cspeed / (H0 * cosm.hubble);
    }
  else if(fabs(k0) < 0.01)
    {
      r_z = ss * (cspeed) / (H0 * cosm.hubble);
    }
  else if(k0 > 0.0 && fabs(k0) > 0.01)
    {
      r_z = sin(sqk0 * ss) / sqk0 * cspeed / (H0 * cosm.hubble);
    }

  return r_z;
}

double angular_distance(Cosmology cosm, double z)
{
  double d_A;
  double k0, sqk0;
  double ss;

  d_A = 0.0;

  k0 = cosm.omega + cosm.lambda - 1.0;
  sqk0 = sqrt(fabs(k0));
  ss = kai_distance(cosm, z);

  if(k0 < 0.0 && fabs(k0) > 0.01)
    {
      d_A = sinh(sqk0 * ss) / sqk0 * cspeed / (H0 * cosm.hubble) / (1.0 + z);
    }
  else if(fabs(k0) < 0.01)
    {
      d_A = ss * (cspeed) / (H0 * cosm.hubble) / (1.0 + z);
    }
  else if(k0 > 0.0 && fabs(k0) > 0.01)
    {
      d_A = sin(sqk0 * ss) / sqk0 * cspeed / (H0 * cosm.hubble) / (1.0 + z);
    }

  return d_A;
}

double luminosity_distance(Cosmology cosm, double z)
{
  double d_L;

  d_L = angular_distance(cosm, z) * (1.0 + z) * (1.0 + z);

  return d_L;
}
/*int main(void)
{
  Cosmology cosm;
  double H0, H0inv, tinit, tfin, tnow, dt, anow, zinit;
  int i;

  cosm.omega = 0.3;
  cosm.lambda = 0.7;
  cosm.hubble = 0.7;
  H0 = 3.24e-18 * cosm.hubble;
  H0inv = 1.0 / H0 / yr;
  
  zinit = 49.0;
  dt = 1.0e7 / H0inv;
  tinit = ztotime(zinit, cosm);
  tfin = ztotime(0.0, cosm);
  tnow = tinit;
  
  //printf("%e %e %e %e\n", tinit, tfin, dt, H0inv);
  
  i = 0;
  
  while(tnow <= tfin)
    {
      tnow = tinit + dt * (double) i;
      anow = timetoa(tnow, cosm);
      
      printf("%e\n", anow);
      
      i++;
    }

  return 0;
}*/
