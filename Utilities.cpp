#define _USE_MATH_DEFINES

#include <cmath>
#include <stdio.h>
#include <string.h>
//#include <math.h>       /* atan */

/*
#include "include\R.h"
#include "include\Rmath.h"
#include "include\Rinternals.h"

#include <R_ext/Applic.h> 
#include <R_ext/Utils.h>
#include <R_ext/Error.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h> 
*/

#define M_PI 3.14159265358979323846



/****************************************************************
 * max/min -- find the maximum/minimum element in an unsorted vector
 *
 * Parameters: vector -- pointer to the vector in which to look for
 *      the maximum (double)
 *
 *      max -- pointer to the variable to store the maximum (double)
 *
 * Returns
 *    integer index of element that is the maximum
 ***************************************************************/

int max(double *vector,
	int num_elements){

  int location  = 0; // index of maximum element 
  double maximum; // the value of the maximum

  maximum = *(vector + location);

  for (int i = 1; i < num_elements; ++i){
  
    if(*(vector + i) > maximum){
      maximum = *(vector + i);
      location = i;
    }
  }
 
  return(location);
}


int min(double *vector,
	int num_elements){

  int location  = 0; // index of minimum element 
  double minimum; // the value of the minimum

  minimum = *(vector + location);

  for (int i = 1; i < num_elements; ++i){
  
    if(*(vector + i) < minimum){
      minimum = *(vector + i);
      location = i;
    }
  }
 
  return(location);
}



// min and max functions for two numbers a and b
static inline double max2(double a, double b){
  
  return((a > b) ? a : b);

}

static inline double min2(double a, double b){
  
  return((a < b) ? a : b);

}


// Calculate day length using Sinclair calculations
static inline double dayLength(double latitude, int dofy){

  // maybe replace 3.1416 with M_PI later for more accuracy

  int doy = dofy;
  double lat = latitude * (3.1416 / 180);
  double sunset = -6 * 3.1416 / 180;
  double decl = (-23.4 * cos(2.0 * 3.1416 * (doy + 10) / 365)) * 3.1416 / 180;
	 
  double tmp1 = (sin(sunset) - sin(lat) * sin(decl)) / (cos(lat) * cos(decl));
  double tmp2 = tmp1 / pow(1 - pow(tmp1,2), 0.5);

  double dleng = 24 / 3.1416 * (3.1416 / 2 - atan(tmp2));
    
  return(dleng);
}


// relative solar radiation
static inline void relSolRad(double dlength, 
			   double *rel_srad_out){
    
  int sunrise = roundf(12 - 0.5 * dlength);
  int sunset = roundf(12 + 0.5 * dlength);
  
  double rel_srad = 0;

  for(int h = 0; h < (sunset - sunrise + 1); ++h){

    rel_srad = max2(0, M_PI / (2.0 * dlength) * sin(M_PI * (h / dlength )));
        
    rel_srad_out[sunrise + h - 1] = rel_srad;

  }    
}



// Calculate VPD from tmax and tmin. Output VPD and DELT
static inline void hourlyVPD(double tmax, 
			     double tmin,
			     double *vpd_out,
			     double *tem_out){
    
  tmin >= tmax && (tmax = tmin);

  double tsamp = (tmax - tmin) / 2.0;

  double vpmn = 6.107 * exp(17.269 * tmin / (237.3 + tmin)); // in mbar or hPa 

  double tem_h, vp_h, vpd_h;

	
  for(int h = 1; h <= 24; ++h){
    tem_h = (tmax + tmin) / 2.0 - tsamp * sin( 2.0 * M_PI * (h + 3) / 24);
    vp_h = 6.107 * exp(17.269 * tem_h / (237.3 + tem_h));
    vpd_h = vp_h - vpmn;

    vpd_out[h-1] = vpd_h;
    tem_out[h-1] = tem_h;
  }
}

static inline void hourlyRueMax(double rueMax,
				double *tem_h,
				double *rueMax_h){

  for(int h = 1; h <= 24; ++h)
    rueMax_h[h-1] = max2(min2(0.069 * (tem_h[h - 1] - 6.6), 1), 0) * rueMax;

}


// calculates leaf number
static inline double leafNum(double dev_state_termleaf, 
			     double co_leaf1,
			     double co_leaf2){

  return (co_leaf1 * exp(dev_state_termleaf * co_leaf2));

}

// calculates leaf area
static inline double leafArea(double leaf_number,
			      double amax,
			      double tlno){
  
  double lnm = 3.53 + 0.46 * tlno;
  double l_num;
  double pla = 0;
        
  l_num = (leaf_number <= (tlno - 2)) ? ( (leaf_number) + 2) : tlno;
  
  
  for(int i = 1; i <= ceil(l_num); ++i){
    
    if(i >= floor(l_num))
      pla += amax * exp( -0.0344 * pow(i - lnm, 2) + 0.000731 * pow(i - lnm, 3)) *
	(leaf_number - floor(leaf_number));
    else
      pla += amax * exp( -0.0344 * pow(i - lnm, 2) + 0.000731 * pow(i - lnm, 3));

  }
  return(pla);
}


static inline void carbonWater(double coefext, 
			       double *rue_max_h, 
			       double lai,
			       double srad,
			       double *r_srad, 
			       double *vpd_h, 
			       double dev_state_seed,
			       double tecoeff, 
			       double vpd_bkp, 
			       double vpd_slope,
			       double w_supply,
			       double *phs_transp_demand){

  // set values in output array to zero
  phs_transp_demand[0] = 0.0;
  phs_transp_demand[1] = 0.0;
  phs_transp_demand[2] = 0.0;

  double frac_int_rad = 1 - exp(-coefext * lai);

  double rue_act = 0;
        
  vpd_bkp = vpd_bkp * 10;  	// Convert kPa to hPa

  double phs_h, transp_h, vpd_p, w_demand_h, s_d_ratio;
    
  double w_supply_h = w_supply;

  // limited transpiration here
  for(int h = 1; h <= 24; ++h){
            
    rue_act = (dev_state_seed <= 722) ? rue_max_h[h-1] : (0.75 * rue_max_h[h-1]);

    w_supply_h = w_supply_h - transp_h;

    if(w_supply_h == 0)
      w_supply_h += 0.000001;

    if(vpd_h[h - 1] <= vpd_bkp){
                
      w_demand_h = rue_act * (r_srad[h - 1] * srad) * frac_int_rad * vpd_h[h - 1] / (tecoeff * 1000);

      s_d_ratio = max2(0.0, min2(1.0, w_supply_h/w_demand_h));

      phs_h = rue_act * (r_srad[h - 1] * srad) * frac_int_rad * s_d_ratio;


      // 0.09 hPa because vpd is hPa. teCoeff maize = 9
      // Pa. factor 1000 is to convert g mass to mm 10^4
      // cm2/m2 10mm/cm 1 cm3/g
      transp_h = phs_h * vpd_h[h - 1] / (tecoeff * 1000);
    } 
    else{
      vpd_p = vpd_bkp + vpd_slope * (vpd_h[h - 1] - vpd_bkp);

      w_demand_h = rue_act * (r_srad[h - 1] * srad) * frac_int_rad * vpd_p / (tecoeff * 1000);
      w_demand_h += 0.000001;

      s_d_ratio = max2(0.0, min2(1.0, w_supply_h/w_demand_h));      

      phs_h = rue_act * (r_srad[h - 1] * srad) * frac_int_rad * s_d_ratio * vpd_p / vpd_h[h - 1];

      transp_h = phs_h * vpd_p / (tecoeff * 1000);
    }
    //Rprintf("frac_int_rad %f srad %f r_srad %f sdration %f phs_h %f\n", frac_int_rad, srad, r_srad[h-1],s_d_ratio,phs_h);
    

    phs_transp_demand[0] += phs_h;
    phs_transp_demand[1] += transp_h;
    phs_transp_demand[2] += w_demand_h;
            
  }
}


//-------------------------------
// development rates

// planting emergence
static inline double devRateEme(double tmean, 
				double tbe){
  return(max2(0.0, tmean - tbe));
}


// planting to vegetative state (?)
static inline double devRateVeg(double tmean,
				double tbv){
	  return(max2(0.0, tmean - tbv));
}


// planting to seed set (?)
static inline double devRateSeed(double tmean, 
				 double tbr){
  return(max2(0.0, tmean - tbr));
}


// planting to root growth (?)
static inline double devRateRoot(double rootmax, 
				 double tmax, 
				 double tmin, 
				 double tbroot){
  return(rootmax); // cm/day
}


// daily GDU  
static inline double dGDU(double tmax, double tmin){

  // convert to Fahrenheit
  tmax = (tmax * 1.8) + 32;
  tmin = (tmin * 1.8) + 32;

  tmax = min2(tmax, 95.0);
  tmin = max2(tmin, 50.0);

  return((tmax + tmin) / 2 - 50);

}



// evaporation function
static inline double evaporation(double dyse, 
				 double tmax,
				 double tmin, 
				 double srad,
				 double prcp,
				 double w_vpd,
				 double lai,
				 double esw,
				 double mulch_effect){
		
  double t_tmean = (tmax + tmin) / 2;
  double delt = exp(21.255 - 5304 / (273 + t_tmean)) * (5304 / pow(273 + t_tmean, 2));
                
  // min evaporation 1.5 mm Muchow et al. 83:1052-1059 (1991) 
  double soil_evp = max2(1.5,
			 (exp(-0.4 * lai) * srad * delt + 0.68 * 0.4 * w_vpd) / 
			 (delt + 0.68 ) * 239 / 538 * mulch_effect); 	
		
  double esw_tmp = esw - soil_evp;


  // Muchow et al. 83:1052-1059 (1991) 		
  if(!(dyse == 1 && esw_tmp >= 0))
    soil_evp *= (sqrt(dyse + 1) - sqrt(dyse)); 
      
  (esw_tmp <= 0) && (soil_evp = esw);

  // only output is soil evaporation to update soil
  // water balance
  return(soil_evp);
}

// silk elongation rate
static inline double ser(double *co_ser, 
			 double ftsw){
  
  return (co_ser[0] + (co_ser[1] - co_ser[0]) / (1 + exp((co_ser[2] - ftsw) / co_ser[3]))) / 
    (co_ser[0] + (co_ser[1] - co_ser[0])/(1 + exp((co_ser[2] - 1)/co_ser[3])));
  
}




// EOF
  
  
