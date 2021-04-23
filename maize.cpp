/* maizeGrainYieldModel

Frank Technow, 9/5/2017, version 2.1

change log: 

version 2.1
- implemented silking module

*/
//#define _USE_MATH_DEFINES

//#include <cmath>

#include <stdio.h>
#include <string.h>
//#include <math.h>       /* atan */
/*
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h> 
#include <R_ext/Utils.h>
#include <R_ext/Error.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h> 

#include "maizeGrainYieldModel_functions_2.2.h"
*/
#include "Utilities.cpp"


/*-------------------------------------------------------------------

  weather -- matrix with num_days_year rows and 4 columns: tmax
  (col1), tmin (col 2), srad (col3), prcp (col4) irri (col5)

  genetic_paras -- array with

  [0] amx (amx) = area largest leaf (800 cm2),
  [1] rue (rue_max) = radiation use efficiency (1.6)
  [2] tln (tlno) = total number leafs (19)
  [3] bkp (vpd_bkp) = vapor pressure deficit breakpoint (2)
  [4] slp (vpd-slp) = vapor pressure deficit slope (0.3)
  [5] shd (co_shd) = Thermal time between end of leaf expansion and shedding (20)
  [6] alf (co_a_lf) = leaf appearance rate parameter 1 (2.5)
  [7] blf (co_b_lf) = leaf appearance rate parameter 2 (0.00265)
  [8] ser (co_maxSer) = maximum silk elongation rate (1.5)
  [9] nrg (co_nrings) = number of rings per ear (40)
  [10] krg (co_kerRing) = kernels per ring (18)
  [11] hul (husk_length) = husk length (200 mm)
  [12] cs1 (co_ser[0]) = silk elongation rate parameter 1 (-0.155)
  [13] cs2 (co_ser[1]) = silk elongation rate parameter 2 (3.28)
  [14] cs3 (co_ser[2]) = silk elongation rate parameter 3 (0.56)
  [15] cs4 (co_ser[3]) = silk elongation rate parameter 4 (0.16)
  [16] kkw (kkw) = kernel weight flex (1.2)
  [17] sdv (swdf_veg) = supl. demand water deficit fact: sensitivity at veg state ([0,1] 0.67)
  [18] sds (swdf_sens) = supl. demand water deficit fact: sensitivity at senescent state ([0,1] 0.67)
  [19] sid (silkInidd) = number of silks initiated daily (60)

  env_paras -- 
  [0]  num_days_year (nrow weather), 
  [1]  latitute
  [2]  planting date,
  [3]  ppop
  [4]  mulch
  [5]  max length of growing season
  [6]  s_sat
  [7]  s_dul 
  [8]  s_ll
  [9] s_kl
  [10] s_layers
  [11] dlayer
  [12] swcon
  [13] isw

  out_end -- vector with state of many variables at termination of simulation

*/

// model_selection[0] - silking model (0 = empirical, 1 = process based)

void maizeGrainYieldModel(double *weather,
			  double *env_paras,
			  double *genetic_paras,
			  double *out_end
			  ){
			  
  //--------------------------------------------------
  // counters and indexes
  
  // day of the year
  int dofy = 0;
  
  // growing degree days
  double gdu = 0;

  // Simulation year
  int count_daily = 1;


  int h; // counter for loops going over hours 

  int z; // counter for loops going over soil layers

  int i = 0; // multi purpose counter

  //----------------------------------------------------------
  // misc. environment parameters

  int num_days_year = (int) env_paras[0];
  double latitude =  env_paras[1];
  int planting_date =  (int) env_paras[2];

  // plant population
  double plant_population = env_paras[3];

  double mulch_effect = env_paras[4];


  // terminate simulation if longer than that many days
  int count_daily_max =  (int) env_paras[5];

  //-----------------------------------------------------------
  // Soil

  // thermal time for emergence (Genetic & management-sowing depth)
  double eme_in = -87;
    
  // soil layer depth in environment
  double dlayr = env_paras[11];

  // saturation
  double s_sat = env_paras[6] * dlayr;

  // drain upper limitation
  double s_dul = env_paras[7] * dlayr;

  // lower limitation
  double s_ll = env_paras[8] * dlayr;
  
  double lai_fl =0;


  // drainage coefficient (swc-s.dul)*swcon -- daily water drainage if
  // sat > swc > dul
  double s_swcon = env_paras[12];

  // constant for all layers -- Valid for Viluco only. Dardanelli et
  // al. 1997
  double s_kl = env_paras[9];

  // number of soil layers
  int s_layers = (int) env_paras[10];

  int nsl;
  double twuse = 0;
  double tswc = 0;

  // initial soil water in fraction of total
  double pct_isw = env_paras[13];

  double insw = s_ll + pct_isw * (s_dul - s_ll);

  // vector soil water content (see s.sat, s.dul, s.ll multiplied by
  // dlayr so in mm)
  double sw[s_layers];
  
  // soil water holding capacity cm3/cm3. Use a function of depth to
  // account for variation in LL?
  double s_whc = (s_dul - s_ll) / dlayr;

  // soil depth mm 
  double s_depth = dlayr * s_layers;

  // initialize top soil water
  double esw = 0.0;

  //----------------------------------------------------
  // genotype specific plant parameters
  double amx = genetic_paras[0];  
  double rue_max = genetic_paras[1];  
  double tlno = genetic_paras[2];


  
  // limited transpiration breakpoint (kPa)
  double vpd_bkp = genetic_paras[3];

  // limited transpiration relative slope (1: no limitation; 0: full
  // limitation as in Sinclair et al & Messina et al
  double vpd_slope = genetic_paras[4];

  // Thermal time between end of leaf expansion and shedding -- 67 DD
  //double co_shd = 20;
  double co_shd = genetic_paras[5];

  //----------------------------------------------------
  // genetic plant parameters treated as constants

  // base temperature  emergence
  double tbeme = 8.0;

  // base temperature root growth
  double tbroot = 8.0;

  // base temperature seed growth
  double tbrep = 0.0;		    

  // base temperature vegetative growth
  double  tbv = 8;				       

  // coefficient of extinction -- may use sin(G) function.
  double coefext = 0.4;

  // te coefficient -- Muchow, Sinclair, Bennett...
  double tecoeff = 0.09;

  // Days between silking and onset of increase HI --also between
  // silking and kernel number determination
  // Thermal time for use in earbio model (3 days * 22 dd per day)
  double co_lag = 66; 
  double co_lagphase = 170; // SM

  // Thermal time duration grain fill (tbrep) Muchow & Sinclair
  double gfill_dur = genetic_paras[20];

  // rooting rate cm/d at optimal conditions (maybe genotype specific)
  double rootmax = 22.0;

  // leaf appearance rate (maybe genotype specific, ASI, SHED)
  //double co_A_lf = 2.5;
  double co_A_lf = genetic_paras[6];

  // leaf appearance rate (maybe genotype specific, ASI, SHED) was
  // 0.00225 !!!!!! confirm !!!!!!
  //double co_B_lf = 0.00265;
  double co_B_lf = genetic_paras[7];

  // input from Muchow & Sinclair from Carberry
  double co_A_ps = 0.00161;

  // psen (senescense) at physiological maturity
  double psen_pm = 0.8;

  // supply demand water deficit factor: sensitivity to senescense 
  double swdf_sens = genetic_paras[18];

  // supply demand water deficit factor: sensitivity [0,1] 0.67 leaf
  // expansion Jones & Kiniry. 1 photosynthesis
  double co_swdf_veg = genetic_paras[17];
  double co_swdf_ear = 1;
  
  // SM constants
  //double co_ser[4] = {-0.155, 3.28, 0.56, 0.16};
  double co_ser[4] = {genetic_paras[12], genetic_paras[13], genetic_paras[14], genetic_paras[15]};
  double co_maxSer = genetic_paras[8]; // mm/C
  double co_nrings = genetic_paras[9];
  double co_kerRing = genetic_paras[10];
  double co_silkear = co_nrings * co_kerRing;
  double silkInidd = genetic_paras[19]; //edmeades 1993 carcova 2003
  
  double husk_length = genetic_paras[11];
  double silkSen = 0.03;   // Anderson et al (2004)
  double silkSenThld = 350; // Anderson et al (2004)
  double alpha = 9.71 * 1000000; // max pollen number (model by Uribelarrea et al (200?))
  double beta = 2.66;            // spread pollen distribution (w days when pollen is 50% max)
  double beta_negsq = 1/(beta * beta);
  double tau = 3.62; // day of max pollen number (from begin pollination)

  // Keep this algorithm sequence. Calculate co.B.ps
  double tt_fl = (log(tlno) - log(co_A_lf)) / co_B_lf;
  double tt_total = tt_fl + co_shd + gfill_dur;
  double co_B_ps = (log(psen_pm) - log(co_A_ps)) / tt_total;

  //-----------------------------------------------------------
  // harvest index (HI) determination parameters

  // appearance rate -- fitted From Cooper et al. 2014
  double ks = 1.4;
  
  // maximum harvest index (maybe genotype specific)
   double maxHI = min2(0.45 + (0.58-0.45)/(990-588)*(co_silkear-588), 0.58); //tie maxHI with potential silk number
                                                                 // 588=42*14; 990=55*18
  // Daily increase in HI g m-2 d-1 -- note that it will be initialized
  // to maxHI /gfill.dur -- else overparameterized

  // Calculated so that maximum HI is only reached if the crop fully
  // matures -- consistent parameterization New implementation--
  // calculated as a function of maxHIa pre onset of grain fill.
  double dmac = maxHI / gfill_dur;

  // Kernel weight flex -- From Borras et al. 2004 FCR 86:131-146
  double kkw = genetic_paras[16];

  // leaf number at which onset of ear development is simulated. 
  double onset_egr = 14;
  double onset_egr_offset[3] = {-1.0,0.0,1.0};


  //--------------------------------------------------
  // potential ear growth preflowering

  // Thermal time at end of leaf expansion (pre-shed)
  tt_fl = (log(tlno) - log(co_A_lf)) / co_B_lf;

  // Thermal time at V stage
  double tt_v = (log(onset_egr) - log(co_A_lf)) / co_B_lf;

  // Thermal time V stage to onset of HI increase.
  double tt_hi = tt_fl + co_shd + co_lag - tt_v;

  // potential ear growth rate from Vn to 66 dd after shedding (g
  // °Cd-1). Target mass (3) from Cooper et al AJAR Review
  double pot_egr = (log(3) - log(0.01)) / tt_hi;

  //----------------------------------------------------
  // Initialize variables
  //----------------------------------------------------

  // Weather
  double tmax = 0.0;
  double tmin = 0.0;
  double srad = 0.0;
  double prcp = 0.0;
  double irri = 0.0;
  
  double tmean = 0.0;
  double dleng = 0.0;

  // hourly relative solar radiation 
  double r_srad[24] = { 0 }; 

  // hourly vapor pressure deficit
  double vpd_h[24] = { 0 };

  // hourly temperature (?)
  double tem_h[24] = { 0 };

  // hourly rue max (?)
  double rue_max_h[24] = { 0 };

  // diurnal mean vpd 
  double w_vpd = 0.0;


  // Development
  double d_eme = 0.0;
  double d_veg = 0.0;
  double d_seed = 0.0;

  double eme = eme_in;
  double dev_state_termleaf = 0.0;
  double dev_state_ff = 0.0;
  double dev_state_laghi = 0.0;
  double dev_state_seed = 0.0;
  double dev_root = 200.0;       // root depth mm
  double dev_state_ear = 0.0;
  double dev_state_ear_m1 = 0;
  double dev_state_ear_p1 = 0;
  double dev_state_pollen = 0.0; // from SM
  double dev_state_lag = 0.0; // from SM

  double d_root = 0.0;

  // all Boolean
  int dev_stage_eme = 0;
  int dev_stage_termleaf = 0;
  int dev_stage_ff = 0;
  int silk_set = 0;
  int dev_stage_bhi = 0;
  int dev_stage_pm = 0;
  int maxHIa_calc = 0;
  int dev_stage_lag = 0;  // from SM
		
  // hardcoded to 3 plants for now!
  int num_ears = 3;
  double ear_contributions[num_ears];
  ear_contributions[0]= 0.4;
  ear_contributions[1]= 0.4;
  ear_contributions[2]= 0.2;
  
  int silk_day = 0;
  int shd_day = 0;

  double gdu_shd, gdu_silk;
			
  // Leaf
  double leaf_number = 1.0;
  double pla = 0.0;
  double dila = 0.0;
  double pla_wde = 0.0;
  double tlai = 0.0;
  double sens_lai = 0.0;
  double lai = 0.0;
  double dev_sen = 0.0;
  double tlai_tm1 = 0.0;
  double pla_tm1 = 0.0;
  double dlai = 0.0;
  double psen = 0.0;
  double sens_lai_previous = 0.0;
  double sens_lai_daily = 0.0;
  double psen_water = 0.0;
  double sens_lai_water = 0.0;
  
  
  
  // Mass (Vegetative & Reproductive)
  double mass_day = 0.0;
  double veg_mass = 1.0;
  double ear_count = 0; // not used
  double actEgr = 0; // not used
  double potEbio = 0.0;
  double potEbio_m1 = 0.0;
  double potEbio_tm1 = 0.0;
  double potEbio_tm1_m1 = 0.0;
  double potEbio_p1 = 0.0;
  double potEbio_tm1_p1 = 0.0;
  double d_ear = 0.0;
  double d_ear_m1 = 0.0;
  double d_ear_p1 = 0.0;

  //--------------------------------------------------
  // silk growth 

  double dbp = 0; //days from begining of pollination
  double dpol_silk = 0; // amount of pollen per silk
  double dpollen[num_ears]; // pollen grains per plant
  double fert_factor[num_ears];
  memset(dpollen, 0, num_ears * sizeof(double));
  memset(fert_factor, 0, num_ears * sizeof(double));

  // number, length, age and receptivity of silks -- three plants, maximum 50 cohorts
  // (plants in rows)
  double silk_num[num_ears * 50]; 
  double silk_len[num_ears * 50];
  double silk_age[num_ears * 50];
  double silk_rec[num_ears * 50];

  memset(silk_num, 0, num_ears * 50 * sizeof(double));
  memset(silk_len, 0, num_ears * 50 * sizeof(double));
  memset(silk_age, 0, num_ears * 50 * sizeof(double));
  memset(silk_rec, 0, num_ears * 50 * sizeof(double));

  // number of silk cohorts (each day in which silks are initiated is a cohort)
  int silk_cohort[num_ears];
  memset(silk_cohort, 0, num_ears * sizeof(int));  

  // total number of silks
  double silk_num_total[num_ears];
  memset(silk_num_total, 0, num_ears * sizeof(double));  

  // temp value used for calculating pollen availability
  double silk_num_pollen = 0;

  // number of embryos (three plants, max 50 cohorts, plants in rows)
  double embryo_num[num_ears * 50];
  memset(embryo_num, 0, num_ears * 50 * sizeof(double));
  double embryo_num_total = 0; // total number of embryos
  double silk_num_total_all = 0; //total number of silks
  
  // daily increments of initiated silks, silk length
  double dsilk_num[num_ears];
  memset(dsilk_num, 0, num_ears * sizeof(double));
  double dsilk_len[num_ears];
  memset(dsilk_len, 0, num_ears * sizeof(double));

  /* double ebio[num_ears] = {0.0}; */
  /* memset(dsilk_len, 0, num_ears * sizeof(double)); */
  /* double ebio_m1 = 0.0; */
  /* double ebio_p1 = 0.0; */

  double SNa = 1.0;
  double SNa_m1 = 0.0;
  double SNa_p1 = 0.0;
  double SNa_0 = 0.0;
  double gg = 1.0;
  double maxHIa = maxHI;
  double seed_mass = 0.0;
  double seed_growth_rate = 0.0;
  double hi = 0.001;           // initial harvest index
  double d_laghi = 0.0;


  //----------------------------------------
  // Water
    
  double w_demand = 1;
  double w_supply = 1;
  double swdf_veg = 1.0;
  double swdf_ear = 0;


  // dev.root (root depth) in cm x 10 to convert to mm
  double ttsw = s_whc * dlayr * s_layers;  

  // dev.root in cm x 10 to convert to mm
  double atsw = ttsw * pct_isw;

  // dev.root in cm x 10 to convert to mm
  double ftsw = atsw / ttsw; 
  double ftsw2 = 0;


  // days from last precipitation or irrigation event          
  double dyse = 1.0;
  double soil_evp = 0.0;
  double transp = 0.0;
  double transp_layer = 0.0;

  double ftsw_shd = 0;

 
  double pot_rwu[s_layers];
  double rwu[s_layers];

  double drain[s_layers];
  double hold[s_layers];
  double flux[s_layers + 1];

  
  // holds output from carbonWater function. Element 1 is
  // photosynthesis, element 2 the transpiration, element 3 the water
  // demand
  double phs_transp_demand[3]; 



  for(z = 0; z < s_layers; ++z){
    sw[z] = insw;
    pot_rwu[z] = 0;
    rwu[z] = 0;
    drain[z] = 0;
    hold[z] = 0;
    flux[z] = 0;
  }

  flux[s_layers] = 0;

  //------------------------------------------------------
  // loop over days within year
    
  while(1){	
    

    // day of year
    dofy = planting_date + count_daily;

    //-----------------------------------------------------
    // Weather calculations and stress factors

    // get weather inputs for day of year
    tmax = weather[count_daily];
    tmin = weather[(count_daily) + num_days_year];
    srad = weather[(count_daily) + num_days_year * 2];
    prcp = weather[(count_daily) + num_days_year * 3];
    irri = weather[(count_daily) + num_days_year * 4];
  
    tmean = (tmax + tmin) / 2; 
    dleng = dayLength(latitude, dofy);
    hourlyVPD(tmax, tmin, vpd_h, tem_h);
    hourlyRueMax(rue_max, tem_h, rue_max_h);
    relSolRad(dleng, r_srad);

    gdu += dGDU(tmax, tmin);

    // diurnal mean vpd 
    w_vpd = (vpd_h[max(vpd_h, 24)] - vpd_h[min(vpd_h, 24)]) * 0.75;

    //---------------------------------------------------------
    // Development
        
    // Development stage planting, emergence, end leaf expansion
    // and onset seed growth (-- is seed growth when HI start
    // growing)
    if(eme >= 0) dev_stage_eme = 1;
    if(leaf_number >= tlno) dev_stage_termleaf = 1;
    if(dev_state_ff >= co_shd) dev_stage_ff = 1;
    if(dev_state_lag >= co_lagphase) dev_stage_lag = 1;
    if(dev_stage_lag == 1) dev_stage_bhi = 1;
    if(dev_state_seed >= gfill_dur) dev_stage_pm = 1;
    
    // Development rate
    d_eme = devRateEme(tmean, tbeme);
    
    if(dev_stage_eme == 1) d_veg = devRateVeg(tmean, tbv); // function maize
    if(dev_stage_bhi == 1) d_seed = devRateSeed(tmean,tbrep);
    d_root = devRateRoot(rootmax, tmax, tmin, tbroot);

    // or one day if want to use Muchow & Sinclair 1 day.
    d_laghi = devRateVeg(tmean, tbv);
    

    //---------------------------------------
    // Development state
    
    // Emergence
    if(dev_stage_eme == 0) eme += d_eme;
        
    // terminal leaf expansion (cultivar specific)
    if(dev_stage_eme == 1 && dev_stage_termleaf == 0)
      dev_state_termleaf += d_veg;


    // Shedding
    if(dev_stage_termleaf == 1 && dev_stage_ff == 0)
      dev_state_ff += d_veg;
       
    // shedding day
    if(dev_state_ff >= co_shd && dev_stage_ff == 0){
      shd_day = dofy;
      gdu_shd = gdu - dGDU(tmax, tmin);
    }
 
    // lag phase - shedding to onset HI or seed growth
    if(dev_stage_ff == 1 && dev_stage_bhi == 0)
      dev_state_laghi += d_laghi;

    // seed growth
    if(dev_stage_bhi == 1) 
      dev_state_seed =  min2(gfill_dur, dev_state_seed + d_seed);

    // Root depth
    if(dev_stage_bhi == 0 && dev_root <= s_depth)
      dev_root = min2(dev_root + d_root, s_depth); // root depth in cm
    
    // senescence
    if(dev_stage_eme == 1)
      dev_sen += d_veg;

    /* // ear growth */
    /* if(leaf_number >= onset_egr && dev_stage_bhi == 0) */
    /*   dev_state_ear = min2(tt_hi, dev_state_ear + d_veg); */

    /* if(leaf_number >= onset_egr_m1 && dev_stage_bhi == 0) */
    /*   dev_state_ear_m1 = min2(tt_hi, dev_state_ear_m1 + d_veg); */

    /* if(leaf_number >= onset_egr_p1 && dev_stage_bhi == 0) */
    /*   dev_state_ear_p1 = min2(tt_hi, dev_state_ear_p1 + d_veg); */

    // leaf Area vegetative stages				
    if(dev_stage_eme == 1 && dev_stage_termleaf == 0){
      tlai_tm1 = tlai;              // tlai t-1 total lai
      pla_tm1 = pla;

      // new function here based on cumulative thermal time and
      // 2 geno coef
      leaf_number = min2(tlno, leafNum(dev_state_termleaf, co_A_lf, co_B_lf)); 

      pla = leafArea(leaf_number, amx, tlno);
      dila = pla - pla_tm1;

      // Jones & Kiniry / Ritchie stress function	
      pla_wde = max2(0, dila * swdf_veg);	
      dlai = plant_population * pla_wde / 10000;
      tlai = tlai_tm1 + dlai;     
      lai = lai + dlai;
    }

    // proportion of senescence relative to total lai Muchow &
    // Carberry (1989) w/param f(CRM)
    sens_lai_previous = sens_lai;
    psen = co_A_ps * exp(co_B_ps * dev_sen);
    sens_lai = tlai * psen;    
    sens_lai_daily = sens_lai - sens_lai_previous;
            
    psen_water = swdf_sens*(1-min2(1,w_supply/w_demand));   //proportion of senescence relative to swdf
    sens_lai_water = tlai * psen_water;
    
    lai = max2(0.0, lai - max2(sens_lai_daily, sens_lai_water));
    
    //Rprintf("lai = %f s.d.ratio= %f sens_lai= %f\n", lai,w_supply/w_demand,max2(sens_lai_daily, sens_lai_water)); 
    
    
    //--------------------------------------------------
    // Growth


    // Vegetative growth	
    if(dev_stage_eme == 1 && dev_stage_bhi == 0){
            
      carbonWater(coefext, rue_max_h, lai, srad, r_srad,
		  vpd_h, dev_state_seed, tecoeff, vpd_bkp, 
		  vpd_slope, w_supply, phs_transp_demand);
	    
      veg_mass += phs_transp_demand[0];
      
    }

    // SM BEGIN
	
	// Here is the calculation of pollen availability. Note
	// that SHD refers to the population of plants and the
	// calculation here is within the first cohort of
	// plants. This is necessary to have the calculation of
	// pollen available starting from the first plants silking
	// and shedding pollen.  update pollen availability later
	// when calculating pollen per silk for other cohorts
	// (substract 2 pollen per silk from pollen available)
	if(dev_stage_ff == 1 && dev_stage_lag == 0){
		dev_state_pollen += d_veg;
		// days from begin pollination
		dbp = dev_state_pollen / (2 * co_shd) * tau;
		// pollen grains per plant -- Uribelarrea et al. 2004
   //	dpollen[0] = alpha * exp(-2 * pow(dbp - tau, 2) * beta_negsq) / (beta * sqrt(M_PI_2));
   	dpollen[0] = alpha * exp(-2 * pow(dbp - tau, 2) * beta_negsq) / (beta * sqrt(M_PI_2));

	}	
	
	//Rprintf("dev_stage_ff = %f\n", dev_stage_ff); 
	//Rprintf("dev_stage_lag = %f\n", dev_stage_lag); 
	//Rprintf("leaf_number = %f\n", leaf_number); 
	//Rprintf("dev_state_ff = %f\n", dev_state_ff); 
	//Rprintf("dev_state_lag = %f\n", dev_state_lag); 
	
	
	
	
	

    // Reproductive growth pre-onset seed growth. Three ears are
    // simulated to represent a plot. Each initiated at leaf.number >=
    // onset.egr - 1, leaf.number >= onset.egr and leaf.number >=
    // onset.egr + 1 to cover about 6 days of silking
 
    for(int ear = 0; ear < num_ears; ++ear){
		if(leaf_number >= (onset_egr + onset_egr_offset[ear]) && dev_stage_lag == 0){
	
			dsilk_num[ear] = 0;
			if(silk_num_total[ear] < co_silkear){
	  
				++silk_cohort[ear];   // new silk cohort	    
				dsilk_num[ear] = min2(silkInidd * d_veg, co_silkear - silk_num_total[ear]);
				silk_num[(silk_cohort[ear]-1) * num_ears + ear] = dsilk_num[ear];
				/* silk_len[(silk_cohort[ear]-1) * num_ears + ear] = 0; */
				/* silk_age[(silk_cohort[ear]-1) * num_ears + ear] = 0; */

				silk_num_total[ear] += dsilk_num[ear];
			}
			
			// elongate silks by cohort
			dsilk_len[ear] = co_maxSer * d_veg * ser(co_ser, ftsw);
			
			for(int co = 0; co < silk_cohort[ear]; ++co){
				silk_len[co * num_ears + ear] += dsilk_len[ear];
			}

			// elongate and age silks by cohort after emergence (= first cohort silks out)
			if(silk_len[0 * num_ears + ear] >= husk_length){
			        
			        
				// record silk date of first ear (silking mid-point)
				if(ear == 0 && silk_day <= 0){ 
					silk_day = dofy;
					gdu_silk = gdu - dGDU(tmax, tmin);
				}

				// accumulate degree days towards lagphase (only for second ear = mid silking)
				if(ear == 1){
					dev_state_lag += d_veg;
				}

				// determine pollination based on silk exposure, pollen and receptivity by cohort
				for(int co = 0; co < silk_cohort[ear]; ++co){

					if((silk_len[co * num_ears + ear] >= husk_length) && 
						(embryo_num[co * num_ears + ear] < silk_num[co * num_ears + ear])){
					        
					        

						silk_age[co * num_ears + ear] += d_veg; 
						silk_rec[co * num_ears + ear] = 1 / (1 + exp(silkSen * (silk_age[ear] - silkSenThld))); 
                  double d = exp(3.0);
		
						silk_num_pollen = max2(silk_num[co * num_ears + ear] - embryo_num[co * num_ears + ear], 0.001);

						// amount of pollen per silk
						dpol_silk = dpollen[ear] / silk_num_pollen;
						fert_factor[ear] = (dpol_silk >= 2) ? (1.0) : 1 - (2-dpol_silk)/2;
		
						embryo_num[co * num_ears + ear] = embryo_num[co * num_ears + ear] +
						silk_num_pollen * silk_rec[co * num_ears + ear] * fert_factor[ear];

						// update pollen to estimate available pollen for next silk cohort
						dpollen[ear] = max2(dpollen[ear] - 2 * silk_num_pollen, 0);
						
					} // if silks > husks && no embryos yet
				} // for cohorts				
			} // silk length > husk length
			
			// update pollen to estimate available pollen for next ear
			if(ear < (num_ears - 1))
				dpollen[ear + 1] = dpollen[ear];
				
		} // onset eargrowth 
    } // for ears

      
    if(dev_stage_lag == 1 && maxHIa_calc == 0){

    // numboer of silks from all ears and cohorts
    silk_num_total_all = 0;
    for(int ear = 0; ear < num_ears; ++ear){
                    silk_num_total_all += (silk_num_total[ear] * ear_contributions[ear]);
      }      
            
      // numboer of embryos from all ears and cohorts
      embryo_num_total = 0;
      for(int ear = 0; ear < num_ears; ++ear){
      	for(int co = 0; co < silk_cohort[ear]; ++co){
      	  embryo_num_total += (embryo_num[co * num_ears + ear] * ear_contributions[ear]);
      	}}
      	

      
      maxHIa = 0;
      if(embryo_num_total > 1){
	double gg1 = (((1/maxHI) - 1) * co_silkear / (embryo_num_total * kkw)) + 1;
	maxHIa = min2(maxHI,
		      1/gg1);
      }

      dmac = maxHIa / gfill_dur;
      maxHIa_calc = 1;
      
    }       // SM END
  
    // Reproductive growth seed fill

    if(dev_stage_lag == 1){
            if(lai_fl==0){
                    lai_fl = lai;
                    //Rprintf("lai %f lai_fl\n", lai,lai_fl);

            }

      carbonWater(coefext, rue_max_h, lai, srad, r_srad,
		  vpd_h, dev_state_seed, tecoeff, vpd_bkp, 
		  vpd_slope, w_supply, phs_transp_demand);
      
					
      for(h = 1; h <= 24; ++h){

	seed_growth_rate = max2(0, ((seed_mass / hi) - veg_mass - seed_mass) / (1.0 - 1.0 / hi));
	veg_mass += ((phs_transp_demand[0] / 24) - seed_growth_rate);
	seed_mass += seed_growth_rate;

	// recall dmac is in dd -- multiply by d.seed (today's
	// degree days to get dmac today's in g m-2 d-1
	hi += (dmac / 24 * d_seed);
                
	if(hi >= maxHIa){
	  hi = maxHIa;
	  break;  // terminate simulation
	}
      }  // close hourly loop
    } // close reproductive if



    //--------------------------------------------------------
    // Soil water balance

			
    // Transpiration
    if(dev_stage_eme == 1)
      transp = phs_transp_demand[1];
    else
      transp = 0.0;

    // Infiltration (pending implementation of USDA CN method)
				
    // Drainage Jones & Kiniry
    
    for(z = 0; z < s_layers; ++z){
     drain[z] = 0;
      hold[z] = 0;
      flux[z] = 0;
      rwu[z] = 0;
    }
    flux[s_layers] = 0;

    for(z = 0; z < s_layers; ++z){

      if(z == 0)
	flux[z] = irri + prcp;
      
      hold[z] = (s_sat - sw[z]);

      if(flux[z] >= hold[z]){
	sw[z] = s_sat;
	drain[z] = (sw[z] - s_dul) * s_swcon;
	sw[z] -= drain[z];
	flux[z + 1] = flux[z] - hold[z] + drain[z];
      }
      else{
	sw[z] += flux[z];
	
	if(sw[z] <= s_dul) 
	  drain[z] =  0;
	else 
	  drain[z] = (sw[z] - s_dul) * s_swcon;
      
	sw[z] -=  drain[z];
	flux[z+1] = drain[z];
      }


    }  // close loop drainage by soil layer
  
    //-------------------------------------------------------
    // calculate Evaporation
    esw = sw[0] - s_ll;

    if(prcp >= 1.5 || irri >= 1.5)
      dyse = 1.0;
    else
      dyse = dyse + 1.0;
  
    soil_evp = evaporation(dyse, tmax, tmin, srad, prcp, w_vpd, lai, esw, mulch_effect);
  

    //-------------------------------------------------------
    // Update water content (mm) based on evaporation

    sw[0] -= soil_evp;

    //  number or soil layers with active roots. Root water uptake
    // already account for soil depth occupancy
    nsl = min2(s_layers,ceil(dev_root / dlayr));  

    // calculate ppotential Root water uptake --Supply free.water<-
    // if(sw>dul) w-dul else 0 . dW/dt = kl * (W -Wll) + free.water
  
    for(z = 0; z < nsl; ++z){

      if(sw[z] <= s_dul)
	pot_rwu[z] = (sw[z] - s_ll) * s_kl;
      else 
	pot_rwu[z] = (s_dul - s_ll) * s_kl + sw[z] - s_dul;

      if(z == (nsl - 1))
	pot_rwu[z] = pot_rwu[z] * (dev_root - (nsl-1)*dlayr) / dlayr;
    }

    //-------------------------------------------------------
    // Update water content (mm) based on Transpiration,

    transp_layer = transp;
    for(z = 0; z < nsl; ++z){

      rwu[z] = min2(pot_rwu[z],transp_layer);

      sw[z] -= rwu[z];
      
      transp_layer = max2(0, transp_layer - rwu[z]);
      twuse = max2(0,twuse - rwu[z]);

    }  // close update water content by layer

    // total soil water content
    tswc = 0.0;
    for(z = 0; z < s_layers; ++z)
      tswc += sw[z];
  
    // total water use
    twuse = transp + soil_evp;
	  
    w_supply = 0;
    for(z = 0; z < s_layers; ++z)
      w_supply +=  pot_rwu[z];


    // informative metrics
    if(dev_root <= s_depth){

      // d.root (increase in soil/root depth in cm) x 10 to
      // convert to mm
      ttsw += d_root * s_whc;

      // new available water from soil exploration/root growth
      // and given initial condition
      atsw += d_root * s_whc * pct_isw;
    }

    // precipitation and runoff
    atsw = min2(ttsw + 4, atsw + prcp + irri);

    // subtracts evaporation and transpiration
    atsw = max2(0.0, atsw - soil_evp - transp); 
    ftsw = max2(0.0, atsw / ttsw);
    //ftsw2 = w_supply/ttsw;
    //Rprintf("ftsw %f ftsw2 %f\n", ftsw, ftsw2);
    //Rprintf("dev_root %f s_depth %f\n", dev_root, s_depth);
    
    //----------------------------------------
    // stress factors

    w_demand = (dev_stage_eme == 1) ? (phs_transp_demand[2]) : 0.1;

    (dev_stage_eme == 1) && (swdf_veg = max2(0, min2(1,w_supply/w_demand * co_swdf_veg))); 
    (dev_stage_eme == 1) && (swdf_ear = max2(0, min2(1,w_supply/w_demand * co_swdf_ear))); 
    //Rprintf("vegmass = %f\n", veg_mass);
    
    //------------------------------------------
    // Simulation termination control
      

    if(dev_stage_pm == 1 || 
       hi >= maxHIa ||
       count_daily >= count_daily_max){
    
      //------------------------------------------
      // Output
      out_end[0] = veg_mass;
      out_end[1] = seed_mass;
      out_end[2] = gdu_silk;
      out_end[3] = gdu_shd;
      out_end[4] = count_daily;
      out_end[5] = silk_num_total_all;
      out_end[6] = embryo_num_total;
      out_end[7] = hi;
      out_end[8] = maxHIa;
      out_end[9] = lai_fl;
      

      // terminate simulation
      break;

    }
  
    // a new day
    ++count_daily;

    // check for an interrupt 
    //R_CheckUserInterrupt();

  }  // end daily loop
}

// EOF
