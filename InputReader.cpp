#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <sstream>
using namespace std;

int ReadMet(string filename, double *weather, int nDays)
{
    ifstream file(filename);
    string str;
    const char s[2] = ",";
    char *token;
    int indx = 0;
    double val;

    while (getline(file, str))
    {
        /* get the first token */
      
        int col = 0;
        const char* ptr = str.c_str();
        token = strtok((char *)ptr, s);

        /* walk through other tokens */
        while (token != NULL)
        {
            val = stod(token);
            weather[indx + col * nDays] = val;
            token = strtok(NULL, s);
            col++;
        }
        indx++;
    }

    return 0;
}

void GetGeneticParams(double *genetic_paras)
{
    genetic_paras[0] = 800;     //amx (amx) = area largest leaf (800 cm2),
    genetic_paras[1] = 1.6;     //rue (rue_max) = radiation use efficiency (1.6)
    genetic_paras[2] = 19;      //tln (tlno) = total number leafs (19)
    genetic_paras[3] = 2.0;     //bkp (vpd_bkp) = vapor pressure deficit breakpoint (2)
    genetic_paras[4] = 0.3;     //slp (vpd-slp) = vapor pressure deficit slope (0.3)
    genetic_paras[5] = 20.0;    //shd (co_shd) = Thermal time between end of leaf expansion and shedding (20)
    genetic_paras[6] = 2.5;     // alf (co_a_lf) = leaf appearance rate parameter 1 (2.5)
    genetic_paras[7] = 0.00265; //blf (co_b_lf) = leaf appearance rate parameter 2 (0.00265)
    genetic_paras[8] = 1.5;     //ser (co_maxSer) = maximum silk elongation rate (1.5)
    genetic_paras[9] = 40.0;    //nrg (co_nrings) = number of rings per ear (40)
    genetic_paras[10] = 18;     //krg (co_kerRing) = kernels per ring (18)
    genetic_paras[11] = 200;    //hul (husk_length) = husk length (200 mm)
    genetic_paras[12] = -0.155; // cs1 (co_ser[0]) = silk elongation rate parameter 1 (-0.155)
    genetic_paras[13] = 3.28;   //cs2 (co_ser[1]) = silk elongation rate parameter 2 (3.28)
    genetic_paras[14] = 0.56;   //cs3 (co_ser[2]) = silk elongation rate parameter 3 (0.56)
    genetic_paras[15] = 0.16;   //cs4 (co_ser[3]) = silk elongation rate parameter 4 (0.16)
    genetic_paras[16] = 1.2;    //kkw (kkw) = kernel weight flex (1.2)
    genetic_paras[17] = 0.67;   //sdv (swdf_veg) = supl. demand water deficit fact: sensitivity at veg state ([0,1] 0.67)
    genetic_paras[18] = 0.67;   //sds (swdf_sens) = supl. demand water deficit fact: sensitivity at senescent state ([0,1] 0.67)
    genetic_paras[19] = 60.0;   //sid (silkInidd) = number of silks initiated daily (60)
    genetic_paras[20] = 1300.0;   //gfd (tbrep) = // Thermal time duration grain fill (tbrep) Muchow & Sinclair
}

void     GetEnvironParams(double *env_paras)
{
    // ENV1     
    env_paras[0] = 240;          // num_days_year (nrow weather),
    env_paras[1] = 42.79374482;  // latitute
    env_paras[2] = 125;          //planting date,
    env_paras[3] = 8.14;         // ppop
    env_paras[4] = 0;            //mulch
    env_paras[5] = 140;          //max length of growing season
    env_paras[6] = 0.35770744;   //s_sat
    env_paras[7] = 0.285193965;  //s_dul
    env_paras[8] = 0.134497798;  //s_ll
    env_paras[9] = 0.08;         // s_kl
    env_paras[10] = 15;          //s_layers
    env_paras[11] = 100;         //dlayer
    env_paras[12] = 0.31;        //swcon
    env_paras[13] = 0.810749022; // isw
}

