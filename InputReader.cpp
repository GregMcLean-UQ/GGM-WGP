#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

using namespace std;

int ReadMet(string filename, double *weather)
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
      

       const char* ptr = str.c_str();
        token = strtok((char *)ptr, s);

        /* walk through other tokens */
        while (token != NULL)
        {
            val = stod(token);
            weather[indx++] = val;
            token = strtok(NULL, s);
        }
    }

    return 0;
}

void GetGeneticParas(double *genetic_paras)
{

}
