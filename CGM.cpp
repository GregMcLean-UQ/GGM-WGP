// CGM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string.h>
#include "maize.cpp"
#include "InputReader.cpp"

using namespace std;
#include <chrono>
using namespace std::chrono;

int main()
{

    double env_paras[14];
    double genetic_paras[20];
    double weather[366 * 5];
    double out_end[10];

    string filename = ".\\Data\\metEnv1.csv";
    filename =  "C:\\Temp\\metEnv1.csv";

    // Read met data into weather.
    ReadMet(filename, weather);
    GetGeneticParams(genetic_paras);
    GetEnvironParams(env_paras);

    auto start = high_resolution_clock::now();
    maizeGrainYieldModel(weather, env_paras, genetic_paras, out_end);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << duration.count() << endl;


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
