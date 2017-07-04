//
//  growth.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 03/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "growth.hpp"
#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

void uniformGrowth(vector<vector<double> > &allParticles, double deltaG)
{
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        allParticles[i][1] += deltaG;
    }
}


void linerarSurfGrowth(vector<vector<double> > &allParticles, double deltaM0, double lp0)
{
    double pid = 3.14159265;
    double v0 = lp0;
    double surf0 = pid*pow(v0*6/pid,0.6666);
    
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        double vi = allParticles[i][1];      // li is actually homogeneous to a volume
        double surfi = pid*pow(vi*6/pid,0.6666);
        double deltaMi = surfi / surf0 * deltaM0;
        allParticles[i][1] += deltaMi;
    }
    cout << endl;
}

void linerarSurfOxi(vector<vector<double> > &allParticles, double deltaM0, double lp0)
{
    double pid = 3.14159265;
    double v0 = lp0;
    double surf0 = pid*pow(v0*6/pid,0.6666);
    
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        double vi = allParticles[i][1];      // li is actually homogeneous to a volume
        double surfi = pid*pow(vi*6/pid,0.6666);
        double deltaMi = surfi / surf0 * deltaM0;
        allParticles[i][1] += deltaMi;
    }
    cout << endl;
}

double dotOxi(vector<vector<double> > const& allParticles, double lp0, double deltaL)
{
    int i(0);
    int countDotOxi(0);
    for(i=0; i<allParticles.size(); i++)
    {
        if(allParticles[i][1]<(lp0-deltaL/2))
            countDotOxi--;
    }
    return countDotOxi;
}

