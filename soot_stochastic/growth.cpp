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


void linerarSurfGrowth(vector<vector<double> > &allParticles, double deltaM0, double lp0, double maxValL, double deltaL)
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
        if(allParticles[i][1]<(maxValL+deltaL/2-deltaMi))
        {
            allParticles[i][1] += deltaMi;
        }
        
    }
    cout << endl;
}

void linerarSurfOxi(vector<vector<double> > &allParticles, double deltaM0, double lp0, double deltaL)
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
        if(allParticles[i][1]>(lp0-deltaL/2-deltaMi))
        {
            allParticles[i][1] += deltaMi;
        }
        
    }
    cout << endl;
}



void surfGrowthAging(vector<vector<double> > &allParticles, double deltaM0, double lp0, double maxValL, double deltaL, double ageFactor)
{
    double pid = 3.14159265;
    double v0 = lp0;
    double surf0 = pid*pow(v0*6/pid,0.6666);
    
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        double surfReactivityi(1);
        surfReactivityi = 1 - ageFactor * allParticles[i][1] / maxValL; // if ageFactor = 1, surf reactivity ranges from 1 - lp0 / maxValL (0.99 for l = lp0) and 0 for l = maxVal. 0.5 for l = maxVal / 2
        
        double vi = allParticles[i][1];      // li is actually homogeneous to a volume
        double surfi = pid*pow(vi*6/pid,0.6666);
        double deltaMi = surfi / surf0 * deltaM0 * surfReactivityi;  // effect of aging added
        if(allParticles[i][1]<(maxValL+deltaL/2-deltaMi))
        {
            allParticles[i][1] += deltaMi;
        }
        
    }
    cout << endl;
}
