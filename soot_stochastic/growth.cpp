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

void uniformGrowth(vector<vector<double> > &allParticles, double deltaG, double timePerIt)
{
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        allParticles[i][1] += deltaG*timePerIt;
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



double outOfBoundSource(vector<vector<double> >const& allParticles, double nT, double maxValL, double lp0, double deltaL)
{
    double dotG(0);
    int countNp(0);
    int countOut(0);
    int i(0);
    for(i=0; i< allParticles.size(); i++)
    {
        countNp++;
        if((allParticles[i][1] >= (maxValL + deltaL/2))||(allParticles[i][1] < (lp0 - deltaL/2)))
        {
            countOut++;
        }
    }
    double Np = countNp;
    double out = countOut;
    double ratioOut = out/ Np;
    dotG = -ratioOut * nT;                  // calcul of dotH, dotG, dotAt, before advancing t.
    return dotG;
}


void advanceGrowthPdf(vector<vector<double> >& allParticles, double nT, double maxValL, double lp0, double deltaL, vector<vector<double> >const& lAndNpL)

{
    //calculation of ratioOut and construction of ranksToReallocate vector
    int countNp(0);
    int countOut(0);
    vector<int> ranksToReallocate;
    
    int i(0);
    for(i=0; i< allParticles.size(); i++)
    {
        countNp++;
        if((allParticles[i][1] >= (maxValL + deltaL/2))||(allParticles[i][1] < (lp0 - deltaL/2)))
        {
            countOut++;
            ranksToReallocate.push_back(i);
        }
    }
    cout << "total number of particles = " << countNp << endl;
    double Np = countNp;
    double out = countOut;
    double ratioOut = out/ Np;
   
    
    //constructing vector with the deltaNpl values for each li of the domain
    vector<int> deltaNpInt;
    i=0;
    for(i=0; i<lAndNpL.size(); i++)
    {
        double npl1(0);
        double npl2(0);
        npl1 = lAndNpL[i][1];
        npl2 = npl1 / (1-ratioOut);
        int rounded1 = rounding(npl1);
        int rounded2 = rounding(npl2);
        int deltaNpl = rounded2 - rounded1;
        deltaNpInt.push_back(deltaNpl);
    }
    
    
    // constructing vector of values to realloc
    vector<double> valuesToRealoc;
    
    i=0;
    double li(0);
    for(i=0; i<lAndNpL.size(); i++)
    {
        li = lAndNpL[i][0];
        int deltaNpli = deltaNpInt[i];
        
        if(deltaNpli>0)                  // for now we treat only this case. if deltaNpli < 0 it means that particles entered the domain by oxidation of bigger particles than maxValL. We are not considering this for now. For lower than lp0 only nucleation.
        {
            
            int j=0;
            for(j=0; j< deltaNpli; j++)
            {
                valuesToRealoc.push_back(li);
            }
        }
    }
    
    
    // count valuesToRealoc
    i=0;
    int countValuesToRealoc(0);
    for(i=0; i<valuesToRealoc.size(); i++)
    {
        countValuesToRealoc++;
    }
    cout << "growth valuesToRealoc = " << countValuesToRealoc << endl;
    cout << "growth ranksToRealoc = " << countOut << endl;

    // Now reallocating!
    i=0;
    int min = countValuesToRealoc;              // due to rounding error
    if(countOut < countValuesToRealoc)
        min = countOut;
    
    for(i=0; i<min; i++)
    {
        int rankParticle(0);
        rankParticle = ranksToReallocate[i];
        
        double newVal(0);
        newVal = frand_a_b(valuesToRealoc[i]-deltaL/2.0, valuesToRealoc[i]+deltaL/2.0);  // a random value in the interval of valuesToRealoc[i] is chosen to reset the particle
        
        allParticles[rankParticle][1] = newVal;
        
    }
    
}



void linearGrowth(vector<vector<double> > &allParticles, double timePerIt)
{
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        allParticles[i][1] += allParticles[i][1] * timePerIt;  // test case agglo+growth. G = x (size coordinate)
    }
}

