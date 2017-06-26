//
//  agglomeration.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 26/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "agglomeration.hpp"
#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

int geometricNpLstar(double lStar, vector<vector<double> > allParticles)
{
    double infBorn = lStar*0.75;
    double supBorn = lStar*1.5;
    int i(0);
    int count(0);
    for(i=0; i<allParticles.size(); i++)
    {
        if((allParticles[i][1]>=infBorn)&(allParticles[i][1]<supBorn))
            count++;
    }
    return count;
}

double geometricNLstar(double lStar, vector<vector<double> > allParticles, double nT)
{
    double infBorn = lStar*0.75;
    double supBorn = lStar*1.5;
    int i(0);
    int count(0);
    int Np(0);
    double nLstar(0);
    
    for(i=0; i<allParticles.size(); i++)
    {
        Np++;
        if((allParticles[i][1]>=infBorn)&(allParticles[i][1]<supBorn))
            count++;
    }
    nLstar = count*nT/Np;
    return nLstar;
}


double beta(double l1, double l2)   // should depend on T and Knudsen not on c
{
    double K(1);
    double pi(3.1415926);
    double kB(1.38064852e-23);
    double rhoSoot(1800);
    K = 2.2*pow(((pi*kB*0.5)/(2*rhoSoot)),0.5);
    double V1 = pi*pow(l1,3)/6;
    double V2 = pi*pow(l2,3)/6;
    
    double betaCalc(1);
    betaCalc = K*pow((1/V1+1/V2),0.5)*pow((l1+l2),2);
    return betaCalc;
}

double geoAggloTotSource(vector<vector< double> > const& allParticles, double maxValL, double lp0, double a, double nT)
{
    double lStar(lp0);
    vector<double> lStarVector;
    while (lStar<=maxValL)
    {
        lStarVector.push_back(lStar);
        lStar = lStar*2;
    }
    
    double dotAt(0);
    int i(0);
    int j(0);
    for(i=0; i<lStarVector.size(); i++)
    {
        double l1 = lStarVector[i];
        for(j=0; j<lStarVector.size(); j++)
        {
            double l2 = lStarVector[j];
            dotAt = dotAt - a * 0.5*beta(l1,l2)*geometricNLstar(l1, allParticles, nT) * geometricNLstar(l2, allParticles, nT);
        }
    }
    return dotAt;
}
