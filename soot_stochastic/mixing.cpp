//
//  mixing.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "mixing.hpp"
#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;


vector<vector<double> > initParticles(int Np0, int Np1, double c0, double c1, double l0, double l1)
{
    vector<vector<double> > allParticles;
    int i(0);
    vector<double> cL0;
    cL0.push_back(c0);
    cL0.push_back(l0);
    
    vector<double> cL1;
    cL1.push_back(c1);
    cL1.push_back(l1);
    
    for(i=0; i< Np0; i++)
    {
        allParticles.push_back(cL0);
    }
    int j(0);
    for(j=0; j< Np1; j++)
    {
        allParticles.push_back(cL1);
    }
    return allParticles;
}


void mix(vector<vector<double> >& allParticles, double deltaT, double tau, double t)
{
    int Nptot(0);
    int j(0);
    for (j=0; j<allParticles.size();j++)  // compteur Nptot
    {
        Nptot++;
    }
    
    // construction of the vector randomL of random picked particles for mixing
    double particlesPicked = deltaT/tau*Nptot;   // number of picked particles depends on tau
    int intPicked = floor(particlesPicked);      // number of picked particles must be an integer
    int nbPicked = intPicked - (intPicked%2);      // number of picked particles must be an even number
    int maxVal(0);
    int k;
    for(k=0; k<allParticles.size(); k++)
    {
        maxVal++;                                // compteur maxVal = allParticles.size
    }
    
    maxVal = maxVal - 1;
    vector<int> randomL;
    randomL = randomList(t,nbPicked, maxVal);
    
    int i;
    for(i=0; i<randomL.size(); i=i+2)
    {
        int rank1 = randomL[i];
        int rank2 = randomL[i+1];
        double newVal = (allParticles[rank1][0] + allParticles[rank2][0]) /2;
        allParticles[rank1][0] = newVal;
        allParticles[rank2][0] = newVal;
    }
}


void printParticles(vector<vector<double> > const& allParticles, double t)
{
    cout << "t = " << t << endl << endl;
    int i;
    for(i=0; i<allParticles.size();i++)
    {
        cout << "c= " << allParticles[i][0] << "  l= " << allParticles[i][1] << "   ";
    }
    cout << endl <<endl;
}


