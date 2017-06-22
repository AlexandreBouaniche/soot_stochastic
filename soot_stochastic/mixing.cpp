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


vector<double> initParticles(int Np0, int Np1, double c0, double c1)
{
    vector<double> allParticles;
    int i(0);
    for(i=0; i< Np0; i++)
    {
        allParticles.push_back(c0);
    }
    int j(0);
    for(j=0; j< Np1; j++)
    {
        allParticles.push_back(c1);
    }
    return allParticles;
}


void mix(std::vector<double>& allParticles, double deltaT, double tau, double t)
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
    int maxVal = allParticles.size();
    maxVal = maxVal - 1;
    vector<int> randomL;
    randomL = randomList(t,nbPicked, maxVal);
    
    int i;
    for(i=0; i<randomL.size(); i=i+2)
    {
        int rank1 = randomL[i];
        int rank2 = randomL[i+1];
        double newVal = (allParticles[rank1] + allParticles[rank2]) /2;
        allParticles[rank1] = newVal;
        allParticles[rank2] = newVal;
    }
    
    /*
    for(i=0; i<mixEvIt; i++)
    {
        double newParticle;
        newParticle = (allParticles[i] + allParticles[Nptot-i-1]) / 2;
        allParticles[i] = newParticle;
        allParticles[Nptot-i-1] = newParticle;
        cout << "mix event number " << i << " new value: " << newParticle << endl;
    }
     */
    
    
}


void printParticles(std::vector<double> const& allParticles, double t)
{
    cout << "t = " << t << endl << endl;
    int i;
    for(i=0; i<allParticles.size();i++)
    {
        cout << allParticles[i] << "   ";
    }
    cout << endl <<endl;
}


