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


vector<vector<double> > initAllParticles(vector<vector<double> > initVector)
{
    vector<vector<double> > allParticles;
    int i(0);
    for(i=0; i<initVector.size(); i++)
    {
        double ci = initVector[i][0];
        double li = initVector[i][1];
        double dnpi = initVector[i][2];
        int npi = floor(dnpi);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            allParticles.push_back(initPdfi);
        }
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
    int intPicked = rounding(particlesPicked);      // number of picked particles must be an integer
    int nbPicked = intPicked - (intPicked%2);      // number of picked particles must be an even number
    int maxVal(0);
    int k;
    for(k=0; k<allParticles.size(); k++)
    {
        maxVal++;                                // compteur maxVal = allParticles.size
    }
    
    maxVal = maxVal - 1;
    vector<int> randomL;
    randomL = randomListWithoutDuplicate(t, nbPicked, maxVal);
    
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
        cout << /*"c= " << allParticles[i][0] << */"  l= " << allParticles[i][1] << "   ";
    }
    cout << endl <<endl;
}


// case nucleation + growth

vector<vector<double> > initCustomized()
{
    vector<vector<double> > allParticles;
    int i(0);
    double li = 0.02;
    
    for(i=0; i< 20; i++)
    {
        double ci = 0.5;
        double fli = 10;
        int npi = rounding(50*fli);
        
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li-0.0099, li+0.0099);
            allParticles.push_back(initPdfi);
        }
        li += 0.02;
    }
    
    li = 0.42;
    i=0;
    for(i=0; i<10; i++)
    {
        double ci = 0.5;
        double fli = 100;
        int npi = rounding(50*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li-0.0099, li+0.0099);
            allParticles.push_back(initPdfi);
        }
        li += 0.02;
    }
    
    li = 0.62;
    
    for(i=0; i< 70; i++)
    {
        double ci = 0.5;
        double fli = 10;
        int npi = rounding(50*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li-0.0099, li+0.0099);
            allParticles.push_back(initPdfi);
        }
        li += 0.02;
    }
    
    return allParticles;
}




//cas pure growth
/*
vector<vector<double> > initCustomized()
{
    vector<vector<double> > allParticles;
    int i(0);
    double li = 2.1;
    for(i=0; i<19; i++)
    {
        double ci = 0.5;
        double fli = 1/pow(0.32*3.141592653,0.5)*exp(-500*pow(0.1*li-0.3, 2));
        int npi = rounding(100*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            allParticles.push_back(initPdfi);
        }
        li += 0.1;
    }
    
    li = 6.1;
    for(i=0; i<19; i++)
    {
        double ci = 0.5;
        double fli = 1;
        int npi = rounding(100*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            allParticles.push_back(initPdfi);
        }
        li += 0.1;
    }
    
    li = 10.1;
    for(i=0; i<19; i++)
    {
        double ci = 0.5;
        double fli = 1 - abs(li-11);
        int npi = rounding(100*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            allParticles.push_back(initPdfi);
        }
        li += 0.1;
    }
    
    li = 14.1;
    for(i=0; i<19; i++)
    {
        double ci = 0.5;
        double fli = pow(1 - 100*pow(0.1*li-1.5,2),0.5);
        int npi = rounding(100*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            allParticles.push_back(initPdfi);
        }
        li += 0.1;
    }
    
    return allParticles;
}
*/
