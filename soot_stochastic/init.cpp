//
//  init.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 21/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "init.hpp"
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


vector<double> initGeoMesh(double lp0, double maxValL, int nBins, double geoQ)
{
    vector<double> liVector;
    double li(lp0);
    liVector.push_back(li);
    
    int i(0);
    for(i=1; i<(nBins+1); i++)
    {
        li = lp0 + (maxValL-lp0)*pow(2.0, (double(i-nBins))/geoQ);
        liVector.push_back(li);
    }
    return liVector;
}


vector<double> initGeo2Mesh(double lp0, int nBins, double geo2Q)
{
    vector<double> liVector;
    double li(lp0);
    liVector.push_back(li);
    
    int i(0);
    for(i=1; i<nBins; i++)
    {
        li = li*geo2Q;
        liVector.push_back(li);
    }
    return liVector;
}



// case nucleation + growth
vector<vector<double> > initCustomNuclGrowth()
{
    vector<vector<double> > allParticles;
    int i(0);
    double li = 0.02;
    
    for(i=0; i< 20; i++)
    {
        double ci = 0.5;
        double fli = 10;
        int npi = rounding(5*fli);
        
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li-0.01, li+0.01);
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
        int npi = rounding(5*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li-0.01, li+0.01);
            allParticles.push_back(initPdfi);
        }
        li += 0.02;
    }
    
    li = 0.62;
    
    for(i=0; i< 70; i++)
    {
        double ci = 0.5;
        double fli = 10;
        int npi = rounding(5*fli);
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li-0.01, li+0.01);
            allParticles.push_back(initPdfi);
        }
        li += 0.02;
    }
    
    return allParticles;
}




//case pure growth
vector<vector<double> > initCustomGrowth()
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



// case pure agglomeration

vector<vector<double> > initCustomAgglo(vector<double> liVector, double maxValL)
{
    vector<vector<double> > allParticles;
    int i(0);
    double li = liVector[0];
    for(i=0; i<liVector.size(); i++)
    {
        li = liVector[i];
        double ci = 0.5;
        double fli = exp(-li);             // to determine initial pdf/ PSD
        double deltaLi(1);
        
        if(i==0)
        {
            deltaLi = li+(liVector[i+1]-li)/2.0 - (li-(liVector[i+1]-li)/2.0);
        }
        else if (li==maxValL)
        {
            deltaLi = li+(li-liVector[i-1])/2.0 - (li-(li-liVector[i-1])/2.0);
        }
        else
        {
            deltaLi = li+(liVector[i+1]-li)/2.0 -(li-(li-liVector[i-1])/2.0);
        }
        
        int npi = rounding(1000*fli);   // to determine number of particles indirectly
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        
        
        if(i==0)
        {
            int j(0);
            for(j=0; j<npi; j++)
            {
                initPdfi[1] = frand_a_b(li-(liVector[i+1]-li)/2.0, li+(liVector[i+1]-li)/2.0 );
                allParticles.push_back(initPdfi);
            }
            
        }
        else if (li == maxValL)
        {
            int j(0);
            for(j=0; j<npi; j++)
            {
                initPdfi[1] = frand_a_b(li-(li-liVector[i-1])/2.0, li+(li-liVector[i-1])/2.0);
                allParticles.push_back(initPdfi);
            }
        }
        
        else
        {
            int j(0);
            for(j=0; j<npi; j++)
            {
                initPdfi[1] = frand_a_b(li-(li-liVector[i-1])/2.0, li+(liVector[i+1]-li)/2.0);
                allParticles.push_back(initPdfi);
            }
        }
        
    }
    return allParticles;
}



// case pure agglomeration geo2

vector<vector<double> > initCustomAggloGeo2(vector<double> liVector, double maxValL)
{
    vector<vector<double> > allParticles;
    int i(0);
    double li = liVector[0];
    for(i=0; i<liVector.size(); i++)
    //for(i=0; i<3; i++)
    {
        li = liVector[i];
        double ci = 0.5;
        double fli = exp(-li);             // to determine initial pdf/ PSD
        //double fli = 1.0;
        
        int npi = rounding(1000*fli);   // to determine number of particles indirectly
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li*0.75, li*1.5);
            allParticles.push_back(initPdfi);
        }
    }
    return allParticles;
}



vector<vector<double> > initCustomAggloMass(vector<double> liVector)
{
    vector<vector<double> > allParticles;
    int i(0);
    double li = liVector[0];
    for(i=0; i<liVector.size(); i++)
    //for(i=0; i<3; i++)
    {
        li = liVector[i];
        double lavg = li*1.125;
        double ci = 0.5;
        double fli = exp(-(lavg));             // to determine initial pdf/ PSD
        double deltaLint = 0.75*li;
        //double fli = 1.0;
        
        double mli = fli *deltaLint * lavg;
        cout << "ml["<< i << "] = " << mli << endl;
        
        int npi = rounding(100000*mli);   // to determine number of particles indirectly
        
        vector<double> initPdfi;
        initPdfi.push_back(ci);
        initPdfi.push_back(li);
        
        int j(0);
        for(j=0; j<npi; j++)
        {
            initPdfi[1] = frand_a_b(li*0.75, li*1.5);
            allParticles.push_back(initPdfi);
        }
    }
    return allParticles;
}

