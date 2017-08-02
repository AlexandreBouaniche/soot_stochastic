//
//  advanceMass.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 01/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "advanceMass.hpp"
#include "aggloGeo2.hpp"
#include "aggloGeo.hpp"
#include "tools.hpp"
#include "nucleation.hpp"
#include "growth.hpp"
#include "agglomeration.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;


vector<double> allAlphaCoefMass(vector<vector< double> > const& allParticles, double a, double mT, double h, vector<vector< double> > const& lNplNvl, double timePerIt, vector<double> lVector)
{
    // source terms dotHm, dotAlm must be calculated as a function of m(l, t+growth) but before deltaT. Then for alpha coef we divide by mT(t+deltaT) -> use of mTminusOne
    
    //double dotHm = nuclSourceMass(allParticles, h);
    double dotHm =0;
    
    double alphaHm = dotHm / mT;
    
    double dotAlm0 = wmTotj(0, lNplNvl, timePerIt, lVector, a);  //mT(t-deltat)
    double alphaAlm0 = dotAlm0 /mT;
    
    double alphaLm0 = alphaAlm0 + alphaHm;
    
    vector<double> alphaVector;
    alphaVector.push_back(alphaLm0);
    
    //cout << "alphaLm[0] =  " << alphaLm0 << endl;
    //cout << "dotAlm[0] = "<< dotAlm0 << endl;
    
    int i(0);
    for(i=1; i<lNplNvl.size(); i++)     // begins at i=1 because we already calculated the first term alphaLm0 (corresponds to lp0 and i = 0)
    {
        double dotALmi = wmTotj(i, lNplNvl, timePerIt, lVector, a); //nT(t-deltat)
        double alphaLmi = dotALmi/mT;
        alphaVector.push_back(alphaLmi);
        //cout << "alphaLm[" << i << "] =  " << alphaLmi << endl;
        //cout << "dotAlm["<<i << "] = "<< dotALmi << endl;
    }
    i = 0;
    return alphaVector;
}



void advancePdfMass(vector<double>const& alphaVector, vector<vector< double> >& allParticles, vector<vector< double> > & lNplNvl, double h, double mT, double a,double it, double mTminusOne, double timePerIt, vector<double> lVector)
{
    //count of Np
    int Np(0);
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        Np++;
    }
    cout << "Nptot = " << Np << endl;
    
    // count of li bins (=alphaVector.size() = lNplNvl.size()
    i=0;
    int countLiBins(0);
    for(i=0; i<alphaVector.size(); i++)
    {
        countLiBins++;
    }
    
    
    // calculation of alphaH and alphaAt which are not in alphaVector. All the alphaL* are in alphaVector
    double dotHm(0);
    double dotAtm(0);
    
    //dotHm = nuclSourceMass(allParticles, h);
    dotHm = 0.0;
    
    
    // calculation of sum of the wm to get dotAtm. calculated with lNplNvl. info t-deltaT
    i=0;
    for (i=0; i<lNplNvl.size(); i++)
    {
        dotAtm += wmTotj(i, lNplNvl, timePerIt, lVector, a);
    }
    
    double alphaHm = dotHm/mT;
    double alphaAtm = dotAtm/mT;
    
    double alphaHmplusAtm = alphaHm + alphaAtm;
    
    
    // calculation of np(l*, t+dt)   (called nplDt here) and storing of the corresponding  Delta np = np(l*,t+dt) - np(l*,t) in a new vector: deltaNpInt. New vector because before we were storing integers in a vector of doubles (lNplNvl)
    
    
    vector<int> deltaNpInt;
    
    i = 0;
    int deltaNplSum(0);
    for(i=0; i<lNplNvl.size(); i++)
    {
        double nplDt(0);
        nplDt = lNplNvl[i][1] * (1 - alphaHmplusAtm) + alphaVector[i] * Np;  //calculation of np(l*,t+dt)
        
        int roundednp2 = rounding(nplDt);                 // rounding function
        int roundednp1 = rounding(lNplNvl[i][1]);         // rounding function
        int deltaNpl = roundednp2 - roundednp1;
        
        //cout << "deltaNpl["<< lNplNvl[i][0] << "] = " << deltaNpl << endl;
        
        deltaNpInt.push_back(deltaNpl);                    //filling vector deltaNpInt with DeltaNp(l*) integer
        
        deltaNplSum = deltaNplSum + deltaNpl;              // In theory should be zero. with rounding errors not zero
        
        
    }
    cout << endl << "deltaNplSum = " << deltaNplSum << endl;
    
    
    // there is a rounding error. deltaNplSum should be equal to 0. we will use its value to correct the error during particles reallocation
    
    
    
    vector<int> ranksToReallocate;
    vector<double> valuesToRealoc;
    int countRanksToRealoc(0);
    
    i=0;
    double li(0);
    
    double infborn(0);
    double supborn(0);
    
    for(i=0; i<lNplNvl.size(); i++)
    {
        li = lNplNvl[i][0];
        infborn = 0.75*li;
        supborn = 1.5*li;
        
        int deltaNpli = deltaNpInt[i];
        if(deltaNpli<0)             // here only for li with deltaNpli < 0 (particles "consumed")
        {
            vector<int> allPartLiRanks;        //vector of int: rank in allParticles of the particles of size li
            int countAllPartLi(0);
            int j(0);
            for(j=0; j<allParticles.size(); j++)
            {
                
                if(allParticles[j][1]>=infborn & allParticles[j][1]<supborn)
                {
                    allPartLiRanks.push_back(j);
                    countAllPartLi++;
                }
            }
            
            int minPart = -deltaNpli;
            
            if(-deltaNpli > countAllPartLi)
            {
                cout << "more particles consumed than present for li = " << li << endl;
                minPart = countAllPartLi;
            }
            
            
            if(countAllPartLi>1)
            {
                vector<int> randomL = randomListWithoutDuplicate(it, minPart, (countAllPartLi-1));
                //countAllPartLi = np(li) (integer) = allPartLiRanks.size. We must take this value -1 because the ranks of the vector allPartLiRanks go from 0 to (countAllPartLi - 1) if we don't put -1 RandomL can pick a rank of allPartLiRanks that doesn't exist
                
                //rankLiAndNbToPick[i][1] particles of size li are picked within the np(li) particles of size li.Their rank IN allPartLiRanks is stored in the vector randomL
                
                j=0;
                for(j=0; j<randomL.size(); j++)            // then we have to "traduce" a rank in allPartLiRanks in a rank in allParticles
                {
                    int rankToRe(0);
                    int picked(0);
                    picked = randomL[j];                     // picked is the rank in allPartLiRanks
                    rankToRe = allPartLiRanks[picked];       // "traducing" rank in allPartLiRanks in rankToRe which is a rank in allParticles as allPartLiRanks stores the ranks in allParticles for size li
                    ranksToReallocate.push_back(rankToRe);
                    countRanksToRealoc++;
                    //cout << "rank["<< i << "] to reallocate = " << rankToRe << "   ;";
                }
                //cout << "ranksToReallocate vector built" << endl;
            }
            else if (countAllPartLi == 1)
            {
                ranksToReallocate.push_back(allPartLiRanks[0]);
            }
            
        }
        
        
        if(deltaNpli>0)
        {
            
            int j=0;
            for(j=0; j< deltaNpli; j++)
            {
                valuesToRealoc.push_back(li);
            }
            
            //cout << "valuesToReallocate vector built" << endl;
        }
        
        
    }
    
    
    //test vector valuesToRealoc  -> ok
    
    i=0;
    int countValuesToRealoc(0);
    for(i=0; i<valuesToRealoc.size(); i++)
    {
        countValuesToRealoc++;
    }
    
    cout << "number of values to set = " << countValuesToRealoc;
    cout << endl;
    
    
    
    // Now reallocating!
    i=0;
    int min = countValuesToRealoc;              // due to rounding error
    if(countRanksToRealoc < countValuesToRealoc)
        min = countRanksToRealoc;
    
    for(i=0; i<min; i++)
    {
        double infborni(0);
        double supborni(0);
        li = valuesToRealoc[i];
        infborni = 0.75*li;
        supborni = 1.5*li;
        
        
        int rankParticle(0);
        rankParticle = ranksToReallocate[i];
        double newVal(0);
        newVal = frand_a_b(infborni, supborni); // a random value in the interval of valuesToRealoc[i] is chosen to reset the particle
        allParticles[rankParticle][1] = newVal;
        
    }
    cout << "number of ranks to reset = " << countRanksToRealoc << endl << endl;  //test
    
}


double aggloTotMassSource(vector<vector<double> > const& lNplNvl, double timePerIt, vector <double> const& lVector, double a)
{
    int i(0);
    double dotAtm(0);
    for(i=0; i<lVector.size(); i++)
    {
        dotAtm += wmTotj(i, lNplNvl, timePerIt, lVector, a);
    }
    return dotAtm;
}




double totalMassBins(vector<vector<double> > const& lNpNv)
{
    int i(0);
    double totMass(0);
    for (i=0; i<lNpNv.size(); i++)
    {
        totMass += lNpNv[i][3];    
    }
    return totMass;
}



vector<vector<double> > lNplNvMass(vector<vector<double> > allParticles, vector<double> liVector, double mT)   // 0:size;  1:Npl;  2:Nvl = ml / (size*1.125);  3:ml = Npl / Nptot * mTot
{
    vector<double> lVector = liVector;
    vector<vector<double> > lAndNpl;
    
    int j(0);
    int Np(0);
    for(j=0;j<allParticles.size();j++)
    {
        Np++;                    // count of tot Np for calculating nv from np
    }
    
    int i(0);
    for(i=0; i<lVector.size(); i++)
    {
        int npL(0);
        double li(0);
        double nvL(0);
        double mvL(0);
        li = lVector[i];
        double infborn = li*0.75;
        double supborn = li*1.5;
        
        // count of np(li)
        j=0;
        for(j=0;j<allParticles.size();j++)
        {
            if((allParticles[j][1]>=infborn)&(allParticles[j][1]<supborn))
            {
                npL++;
            }
        }
        
        mvL = double(npL)/double(Np)*mT;      // mass stochastic particles different from number density stochastic particles. One stochastic particle represents a mass mT/Np. Before we had one stochastic particle represents nT/Np particles
        nvL = mvL/(li*1.125);
        lAndNpl.push_back(vector<double>(4,0));
        lAndNpl[i][0] = li;
        lAndNpl[i][1] = npL;
        lAndNpl[i][2] = nvL;
        lAndNpl[i][3] = mvL;
        
    }
    return lAndNpl;
}



double nTfromMass(vector<vector<double> > lNpNvMass)
{
    double nT(0);
    int i(0);
    for (i=0; i<lNpNvMass.size(); i++)
    {
        nT += lNpNvMass[i][2];
    }
    return nT;
}
