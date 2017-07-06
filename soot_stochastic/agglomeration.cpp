//
//  agglomeration.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 26/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "agglomeration.hpp"
#include "tools.hpp"
#include "nucleation.hpp"
#include "growth.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

vector<double> liVector(double lp0, double deltaL, double maxValL)
{
    int i(0);
    double steps = (maxValL - lp0) / deltaL;
    int stepsInt = floor(steps) + 1;
    vector<double> lVector;
    double li(lp0);
    for(i=0; i<stepsInt; i++)
    {
        lVector.push_back(li);
        li = li+deltaL;
    }
    return lVector;
}

int npLstar(double lStar, vector<vector<double> > allParticles, double deltaL)
{
    double infBorn = lStar-deltaL/2;
    double supBorn = lStar+deltaL/2;
    int i(0);
    int count(0);
    for(i=0; i<allParticles.size(); i++)
    {
        if((allParticles[i][1]>=infBorn)&(allParticles[i][1]<supBorn))
            count++;
    }
    return count;
}


// returns a vector<vector<doube> > with 3 columns, 0:li  1:np(li)  2:nv(li)
vector<vector<double> > liNpliNvli(vector<vector<double> > allParticles, vector<double> liVector, double deltaL, double nT)
{
    vector<double> lVector = liVector;
    vector<vector<double> > lAndNpl;
    int i(0);
    for(i=0; i<lVector.size(); i++)
    {
        int npL(0);
        double l(0);
        double nvL(0);
        l = lVector[i];
        npL = npLstar(l, allParticles, deltaL);
        nvL = nvLstar(l, allParticles, nT, deltaL);
        lAndNpl.push_back(vector<double>(3,0));
        lAndNpl[i][0] = l;
        lAndNpl[i][1] = npL;
        lAndNpl[i][2] = nvL;
    }
    return lAndNpl;
}

double nvLstar(double lStar, vector<vector<double> > allParticles, double nT, double deltaL)
{
    double infBorn = lStar-deltaL/2;
    double supBorn = lStar+deltaL/2;
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
    /*
     For now we consider that l1 and l2 are homogeneous to a mass/volume. Let's say lp0 is the volume of a sphere of diameter realLp0 = 1 nm
     
     */
    double K(1);
    double pi(3.1415926);
    double kB(1.38064852e-23);
    double rhoSoot(1800);
    K = 2.2*pow(((pi*kB*0.5*2000)/(2*rhoSoot)),0.5);  // to do: substitute 0.5 by avg c for particles of size l1 and l2. 2000 is an approximation of T if c = 1
    
    /*
    double V1 = pi*pow(l1,3)/6;          // if we consider volume corresponding to diameter. for now we will consider that l is representative of a "size", eg a volume of mass and lc = l* - li
    double V2 = pi*pow(l2,3)/6;
    */
    
    double V1 = l1;
    double V2 = l2;
    double realL1 = pow((6*l1/pi),0.3333);   //realL1 is the diameter the particle would have, to have a volume of "l1"
    double realL2 = pow((6*l2/pi),0.3333);
    
    double betaCalc(1);
    betaCalc = K*pow((1/V1+1/V2),0.5)*pow((realL1+realL2),2);
    //return betaCalc;
    double betaTest =1.0;
    return betaTest;
}



double aggloTotSource(vector<vector< double> > const& allParticles, vector<vector< double> > const& lNplNvl, double a)
{
    double l1(0);
    double l2(0);
    double nvl1(0);
    double nvl2(0);
    double dotAt(0);
    int i(0);
    int j(0);
    for(i=0; i<lNplNvl.size(); i++)
    {
        l1 = lNplNvl[i][0];
        nvl1 = lNplNvl[i][2];
        for(j=0; j<lNplNvl.size(); j++)
        {
            l2 = lNplNvl[j][0];
            nvl2 = lNplNvl[j][2];
            dotAt = dotAt - a * 0.5*beta(l1,l2)*nvl1 * nvl2;
        }
    }
    return dotAt;
}


double dotAlStar(double lStar, vector<vector< double> > const& allParticles, vector<vector< double> > const& lNplNvl, double a, double deltaL, double nT)
{
    double AlStarNeg(0);
    int i(0);
    double ls(lStar);
    double li(0);
    int rankLstar(0);
    
    for(i=0; i<lNplNvl.size(); i++)          // looking for closest value (to lStar) of li of the lNplNvl vector.
    {
        li = lNplNvl[i][0];
        if((ls>=(li-deltaL/2))&(ls<(li+deltaL/2)))
        {
            ls = li;                // taking this closest value for source calculation
            rankLstar = i;          // corresponding rank in lNplNvl
        }
    }
    
    ls = lNplNvl[rankLstar][0];
    
    i=0;
    for(i=0; i<lNplNvl.size(); i++)
    {
        li = lNplNvl[i][0];
        AlStarNeg = AlStarNeg - a * beta(ls,li) * lNplNvl[rankLstar][2] * lNplNvl[i][2];  //negative term of Al*
    }
    
    double AlStarPos(0);
    double lc = ls;
    i=0;
    for(i=0; i<(rankLstar); i++)
    {
        li = lNplNvl[i][0];
        lc = ls-li;
        AlStarPos = AlStarPos + a * 0.5 * beta(lc,li) * nvLstar(lc, allParticles, nT, deltaL) * lNplNvl[i][2];
        //positive term of Al*
        
    }
    
    double AlStarTot(0);
    AlStarTot = AlStarPos + AlStarNeg;   // Al* total
    return AlStarTot;
}


vector<double> allAlphaCoef(vector<vector< double> > const& allParticles, double lp0, double a, double nT, double nTtminusOne, double h, double deltaL, vector<vector< double> > const& lNplNvl)
{
    double dotH = nuclSource(allParticles, h);
    double alphaH = dotH / nT;
    double dotAl0 = dotAlStar(lp0, allParticles, lNplNvl, a, deltaL, nTtminusOne);  //nT(t-deltat)
    double alphaAl0 = dotAl0 /nT;
    
    double alphaL0 = alphaAl0 + alphaH;  // added oxidation of nascent particles
    
    vector<double> alphaVector;
    alphaVector.push_back(alphaL0);
    
    int i(0);
    for(i=1; i<lNplNvl.size(); i++)     // begins at i=1 because we already calculated the first term alphaL0 (corresponds to lp0 and i = 0)
    {
        double l1 = lNplNvl[i][0];
        double dotALi = dotAlStar(l1, allParticles, lNplNvl, a, deltaL, nTtminusOne); //nT(t-deltat)
        double alphaLi = dotALi/nT;
        alphaVector.push_back(alphaLi);
    }
    return alphaVector;
}



void advancePdf(vector<double>const& alphaVector, vector<vector< double> >& allParticles, vector<vector< double> > & lNplNvl, double h, double nT, double a, double deltaL, double t)
{
    //count of Np
    int Np(0);
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        Np++;
    }
    
    // count of li bins (=alphaVector.size() = lNplNvl.size()
    i=0;
    int countLiBins(0);
    for(i=0; i<alphaVector.size(); i++)
    {
        countLiBins++;
    }
    
    
    // calculation of alphaH and alphaAt which are not in alphaVector. All the alphaL* are in alphaVector
    double dotH(0);
    double dotAt(0);
    
    dotH = nuclSource(allParticles, h);
    dotAt = aggloTotSource(allParticles, lNplNvl, a);
    
    double alphaH = dotH/nT;
    double alphaAt = dotAt/nT;
    
    double alphaHplusAt = alphaH + alphaAt;   // added oxidation of nascent particles
    
    
    // calculation of np(l*, t+dt)   (called nplDt here) and storing of the corresponding  Delta np = np(l*,t+dt) - np(l*,t) in a new vector: deltaNpInt. New vector because before we were storing integers in a vector of doubles (lNplNvl)
    
    
    vector<int> deltaNpInt;
    
    i = 0;
    int deltaNplSum(0);
    for(i=0; i<lNplNvl.size(); i++)
    {
        double nplDt(0);
        nplDt = lNplNvl[i][1] * (1 - alphaHplusAt) + alphaVector[i] * Np;  //calculation of np(l*,t+dt)
        
        int roundednp2 = rounding(nplDt);                 // rounding function
        int roundednp1 = rounding(lNplNvl[i][1]);         // rounding function
        int deltaNpl = roundednp2 - roundednp1;
        
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
    for(i=0; i<lNplNvl.size(); i++)
    {
        li = lNplNvl[i][0];
        int deltaNpli = deltaNpInt[i];
        if(deltaNpli<0)             // here only for li with deltaNpli < 0 (particles "consumed")
        {
            vector<int> allPartLiRanks;        //vector of int: rank in allParticles of the particles of size li
            int countAllPartLi(0);
            int j(0);
            for(j=0; j<allParticles.size(); j++)
            {
                if(allParticles[j][1]>=(li-deltaL/2) & allParticles[j][1]<(li+deltaL/2))
                {
                    allPartLiRanks.push_back(j);
                    countAllPartLi++;
                }
            }
            vector<int> randomL = randomListWithoutDuplicate(t, (-deltaNpli), (countAllPartLi-1));
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
        }
        if(deltaNpli>0)
        {
            
            int j=0;
            for(j=0; j< deltaNpli; j++)
            {
                valuesToRealoc.push_back(li);
            }
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
        int rankParticle(0);
        rankParticle = ranksToReallocate[i];
        allParticles[rankParticle][1] = valuesToRealoc[i];
        
    }
    cout << "number of ranks to reset = " << countRanksToRealoc << endl << endl;  //test
    
}

