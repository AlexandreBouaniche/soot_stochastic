//
//  aggloGeo.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 24/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "aggloGeo.hpp"
#include "tools.hpp"
#include "nucleation.hpp"
#include "growth.hpp"
#include "geometric.hpp"
#include "agglomeration.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;




double nvLstarGeo(double lStar, vector<vector<double> > lNplNvl, double maxvalL)
{
    
    // infborn and supborn calculation for each li
    double infborn(0);
    double supborn(0);
    double nLstar(0);
    
    int i(0);
    for(i=0; i<lNplNvl.size(); i++)
    {
        double li = lNplNvl[i][0];
        
        if(i==0)
        {
            infborn = li-(lNplNvl[i+1][0]-li)/2.0;
            supborn = li+(lNplNvl[i+1][0]-li)/2.0;
        }
        else if (li == maxvalL)
        {
            infborn = li-(li-lNplNvl[i-1][0])/2.0;
            supborn = li+(li-lNplNvl[i-1][0])/2.0;
        }
        else
        {
            infborn = li-(li-lNplNvl[i-1][0])/2.0;
            supborn = li+(lNplNvl[i+1][0]-li)/2.0;
        }
        
        
        
        if(lStar>=infborn & lStar < supborn)
        {
            nLstar = lNplNvl[i][2];
        }
    }
    return nLstar;
}



double dotAlStarGeo(double lStar, vector<vector< double> > const& allParticles, vector<vector< double> > const& lNplNvl, double a, double deltaL, double nT, double timePerIt, double maxValL)
{
    double AlStarNeg(0);
    int i(0);
    double ls(lStar);
    double li(0);
    int rankLstar(0);
    
    double infborn(0);
    double supborn(0);
    
    double infbornLstar(0);
    double supBornLstar(0);
    
    for(i=0; i<lNplNvl.size(); i++)          // looking for closest value (to lStar) of li of the lNplNvl vector.
    {
        li = lNplNvl[i][0];
        
        if(i==0)
        {
            infborn = li-(lNplNvl[i+1][0]-li)/2.0;
            supborn = li+(lNplNvl[i+1][0]-li)/2.0;
        }
        else if (li == maxValL)
        {
            infborn = li-(li-lNplNvl[i-1][0])/2.0;
            supborn = li+(li-lNplNvl[i-1][0])/2.0;
        }
        else
        {
            infborn = li-(li-lNplNvl[i-1][0])/2.0;
            supborn = li+(lNplNvl[i+1][0]-li)/2.0;
        }
        
        //cout << "li = " << li << endl;
        //cout << infborn << "   " << supborn << endl << endl;
        
        
        if(lStar>=infborn & lStar < supborn)
        {
            ls = li;                // taking this closest value for source calculation
            rankLstar = i;          // corresponding rank in lNplNvl
            supBornLstar = supborn;
            infbornLstar = infborn;
        }

    }
    
    ls = lNplNvl[rankLstar][0];
    
    i=0;
    for(i=0; i<lNplNvl.size(); i++)
    {
        li = lNplNvl[i][0];
        AlStarNeg = AlStarNeg - a * beta(ls,li, timePerIt) * lNplNvl[rankLstar][2] * lNplNvl[i][2];  //negative term of Al*
    }
    
    double AlStarPos(0);
    i=0;
    for(i=0; i<=(rankLstar); i++)
    {
        li = lNplNvl[i][0];
        
        int j(0);
        for(j=0; j<=rankLstar; j++)
        {
            double lc = lNplNvl[j][0];
            if ((li+lc)>=infbornLstar & (li+lc)<supBornLstar)      //positive term of Al*. Consider distribution from infborn to supborn?
                // next development: sum over all li and lc with stoichiometric coef
            {
                 AlStarPos = AlStarPos + a * 0.5 * beta(lc,li, timePerIt) * nvLstarGeo(lc, lNplNvl, maxValL) * lNplNvl[i][2];
            }
        }
        //positive term of Al*
    }
    double AlStarTot(0);
    AlStarTot = AlStarPos + AlStarNeg;   // Al* total
    
    return AlStarTot;
}


vector<double> allAlphaCoefGeo(vector<vector< double> > const& allParticles, double lp0, double a, double nT, double nTtminusOne, double h, double deltaL, vector<vector< double> > const& lNplNvl, double t, double timePerIt, double maxValL)
{
    // source terms dotH, dotAl, must be calculated in function of n(l,t+growth) but before deltaT. Then for alpha coef divided by nT(t+deltaT) -> use of nTminusOne
    
    double dotH = nuclSource(allParticles, h);
    //double dotH = nuclSourceCustomized(t, deltaL);
    double alphaH = dotH / nT;
    
    double dotAl0 = dotAlStarGeo(lp0, allParticles, lNplNvl, a, deltaL, nTtminusOne, timePerIt, maxValL);  //nT(t-deltat)
    double alphaAl0 = dotAl0 /nT;
    
    double alphaL0 = alphaAl0 + alphaH;
    
    vector<double> alphaVector;
    alphaVector.push_back(alphaL0);
    
    int i(0);
    for(i=1; i<lNplNvl.size(); i++)     // begins at i=1 because we already calculated the first term alphaL0 (corresponds to lp0 and i = 0)
    {
        double l1 = lNplNvl[i][0];
        double dotALi = dotAlStarGeo(l1, allParticles, lNplNvl, a, deltaL, nTtminusOne, timePerIt, maxValL); //nT(t-deltat)
        double alphaLi = dotALi/nT;
        alphaVector.push_back(alphaLi);
        //cout << "alphaL[" << i << "] =  " << alphaLi << endl;
        //cout << i << "   "<< dotALi << endl;
    }
    i = 0;
    return alphaVector;
}



void advancePdfGeo(vector<double>const& alphaVector, vector<vector< double> >& allParticles, vector<vector< double> > & lNplNvl, double h, double nT, double a, double deltaL, double it, double maxValL, double lp0, double nTtminusOne, double timePerIt)
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
    double dotH(0);
    double dotAt(0);
    
    dotH = nuclSource(allParticles, h);
    //dotH = nuclSourceCustomized(t, deltaL);
    dotAt = aggloTotSource(allParticles, lNplNvl, a, timePerIt); // calculated with lNplNvl. info t-deltaT
    
    double alphaH = dotH/nT;
    double alphaAt = dotAt/nT;
    
    double alphaHplusAt = alphaH + alphaAt;
    
    
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
        
        if(i==0)
        {
            infborn = li-(lNplNvl[i+1][0]-li)/2.0;
            supborn = li+(lNplNvl[i+1][0]-li)/2.0;
        }
        else if (li == maxValL)
        {
            infborn = li-(li-lNplNvl[i-1][0])/2.0;
            supborn = li+(li-lNplNvl[i-1][0])/2.0;
        }
        else
        {
            infborn = li-(li-lNplNvl[i-1][0])/2.0;
            supborn = li+(lNplNvl[i+1][0]-li)/2.0;
        }
        
        
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
        
        // additional lines to determine infborn and supborn of valuesToRealoc[i] for the random reallocation within the right intervals - beginning
        double infborni(0);
        double supborni(0);
        int rankInLvector(0);
        li = valuesToRealoc[i];
        
        int j(0);
        while(li<lNplNvl[j][0])
        {
            j++;
        }
        rankInLvector = j;
        
        
        if(rankInLvector==0)
        {
            infborni = lNplNvl[j][0]-(lNplNvl[j+1][0]-li)/2.0;
            supborni = lNplNvl[j][0]+(lNplNvl[j+1][0]-li)/2.0;
        }
        else if (li == maxValL)
        {
            infborni = li-(li-lNplNvl[j-1][0])/2.0;
            supborni = li+(li-lNplNvl[j-1][0])/2.0;
        }
        else
        {
            infborni = li-(li-lNplNvl[j-1][0])/2.0;
            supborni = li+(lNplNvl[j+1][0]-li)/2.0;
        }
        
         // additional lines to determine infborn and supborn for the random reallocation within the right intervals - end
        
        int rankParticle(0);
        rankParticle = ranksToReallocate[i];
        double newVal(0);
        newVal = frand_a_b(valuesToRealoc[i]-infborni, valuesToRealoc[i]+supborni); // a random value in the interval of valuesToRealoc[i] is chosen to reset the particle
        allParticles[rankParticle][1] = newVal;
        
    }
    cout << "number of ranks to reset = " << countRanksToRealoc << endl << endl;  //test
    
}
