//
//  advanceNdf.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 03/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "advanceNdf.hpp"
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

vector<vector<double> > ndf(vector<vector<double> > allParticles, vector<double> liVector, double nT)   // 0: size li  != lavg;    [mass (or volume)]
          // 1: Npl: Number of stochastic particles in the bin;  []
          // 2: Nvl = Ndl * deltaLintervali  [part/volume]
          // could have but not necessary Ndl [part/volume/xcoordinate in size space]
          // could have but not necessary 4:ml = Nvl * lavg = Ndl * deltaLintervali * lavg [mass/volume]
          // [part/volume]   nT = sum of the Nvli = sum of the Ndi * deltaLinti   -> relation used to get Nd from nT and Npl. nT evolved with global source terms dotAt, dotH ...

{
    vector<double> lVector = liVector;
    vector<vector<double> > ndft;
    
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
        double nvl(0);
        
        li = lVector[i];
        double infborn = li*0.75;
        double supborn = li*1.5;
        
        // count of np(li)
        j=0;
        for(j=0;j<allParticles.size();j++)
        {
            if((allParticles[j][1]>=infborn)&&(allParticles[j][1]<supborn))
            {
                npL++;
            }
        }
        
        nvl = double(npL)/double(Np)*nT;    // [part/volume]
        
        ndft.push_back(vector<double>(3,0));
        ndft[i][0] = li;
        ndft[i][1] = npL;
        ndft[i][2] = nvl;
        
    }
    return ndft;
}




vector<double> allAlphaCoefNdf(vector<vector< double> > const& allParticles, double a, double nT, double h, vector<vector< double> > const& ndft, double timePerIt, vector<double> lVector)
{
    // The stochastic particles repartition is weighted by nvl not ml but the source terms are calculated in mass then "traduced" in particles.
    // in advanceMass the stochastic particles repartition the ndf was weighted by mass. -> pb for small particles
    // source terms dotH, dotAl must be calculated as a function of n(l, t+growth) but before deltaT. Then for alpha coef we divide by nT(t+deltaT) -> use of nTminusOne
    
    //double dotH = nuclSourceMass(allParticles, h);
    
    double wmTotAll(0);
    double dotH =h;
    
    double alphaH = dotH / nT;     //  [ /time ]
    
    double lp0 = lVector[0];
    double lavg0 = lp0*1.125;
    
    double dotAl0 = wmTotj(0, ndft, timePerIt, lVector, a) / lavg0;  // [part/volume/time]
    // based on nT(t-deltat)
    
    wmTotAll += wmTotj(0, ndft, timePerIt, lVector, a);
    
    
    double alphaAl0 = dotAl0 /nT;   //  [ /time ]
    
    double alphaL0 = alphaAl0 + alphaH;
    
    vector<double> alphaVector;
    alphaVector.push_back(alphaL0);
    
    
    int i(0);
    for(i=1; i<ndft.size(); i++)     // begins at i=1 because we already calculated the first term alphaL0 (corresponds to lp0 and i = 0)
    {
        double li = lVector[i];
        double lavgi = li * 1.125;
        
        double dotALi = wmTotj(i, ndft, timePerIt, lVector, a) / lavgi; //nT(t-deltat)
        
        wmTotAll += wmTotj(i, ndft, timePerIt, lVector, a);
        
        double alphaLi = dotALi/nT;
        alphaVector.push_back(alphaLi);
    }
    
    cout << "wmTotAll = " << wmTotAll << endl;
    return alphaVector;
}




void advanceNdf(vector<double>const& alphaVector, vector<vector< double> >& allParticles, vector<vector< double> > & ndft, double h, double nT, double a,double it, double timePerIt, vector<double> lVector)
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
    
    //dotH = nuclSourceMass(allParticles, h);
    dotH = h;
    
    
    // calculation of sum of the wm to get dotAt. calculated with ndft. info t-deltaT
    i=0;
    for (i=0; i<ndft.size(); i++)
    {
        double li = lVector[i];
        double lavgi = li*1.125;
        dotAt += wmTotj(i, ndft, timePerIt, lVector, a) / lavgi;  // [part/volume/time]
    }
    cout << "dotAt = " << dotAt << endl;
    double alphaH = dotH/nT;
    double alphaAt = dotAt/nT;
    
    double alphaHplusAt = alphaH + alphaAt;
    
    
    // calculation of np(l*, t+dt)   (called nplDt here) and storing of the corresponding  Delta np = np(l*,t+dt) - np(l*,t) in a new vector: deltaNpInt. New vector because before we were storing integers in a vector of doubles (lNplNvl)
    
    
    vector<int> deltaNpInt;
    
    i = 0;
    int deltaNplSum(0);
    for(i=0; i<ndft.size(); i++)
    {
        double nplDt(0);
        nplDt = ndft[i][1] * (1 - alphaHplusAt) + alphaVector[i] * Np;  //calculation of np(l*,t+dt)
        
        int roundednp2 = rounding(nplDt);                 // rounding function
        int roundednp1 = rounding(ndft[i][1]);         // rounding function
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
    
    for(i=0; i<ndft.size(); i++)
    {
        li = ndft[i][0];
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
                
                if(allParticles[j][1]>=infborn && allParticles[j][1]<supborn)
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



double totalMassNdf(vector<vector<double> > ndft)
{
    double sootVoltot(0);
    double rhoSoot = 1800;     // 1800 [kg/m3soot]   if volume in nm3  then multiply by 1e-(9*3)
    double mtot(0);
    int i(0);
    for(i=0; i< ndft.size(); i++)
    {
        double li = ndft[i][0];
        double lavgi = li*1.125;
        double volumeli = ndft[i][2] * lavgi * 1e-27;     //  [part / m3gas] * [m3soot/part].  lavgi is in [nm3soot] so multiply by 1e-(9*3) to get in [m3soot]
        sootVoltot += volumeli;                   //  [m3soot / m3gas]
    }
    mtot = rhoSoot*sootVoltot;              //  [kg/m3soot] * [m3soot/m3gas] = [kgsoot/m3gas] 
    
    return mtot;
}





vector<double> alphaCustomNuclAllSizes(vector<double> const& dotHvector, vector<vector< double> > const& allParticles, double a, double nT, double h, vector<vector< double> > const& ndft, double timePerIt, vector<double> lVector)
{
    double wmTotAll(0);
    double dotH0 =dotHvector[0];
    
    double alphaH0 = dotH0 / nT;     //  [ /time ]
    
    double lp0 = lVector[0];
    double lavg0 = lp0*1.125;
    
    double dotAl0 = wmTotj(0, ndft, timePerIt, lVector, a) / lavg0;  // [part/volume/time]
    // based on nT(t-deltat)
    
    wmTotAll += wmTotj(0, ndft, timePerIt, lVector, a);
    
    
    double alphaAl0 = dotAl0 /nT;   //  [ /time ]
    
    double alphaL0 = alphaAl0 + alphaH0;
    
    vector<double> alphaVector;
    alphaVector.push_back(alphaL0);
    
    
    int i(0);
    for(i=1; i<ndft.size(); i++)     // begins at i=1 because we already calculated the first term alphaL0 (corresponds to lp0 and i = 0)
    {
        double li = lVector[i];
        double lavgi = li * 1.125;
        
        double dotALi = wmTotj(i, ndft, timePerIt, lVector, a) / lavgi; //nT(t-deltat)
        
        wmTotAll += wmTotj(i, ndft, timePerIt, lVector, a); // just in order to check wmTotAll. Not used for alphaCoef.
        
        double alphaLi = dotALi/nT + dotHvector[i]/nT;  // in this custom case, added nucleation at all sizes.
        alphaVector.push_back(alphaLi);
    }
    
    cout << "wmTotAll = " << wmTotAll << endl;
    return alphaVector;
}




void advanceNdfNuclAllSizes(vector<double> const& dotHvector, vector<double>const& alphaVector, vector<vector< double> >& allParticles, vector<vector< double> > & ndft, double h, double nT, double a,double it, double timePerIt, vector<double> lVector)
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
    double dotHtot(0);
    double dotAt(0);
    
    //dotH = nuclSourceMass(allParticles, h);
    
    dotHtot = 0.0;
    i=0;
    for(i=0; i<dotHvector.size(); i++)
    {
        dotHtot+= dotHvector[i];
    }
    
    // calculation of sum of the wm to get dotAt. calculated with ndft. info t-deltaT
    i=0;
    for (i=0; i<ndft.size(); i++)
    {
        double li = lVector[i];
        double lavgi = li*1.125;
        dotAt += wmTotj(i, ndft, timePerIt, lVector, a) / lavgi;  // [part/volume/time]
    }
    cout << "dotAt sum of Ali = " << dotAt << endl;
    double alphaH = dotHtot/nT;
    double alphaAt = dotAt/nT;
    
    double alphaHplusAt = alphaH + alphaAt;
    
    
    // calculation of np(l*, t+dt)   (called nplDt here) and storing of the corresponding  Delta np = np(l*,t+dt) - np(l*,t) in a new vector: deltaNpInt. New vector because before we were storing integers in a vector of doubles (lNplNvl)
    
    
    vector<int> deltaNpInt;
    
    i = 0;
    int deltaNplSum(0);
    for(i=0; i<ndft.size(); i++)
    {
        double nplDt(0);
        nplDt = ndft[i][1] * (1 - alphaHplusAt) + alphaVector[i] * Np;  //calculation of np(l*,t+dt)
        
        int roundednp2 = rounding(nplDt);                 // rounding function
        int roundednp1 = rounding(ndft[i][1]);         // rounding function
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
    
    for(i=0; i<ndft.size(); i++)
    {
        li = ndft[i][0];
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
                
                if(allParticles[j][1]>=infborn && allParticles[j][1]<supborn)
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
