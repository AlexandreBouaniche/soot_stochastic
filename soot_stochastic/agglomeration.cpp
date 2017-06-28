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

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

vector<double> liVector(double lp0, double deltaL, double maxValL)
{
    int i(0);
    vector<double> lVector;
    double li(lp0);
    for(i=0; i<maxValL; i++)
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
    double realL1 = pow((6*l1/pi),(1/3));   //realL1 is the diameter the particle would have, to have a volume of "l1"
    double realL2 = pow((6*l2/pi),(1/3));
    
    
    double betaCalc(1);
    betaCalc = K*pow((1/V1+1/V2),0.5)*pow((realL1+realL2),2);
    return betaCalc;
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


vector<double> allAlphaCoef(vector<vector< double> > const& allParticles, double lp0, double a, double nT, double h, double deltaL, vector<vector< double> > const& lNplNvl)
{
    double dotH = nuclSource(allParticles, h);
    double alphaH = dotH / nT;
    double dotAl0 = dotAlStar(lp0, allParticles, lNplNvl, a, deltaL, nT);
    double alphaAl0 = dotAl0 /nT;
    
    double alphaL0 = alphaAl0 + alphaH;
    
    vector<double> alphaVector;
    alphaVector.push_back(alphaL0);
    
    int i(0);
    for(i=1; i<lNplNvl.size(); i++)     // begins at i=1 because we already calculated the first term alphaL0 (corresponds to lp0 and i = 0)
    {
        double l1 = lNplNvl[i][0];
        double dotALi = dotAlStar(l1, allParticles, lNplNvl, a, deltaL, nT);
        double alphaLi = dotALi/nT;
        alphaVector.push_back(alphaLi);
    }
    return alphaVector;
}

/*

void allAlphaPdf(vector<double> alphaVector, vector<vector< double> >& allParticles, double maxValL, double lp0, double t, double h, double a, double nT, double deltaL, vector<vector< double> > const& lNplNvl)
{
    
    double dotH(0);
    double dotAt(0);
    
    dotH = nuclSource(allParticles, h);
    dotAt = aggloTotSource(allParticles, lNplNvl, a);
    double alphaH = dotH/nT;
    double alphaAt = dotAt+nT;
    
 
    double li = lp0;
    vector<double> liVector;
    while (li<=maxValL)
    {
        liVector.push_back(li);
        li = li*2;
    }
 
    
    vector<int> particlesToReallocate;
    vector<double> ranksVectorLi;
    int nbToPick(0);
    double toPick(0);
    
    int i(0);
    for(i=0; i<lNplNvl.size(); i++)
    {
        double l1 = lNplNvl[i][0];
        
        toPick = (alphaH + alphaAt)*npLstar(l1, allParticles, deltaL);
        nbToPick = floor(toPick);             // number of particles to pick from np(l*) to put in the vector to reallocate
        
        int j=0;
    
        for(j=0; j<allParticles.size(); j++)
        {
            if ((allParticles[i][1]>=(l1*0.75)) & (allParticles[i][1] < (l1*1.5)))
            {
                ranksVectorLi.push_back(allParticles[i][0]);        // vector with ranks of particles of size li
            }
        }
        int npl1 = lNplNvl[i][1];
        
        vector<int> randomL = randomList(t+l1, nbToPick, npl1);  // countLi = ranksVectorLi.size()
        
        j=0;
        for(j=0; j<randomL.size(); j++)
        {
            int rankInLi(0);
            rankInLi = randomL[j];
            int rankGeneral(0);
            rankGeneral = ranksVectorLi[rankInLi];
            particlesToReallocate.push_back(rankGeneral);  // we have now a vector with all the ranks of the particles to reallocate
        }
        
    }
    
    i=0;
    int countNp(0);
    for (i=0; i<allParticles.size(); i++)
    {
        countNp++;                                   // count Np
    }
    
    
    i=0;
    int countReallocate(0);
    for (i=0; i<particlesToReallocate.size(); i++)
    {
        countReallocate++;                                   // count Np
    }

    
    vector<double> valuesToSet;
    i=0;
    for(i=0; i<lNplNvl.size(); i++)
    {
        toPick = alphaVector[i] * countNp;
        nbToPick = floor(toPick);
        // the particles in the vector to be reallocated were picked randomly. but then to allocate new values we don't need to pick randomly.
        double li(0);
        li = lNplNvl[i][0];
        int k(0);
        for(k=0; k< nbToPick; k++)
        {
            valuesToSet.push_back(li);   // a vector with the values of size we want to set
        }
        
    }
    
    i=0;
    int rankToSet(0);
    for (i=0; i<particlesToReallocate.size(); i++)
    {
        rankToSet = particlesToReallocate[i];
       // cout << "rankToSet = " << rankToSet << "   ";
      //  cout << "valueToSet = " << valuesToSet[i] << "   ";
        
     //   allParticles[rankToSet][1] = valuesToSet[i];      // resetting the value of li of the picked particles
    }
}*/

void advancePdf(vector<double>const& alphaVector, vector<vector< double> >& allParticles, vector<vector< double> > const& lNplNvl, double h, double nT, double a, double deltaL, double t)
{
    double dotH(0);
    double dotAt(0);
    
    dotH = nuclSource(allParticles, h);
    dotAt = aggloTotSource(allParticles, lNplNvl, a);
    double alphaH = dotH/nT;
    double alphaAt = dotAt/nT;
    
    //double alphaH = 0.1;        //for testing
    //double alphaAt = 0.1;       //for testing
    
    vector<vector<double> > lnbToReallocate;
    vector<vector<int> > rankLiAndNbToPick;     // building vectors with info of rankLi and nbToPick
    int i(0);
    for(i=0; i<lNplNvl.size(); i++)
    {
        lnbToReallocate.push_back(vector<double>(2,0));
        lnbToReallocate[i][0] = lNplNvl[i][0];
        lnbToReallocate[i][1] = (alphaH+alphaAt)*lNplNvl[i][1];
        
        rankLiAndNbToPick.push_back(vector<int>(2,0));
        rankLiAndNbToPick[i][0] = i;
        rankLiAndNbToPick[i][1] = floor(lnbToReallocate[i][1]);
    }
    
    /*
    i=0;
    for(i=0; i<lNplNvl.size(); i++)
    {
        cout << "rank li = " << rankLiAndNbToPick[i][0] << "   ";
        cout << "li = " << lnbToReallocate[i][0] << "   ";
        cout << "nbTo = " << rankLiAndNbToPick[i][1] << "   ";
    }
     */
    
    vector<int> allPartLiRanks;
    vector<int> ranksToReallocate;
    i=0;
    double li(0);
    for(i=0; i<lNplNvl.size(); i++)
    //for(i=0; i<1; i++)                 // for testing
    {
        li = lNplNvl[i][0];
        int j(0);
        for(j=0; j<allParticles.size(); j++)
        {
            if(allParticles[j][1]>=(li-deltaL/2) & allParticles[j][1]<(li+deltaL/2))
               {
                   allPartLiRanks.push_back(j);
               }
        }
        vector<int> randomL = randomList(t, rankLiAndNbToPick[i][1], lNplNvl[i][1]);
        
        j=0;
        for(j=0; j<randomL.size(); j++)
        {
            int rankToRe(0);
            int picked(0);
            picked = randomL[j];
            rankToRe = allPartLiRanks[picked];
            ranksToReallocate.push_back(rankToRe);
        }
    }
    
    /*  //print testing ranksToReallocate -> ok
    i=0;
    for(i=0; i<ranksToReallocate.size(); i++)
    {
        cout << "rankToRe = " << ranksToReallocate[i] << "   ";
        int rankPrint = ranksToReallocate[i];
        double liPrint = allParticles[rankPrint][1];
        cout << "li  = " << liPrint << "  ; ";
    }
     */
    
    //count of Np
    int Np(0);
    i=0;
    for(i=0; i<allParticles.size(); i++)
    {
        Np++;
    }
    
    vector<double> valuesToRealoc;
    i=0;
    for(i=0; i<lNplNvl.size(); i++)
    {
        double li = lNplNvl[i][0];
        double alphaLi;
        alphaLi = alphaVector[i];
        double toReset(0);
        toReset = Np*alphaLi;
        int nbToReset = floor(toReset); // number of particles "received" by Li
        
        int j=0;
        for(j=0; j<nbToReset; j++)
        {
            valuesToRealoc.push_back(li);  // filling vector valuesToRealoc with the right number of particles of size li
        }
    }
    
    /*
    //test vector valuesToRealoc  -> ok
    
    i=0;
    cout << "values to realoc  ";
    for(i=0; i<valuesToRealoc.size(); i++)
    {
        cout << valuesToRealoc[i] << "   ";
    }
    cout << endl;
     */
    
    
    // Now reallocating!
    i=0;
    for(i=0; i<ranksToReallocate.size(); i++)
    {
        int rankParticle(0);
        rankParticle = ranksToReallocate[i];
        allParticles[rankParticle][1] = valuesToRealoc[i];
    }
}

