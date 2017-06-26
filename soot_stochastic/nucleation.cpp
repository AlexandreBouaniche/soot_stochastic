//
//  nucleation.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 23/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "nucleation.hpp"
#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

double nuclSource(vector<vector< double> > const& allParticles, double h)
{
    int i(0);
    double dotH(0);
    for(i=0; i<allParticles.size(); i++)
    {
        double ci(0);
        ci = allParticles[i][0];
        //cout << "ci = " << ci << "   ";
        dotH = dotH + h * ci * pow((1-ci),5);
        //cout << "dotHi = " << dotH << "   ";
    }
    //cout << endl;
    return dotH;
}

void LpdfAlphaH(vector<vector< double> >& allParticles, double nT, double dotH, double deltaL, double lp0, double t)
{
    double alphaH = dotH / nT;
    int npL0(0);                             // declaration: number of particles of size l0
    int Np(0);
    int i;
    double li(0);
    
    cout << "alphaH = " << alphaH << endl;
    
    
    for(i=0; i< allParticles.size(); i++)
    {
        Np++;                                // count of Np
        li = allParticles[i][1];
        if(li < (lp0+deltaL))                 // count of n(lo, t)
            npL0++;
    }
    
    //cout << "npL0 = " << npL0 << endl;
    
    int npLstar;
    npLstar = Np - npL0;                      // expression for n(l*) for l* != l0
    
    //cout << "npL* = " << npLstar << endl;
    
    int nbToPick(0);
    nbToPick = floor(npLstar*alphaH);         // number of particles to switch from l* to l0
    //cout << "nbToPick = " << nbToPick << endl << endl;
    
    vector<vector<double> > tempLstar;
    int j;
    int itemp(0);                            // itemp: int to count size of tempLstar for randomList function argument after
    for(j=0; j<allParticles.size(); j++)
    {
        if(allParticles[j][1]>=(lp0+deltaL))
        {
            // temporary vector to separate partcles of size > (l0+deltaL)   putting rank in the vector in the first column and size in the second column
            tempLstar.push_back(vector<double>(2,0));
            tempLstar[itemp][0] = j;
            tempLstar[itemp][1] = allParticles[j][1];
            
            itemp++;
        }
    }
    
    itemp = itemp-1;                      // -1 because randomList takes values from 0 to maxVal INCLUDED. with -1 EXCLUDED (corresponds to size of tempLstar vector)

    vector<int> randomL = randomList((t+lp0), nbToPick, itemp);   // random pick of particles to be reset to l0. parameter time chosen as t+lp0 so that it is different for each interval l* pick. and for each time iteration t.
    
    j=0;
    int rankAll(0);                                    // rank in AllParticles vector
    int rankTemp(0);                                   // rank in tempLstar
    for(j=0; j<randomL.size(); j++)
    {
        rankTemp = randomL[j];
        rankAll = tempLstar[rankTemp][0];
        allParticles[rankAll][1] = lp0;
    }
    
}
