//
//  geometric.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 21/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "geometric.hpp"
#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

vector<vector<double> > geolNplNv(vector<vector<double> > allParticles, vector<double> liVector, double deltaL, double nT, double lp0, double maxvalL)
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
        li = lVector[i];
        double infborn(1);
        double supborn(2);
        
        // infborn and supborn calculation for each li
        if(i==0)
        {
            infborn = li-(liVector[i+1]-li)/2.0;
            supborn = li+(liVector[i+1]-li)/2.0;
        }
        else if (li == maxvalL)
        {
            infborn = li-(li-liVector[i-1])/2.0;
            supborn = li+(li-liVector[i-1])/2.0;
        }
        else
        {
            infborn = li-(li-liVector[i-1])/2.0;
            supborn = li+(liVector[i+1]-li)/2.0;
        }
        
        
        // count of np(li)
        j=0;
        for(j=0;j<allParticles.size();j++)
        {
            if((allParticles[j][1]>=infborn)&(allParticles[j][1]<supborn))
            {
                npL++;
            }
        }
        
        nvL = double(npL)/double(Np)*nT;
        lAndNpl.push_back(vector<double>(3,0));
        lAndNpl[i][0] = li;
        lAndNpl[i][1] = npL;
        lAndNpl[i][2] = nvL;
        
        /*
        cout << "borninf = " << infborn << "   supborn = " << supborn << endl;
        cout << "deltaL["<<i<<"] = " << (supborn - infborn) << endl;
        cout << "lAndNpl["<<i<<"][1] = " << npL << endl;
        cout << "nv["<<i<<"] = " << nvL << endl;
         */
        
    }
    return lAndNpl;
}

