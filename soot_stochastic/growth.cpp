//
//  growth.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 03/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "growth.hpp"
#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

void uniformGrowth(vector<vector<double> > &allParticles, double deltaG)
{
    int i(0);
    for(i=0; i<allParticles.size(); i++)
    {
        allParticles[i][1] += deltaG;
    }
}
