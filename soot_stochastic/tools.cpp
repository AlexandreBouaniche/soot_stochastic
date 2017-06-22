//
//  tools.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "tools.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// for random function
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

/*
 Picks random integer values. Returns a vector of size nbSorted of random integers with values between 0 and maxVal INCLUDED. 
 
 DEPENDS ON t!! For the same t, picks the same list. (Random within the list but two identical lists if calculated during the same run without dependency on t)
 
 nbPicked must be an EVEN number
 
 */
vector<int> randomList(double t, int nbPicked, int maxVal)
{
    int random_integer = 0;
    vector<int> randomList;
    srand(time(NULL)+t); // initialisation de rand
    int i;
    for(i=0; i<nbPicked; i++)
    {
        random_integer = rand()%(maxVal+1);
        randomList.push_back(random_integer);
    }
    return randomList;
}
