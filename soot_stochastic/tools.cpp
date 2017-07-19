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
 Picks random integer values. Returns a vector of size nbPicked of random integers with values between 0 and maxVal INCLUDED.
 
 DEPENDS ON t!! For the same t, picks the same list. (Random within the list but two identical lists if calculated during the same run without dependency on t)
 
 nbPicked must be an EVEN number for mixing
 
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


vector<int> randomListWithoutDuplicate(double t, int nbPicked, int maxVal)
{
    int random_integer = 0;
    vector<int> randomList;
    srand(time(NULL)+t); // initialisation de rand
    int i;
    for(i=0; i<nbPicked; i++)
    {
        random_integer = rand()%(maxVal+1);
        
        while(isIntInlist(random_integer, randomList))
            random_integer = rand()%(maxVal+1);
        
        randomList.push_back(random_integer);
    }
    return randomList;
}



double maxColi(vector<vector<double> > const& matrix, int col)
{
    double max(0);
    max = matrix[0][col];
    int i(0);
    for(i=0; i<matrix.size(); i++)
    {
        if(max<matrix[i][col])
            max = matrix[i][col];
    }
    return max;
}


int rounding(double d)
{
    double dSuperior = d+0.5;
    int rounded = floor(dSuperior);
    return rounded;
}
 

bool isIntInlist(int rank, vector<int> list)
{
    int i=0;
    bool is = 0;
    for(i=0; i< list.size(); i++)
    {
        if(list[i] == rank)
            is = 1;
    }
    return is;
}


double frand_a_b(double a, double b)
{
    return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}
