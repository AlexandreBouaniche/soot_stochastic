//
//  main.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>

#include "mixing.hpp"
#include "tools.hpp"

using namespace std;

int main()
{
    cout << "Soot Stochastic" << endl << endl;
    
    //user inputs
    int Np0 = 10;
    int Np1 = 10;
    double c0 = 0.0;
    double c1 = 1.0;
    int it = 50;
    
    // time and mixing parameters
    double deltaT(1);
    double tau(5);
    
    // initiate particles
    vector<double> allParticles;
    allParticles = initParticles(Np0, Np1, c0, c1);
     
    // advancing t, mixing and printing
    double t(0);
    int j;
    for(j=0; j<it; j++ )
    {
        t = t+deltaT;
        mix(allParticles, deltaT, tau, t);
        printParticles(allParticles, t);
    }
    
    
     
    return 0;
}
