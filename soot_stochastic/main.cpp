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
#include <fstream>

#include "mixing.hpp"
#include "tools.hpp"
#include "streams.hpp"

using namespace std;

int main()
{
    cout << "Soot Stochastic" << endl << endl;
    
    //user inputs
    string pathProject("/Users/bouaniche/Xcode_projects/soot_stochastic");
    int Np0 = 100;
    int Np1 = 100;
    double c0 = 0.0;
    double c1 = 1.0;
    int it = 100;
    
    // time and mixing parameters
    double deltaT(1);
    double tau(5);
    
    // initiate particles
    vector<double> allParticles;
    allParticles = initParticles(Np0, Np1, c0, c1);
     
    // advancing t, mixing, printing and writing an output
    double t(0);
    int j;
    writeFile(pathProject, "/outputs/trajectory.dat", allParticles);
    for(j=0; j<it; j++ )
    {
        t = t+deltaT;
        mix(allParticles, deltaT, tau, t);
        printParticles(allParticles, t);
        updateFile(pathProject, "/outputs/trajectory.dat", t, allParticles);
        
    }
    
    
    return 0;
}
