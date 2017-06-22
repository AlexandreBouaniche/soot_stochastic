//
//  streams.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "streams.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

void writeFile(string pathProject, string pathTarget, vector<double> allParticles)
{
    string finalPath = pathProject.append(pathTarget);
    ofstream stream1(finalPath.c_str());
    if(stream1) // error test
    {
        cout << "stream OK" << endl;
        
        stream1 << "time   ";
        int i;
        for(i=0; i<allParticles.size(); i++)
        {
            stream1 << "p" << i << "   ";
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }

}


void updateFile(string pathProject, string pathTarget, double t, vector<double> allParticles)
{
    string finalPath = pathProject.append(pathTarget);
    ofstream stream1(finalPath.c_str(), ios::app);
    if(stream1) // error test
    {
        cout << "stream OK" << endl;
        
        stream1 << t << "   ";
        int i;
        for(i=0; i<allParticles.size(); i++)
        {
            stream1 << allParticles[i] << "   ";
        }
        stream1 << endl;
        
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }

    
}
