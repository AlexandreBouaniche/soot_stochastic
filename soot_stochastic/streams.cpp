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
#include <sstream>

using namespace std;


void writePdft(string pathProject, string pathTarget, int it, vector<vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column)
{
    string finalPath = pathProject.append(pathTarget);
    
    int tInt = it;
    stringstream ss;
    ss << tInt;
    string strIt = ss.str();
    
    finalPath = finalPath.append(strIt);
    string dat = ".dat";
    finalPath = finalPath.append(dat);
    
    ofstream stream1(finalPath.c_str());
    if(stream1) // error test
    {
        //cout << "stream OK" << endl;
        
        stream1 << "#iteration number = " << it << " pdf bins vertically in column 1"<<endl;
        int j;
        double c=minVal;
        double cGrid;
        cGrid = maxVal/pdfGrid;
        int intGrid;
        double count(0);
        int i(0);
        intGrid = floor(cGrid)+1;
        for(j=0; j<intGrid; j++)
        {
            stream1 << c << "   ";
            
            count = 0;
            for(i=0; i<allParticles.size(); i++)
            {
                if((allParticles[i][column]>=(c-pdfGrid/2))&(allParticles[i][column]<(c+pdfGrid/2)))
                {
                    count++;
                }
            }
            double Pc;
            Pc = count/allParticles.size();
            stream1 << Pc << endl;

            c = c+pdfGrid;
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}

