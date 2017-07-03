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

void writeFile(string pathProject, string pathTarget, vector<vector<double> > allParticles)
{
    string finalPath = pathProject.append(pathTarget);
    ofstream stream1(finalPath.c_str());
    if(stream1) // error test
    {
        //cout << "stream OK" << endl;
        
        stream1 << "#time   ";
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


void updateFile(string pathProject, string pathTarget, double t, vector<vector<double> > allParticles)
{
    string finalPath = pathProject.append(pathTarget);
    ofstream stream1(finalPath.c_str(), ios::app);
    if(stream1) // error test
    {
        //cout << "stream OK" << endl;
        
        stream1 << t << "   ";
        int i;
        for(i=0; i<allParticles.size(); i++)
        {
            stream1 << allParticles[i][0] << "   ";
        }
        stream1 << endl;
        
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}


void writeCpdf(string pathProject, string pathTarget, vector<vector<double> > allParticles, double cSensibility)
{
    string finalPath = pathProject.append(pathTarget);
    ofstream stream1(finalPath.c_str());
    if(stream1) // error test
    {
        //cout << "stream OK" << endl;
        
        stream1 << "#time vertically, c horizontally   " <<endl;
        stream1 << "   #c=";
        int i;
        double c(0.0);
        double cGrid;
        cGrid = 1/cSensibility;
        int intGrid;
        intGrid = floor(cGrid)+1;
        for(i=0; i<intGrid; i++)
        {
            stream1 << c << "   ";
            c = c+cSensibility;
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}


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



void updateCpdf(string pathProject, string pathTarget, double t, vector<vector<double> > allParticles, double cSensibility)
{
    string finalPath = pathProject.append(pathTarget);
    ofstream stream1(finalPath.c_str(), ios::app);
    if(stream1) // error test
    {
        //cout << "stream OK" << endl;
        
        stream1 << t << "   ";
        int i;
        double count(0);
        
        double cGrid;
        cGrid = 1/cSensibility;
        int intGrid;
        intGrid = floor(cGrid)+1;
        int j;
        
        double c(0.0);
        for(j=0; j<intGrid; j++)
        {
            count = 0;
            for(i=0; i<allParticles.size(); i++)
            {
                
                if((allParticles[i][0]>=(c-cSensibility/2))&(allParticles[i][0]<(c+cSensibility/2)))
                {
                    count++;
                }
            }
            
            double Pc;
            Pc = count/allParticles.size();
            
            stream1 << Pc << "   ";
            c = c+cSensibility;
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}

