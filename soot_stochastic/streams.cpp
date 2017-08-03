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


void writePdft(string pathProject, string pathTarget, int it, vector<vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, vector<vector<double> > lAndNpl)
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
        double count(0);
        
        for(j=0; j<lAndNpl.size(); j++)
        {
            stream1 << c << "   ";
            
            count = lAndNpl[j][1];
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


void writeNvt(string pathProject, string pathTarget, int it, vector<vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, double nT, vector<vector<double> > lAndNpl)
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
        double li=lAndNpl[0][0];
        double count(0);
        for(j=0; j<lAndNpl.size(); j++)
        {
            li = lAndNpl[j][0];
            stream1 << li << "   ";
            
            count = lAndNpl[j][1];
            double Pc;
            Pc = count/allParticles.size();
            double nv(0);
            nv = Pc*nT;
            stream1 << nv << endl;
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}



void writeCustomNv(string pathProject, string pathTarget, int it, vector<vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, double nT)
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
        double c=0.0025;
        double time(0.5);
        double nv(0);
        for(j=0; j<200; j++)
        {
            stream1 << c << "   ";
            nv = 100+1000*exp(-10000*pow((time-0.215),2));
            stream1 << nv << endl;
            c = c+pdfGrid;
            time = time - 0.0025;
        }
        
        j=0;
        for(j=0; j<160; j++)
        {
            stream1 << c << "   ";
            nv = 10;
            stream1 << nv << endl;
            c = c+pdfGrid;
        }
        
        j=0;
        for(j=0; j<80; j++)
        {
            stream1 << c << "   ";
            nv = 100;
            stream1 << nv << endl;
            c = c+pdfGrid;
        }
        
        j=0;
        for(j=0; j<390; j++)
        {
            stream1 << c << "   ";
            nv = 10;
            stream1 << nv << endl;
            c = c+pdfGrid;
        }
        
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}



void writeGeoNdt(string pathProject, string pathTarget, int it, vector<vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, double nT, vector<vector<double> > lAndNpl)
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
        double li=lAndNpl[0][0];
        for(j=0; j<lAndNpl.size(); j++)
        {
            li = lAndNpl[j][0];
            double lavg = li*1.125;
            double deltaLinti = 1.5*li - 0.75*li;
            stream1 << lavg << "   ";
            double nv;
            double nd;
            nv = lAndNpl[j][2];    // [part/volume]
            nd = nv / deltaLinti;  // [part/volume/deltaLint]
            stream1 << nd << endl;
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}


void writeCustomAggloCase(string pathProject, string pathTarget, vector<double> lVector, double time)
{
    string finalPath = pathProject.append(pathTarget);
    stringstream ss;
    ss << time;
    string strIt = ss.str();
    
    finalPath = finalPath.append(strIt);
    string dat = ".dat";
    finalPath = finalPath.append(dat);
    
    ofstream stream1(finalPath.c_str());
    if(stream1) // error test
    {
        //cout << "stream OK" << endl;
        
        stream1 << "#time = " << time << " pdf bins vertically in column 1"<<endl;
        
        double beta0 = 1.0;
        int j;
        for(j=0; j<lVector.size(); j++)
        {
            double li = lVector[j];
            double lavg = li*1.125;
            stream1 << lavg << "   ";
            
            double solutioni(0);
            solutioni = 4 / pow((beta0*time+2),2.0) * exp(-2*lavg/(beta0*time+2));
            stream1 << solutioni << endl;
        }
        stream1 << endl;
    }
    else
    {
        cout << "ERROR: Impossible to open the file." << endl;
    }
}
