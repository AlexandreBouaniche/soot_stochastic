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
#include "nucleation.hpp"
#include "agglomeration.hpp"
#include "growth.hpp"

using namespace std;

int main()
{
    cout << "Soot Stochastic" << endl << endl;
    
    //user inputs
    string pathProject("/Users/bouaniche/Xcode_projects/soot_stochastic");
    int Np0 = 5;                // initial stochastic particles at inlet0
    int Np1 = 5;                // initial stochastic particles at inlet1
    int Np2 = 5;
    int Np3 = 5;
    
    double c0 = 0.0;              // initial progress variable at inlet0
    double c1 = 1.0;              // initial progress variable at inlet1
    
    double lp0 = 1.0;             // nascent particles size
    double l0 = 10.0;              // initial soot size at inlet0
    double l1 = 20.0;             // initial soot size at inlet1
    double l2 = 30.0;
    double l3 = 40.0;
    
    int it = 10000;                 // number of iteration
    
    
    double pdfGrid(0.1);    // distance between two c bins for graphic representation of P(c)
    double LpdfGrid(1);      // distance between two l bins for graphic representation of P(l)
    double maxValC(1);       // maximum value of c considered for graphic representation of P(c)
    double maxValL(50);     // maximum value of l considered for graphic representation of P(l)
    double minValC(0);      // minimum value of c for graphic representation of P(c)
    double deltaL(1);        // spacing between two intervals Il*
    
    
    // model parameters
    double h = 0.0e4;              // constant used for source term of nucleation
    double a = 0.0;                 // constant used for source term of agglomeration
    double nT0 = 1e10;             // initial total soot number density
    //double uniformG = 0.2;
    //double linearG = 0.01;
    double linearOxi = -0.01;
    
    // time and mixing parameters
    double deltaT(1);             // iteration step time
    double tau(2);                // characteristic mixing time
    
    // initiate particles
    vector<vector<double> > allParticles;  // col0: ci; col1: li
    vector<vector<double> > initVector;    // col0: ci; col1: li; col2: Npi
    initVector.push_back(vector<double>(3,0));
    initVector[0][0] = c0;
    initVector[0][1] = l0;
    initVector[0][2] = Np0;
    
    initVector.push_back(vector<double>(3,0));
    initVector[1][0] = c1;
    initVector[1][1] = l1;
    initVector[1][2] = Np1;
    
    initVector.push_back(vector<double>(3,0));
    initVector[2][0] = c0;
    initVector[2][1] = l2;
    initVector[2][2] = Np2;
    
    initVector.push_back(vector<double>(3,0));
    initVector[3][0] = c1;
    initVector[3][1] = l3;
    initVector[3][2] = Np3;
    
    allParticles = initAllParticles(initVector);
     
    // Initial state. write in output files
    double t(0);
    double nT = nT0;
    writePdft(pathProject, "/outputs/Cpdf_t/Cpdf", t, allParticles, pdfGrid, minValC, maxValC, 0);
    writePdft(pathProject, "/outputs/Lpdf_t/Lpdf", t, allParticles, LpdfGrid, lp0, maxValL, 1);

    
    // advancing t, mixing (Cpdf), source terms, advancing nT and Lpdf
    int j;
    for(j=0; j<it; j++ )
    {
        t = t+deltaT;                                         // advancing t
        mix(allParticles, deltaT, tau, t);                    // advancing Cpdf = mixing
        //linerarSurfGrowth(allParticles, linearG, lp0, maxValL, deltaL);          // growth proportional to the surface
        linerarSurfOxi(allParticles, linearOxi, lp0, deltaL);          // oxidation proportional to the surface
        
        vector<double> lVector = liVector(lp0, deltaL, maxValL);  // vector with all the li
        vector<vector<double> > lAndNpL;
        lAndNpL = liNpliNvli(allParticles, lVector, deltaL, nT);  // col0: li; col1: npli; col2: nvli
        
        double dotH = nuclSource(allParticles, h);                // calculation of nucleation source term
        double dotAt = aggloTotSource(allParticles, lAndNpL, a);  // calculation of total agglomeration source term
        nT = nT + dotH +dotAt;                                 // advancing nT
        
        cout << "t = " << t << endl;
        cout << "nT = " << nT << "   " << "dotH = " << dotH << "   dotAt = " << dotAt << endl;
        
        
        vector<double> alphaVector = allAlphaCoef(allParticles, lp0, a, nT, h, deltaL, lAndNpL);  // coefs used for advancePdf
        
        advancePdf(alphaVector, allParticles, lAndNpL, h, nT, a, deltaL, t);
        
        //printParticles(allParticles, t);
        writePdft(pathProject, "/outputs/Cpdf_t/Cpdf", t, allParticles, pdfGrid, minValC, maxValC, 0);
        writePdft(pathProject, "/outputs/Lpdf_t/Lpdf", t, allParticles, LpdfGrid, lp0, maxValL, 1);
    }
    
    
    return 0;
}
