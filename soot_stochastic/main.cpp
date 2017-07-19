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
    int Np0 = 500;                // initial stochastic particles at inlet0
    int Np1 = 500;                // initial stochastic particles at inlet1
    int Np2 = 500;
    int Np3 = 500;
    
    double c0 = 0.0;              // initial progress variable at inlet0
    double c1 = 1.0;              // initial progress variable at inlet1
    
    double lp0 = 0.02;             // nascent particles size
    
    double l0 = 0.01;              // initial soot size at inlet0
    double l1 = 0.01;             // initial soot size at inlet1
    double l2 = 0.01;
    double l3 = 0.01;
    
    int itTot = 50;                 // number of iteration
    
    
    double pdfGrid(0.1);    // distance between two c bins for graphic representation of P(c)
    double LpdfGrid(0.02);      // distance between two l bins for graphic representation of P(l)
    double maxValC(1);       // maximum value of c considered for graphic representation of P(c)
    double maxValL(2.0);     // maximum value of l considered for graphic representation of P(l)
    double minValC(0);      // minimum value of c for graphic representation of P(c)
    double deltaL(0.02);        // spacing between two intervals Il*
    
    
    // model parameters
    double h = 1.0;              // constant used for source term of nucleation
    double a = 0.0;                 // constant used for source term of agglomeration
    double nT0 = 1900;             // initial total soot number density
    double uniformG = 0.01;
    //double linearG = 0.02;
    //double linearOxi = -0.0005;
    //double ageFactor = 20;
    
    // time and mixing parameters
    double time(0);                 // time
    double timePerIt = 0.01;        // time per iteration
    //double tau(2);                // characteristic mixing time as a function of iterations.
    
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
    
    //allParticles = initAllParticles(initVector);
    allParticles = initCustomized();
    
    
    // Initial state. write in output files
    double it(0);
    double nT = nT0;
    double nTtminusOne = nT;
    writePdft(pathProject, "/outputs/Cpdf_t/Cpdf", it, allParticles, pdfGrid, minValC, maxValC, 0);
    writePdft(pathProject, "/outputs/Lpdf_t/Lpdf", it, allParticles, LpdfGrid, lp0, maxValL, 1);
    writeNvt(pathProject, "/outputs/Nv_t/Nv", it, allParticles, LpdfGrid, lp0, maxValL, 1, nT);
    vector<double> lVector = liVector(lp0, deltaL, maxValL);  // vector with all the li from lp0 to maxValL. 
    
    //printParticles(allParticles, t);
    writeCustomNv(pathProject, "/outputs/Nv_t/NvRef", 0, allParticles, 0.02, 0.005, maxValL, 1, nT);
    
    // advancing t, mixing (Cpdf), source terms, advancing nT and Lpdf
    int j;
    for(j=0; j<itTot; j++ )
    {
        
        
        vector<vector<double> > lAndNpL;
        lAndNpL = liNpliNvli(allParticles, lVector, deltaL, nT);  // col0: li; col1: npli; col2: nvli. calculated BEFORE ADVANCING nT to nT(t+deltat) ! doesn't "see" particles out of bounds (lp0 and maxValL)
        
        int i=0;
        double li(0.02);
        cout << "time = " << time<<endl;
        for(i=0; i<100; i++)
        {
            int countAllPartLi(0);
            int j(0);
            for(j=0; j<allParticles.size(); j++)
            {
                if(allParticles[j][1]>=(li-0.02/2.0) & allParticles[j][1]<(li+0.02/2.0))
                {
                    countAllPartLi++;
                }
            }
            cout << "coutSmallGrid("<<li << ") = " << countAllPartLi << "   ;";
            li += 0.02;
        }
        cout << endl << endl;
        
        int countNextOut(0);
        int j = 0;
        for(j=0; j<allParticles.size(); j++)
        {
            if(allParticles[j][1]>=2.0 & allParticles[j][1]<2.01)
            {
                countNextOut++;
            }
        }
        cout << "countNextOut = " << countNextOut << endl << endl;

        
        
        
        it = it + 1;                                         // advancing t
        time = it*timePerIt;
        //mix(allParticles, deltaT, tau, t);                    // advancing Cpdf = mixing
        uniformGrowth(allParticles, uniformG);
        //linerarSurfGrowth(allParticles, linearG, lp0, maxValL, deltaL);          // growth proportional to the surface
        //linerarSurfOxi(allParticles, linearOxi, lp0, deltaL);          // oxidation proportional to the surface
        //surfGrowthAging(allParticles, linearG, lp0, maxValL, deltaL, ageFactor);
        
        
        double dotH = nuclSourceCustomized(it, deltaL);          // calculation of nucleation source term
        double dotAt = aggloTotSource(allParticles, lAndNpL, a);  // calculation of total agglomeration source term
        double dotG = outOfBoundSource(allParticles, nT, maxValL, lp0, deltaL);
        
        vector<vector<double> > lAndNpLg;            // updated vector lAndNpLg after growth before nT = nT(t+dt) and advancing pdf
        lAndNpLg = liNpliNvli(allParticles, lVector, deltaL, nT);
        
        advanceGrowthPdf(allParticles, nT, maxValL, lp0, deltaL, lAndNpLg);
        writeNvt(pathProject, "/outputs/Nv_tg/NvG", it, allParticles, LpdfGrid, lp0, maxValL, 1, nT);
        
        nTtminusOne = nT;                                    // storing nT(t-deltat)
        nT = nT + dotH +dotAt + dotG;                                 // advancing nT
        
        cout << endl << endl;
        cout << "it = " << it << "   time = " << time << endl;
        cout << "nT = " << nT << "   " << "dotH = " << dotH << "   dotAt = " << dotAt << endl;
        cout << "dotG = " << dotG << endl;
        
        
        vector<double> alphaVector = allAlphaCoef(allParticles, lp0, a, nT,nTtminusOne, h, deltaL, lAndNpL, it);  // coefs used for advancePdf
        
         
        advancePdf(alphaVector, allParticles, lAndNpLg, h, nT, a, deltaL, it, maxValL, lp0, nTtminusOne);
        
        //printParticles(allParticles, t);
        writePdft(pathProject, "/outputs/Cpdf_t/Cpdf", it, allParticles, pdfGrid, minValC, maxValC, 0);
        writePdft(pathProject, "/outputs/Lpdf_t/Lpdf", it, allParticles, LpdfGrid, lp0, maxValL, 1);
        writeNvt(pathProject, "/outputs/Nv_t/Nv", it, allParticles, LpdfGrid, lp0, maxValL, 1, nT);
    }
    
    
    return 0;
}
