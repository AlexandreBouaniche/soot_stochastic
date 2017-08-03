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
#include <cmath>

#include "mixing.hpp"
#include "tools.hpp"
#include "streams.hpp"
#include "nucleation.hpp"
#include "agglomeration.hpp"
#include "growth.hpp"
#include "init.hpp"
#include "aggloGeo.hpp"
#include "aggloGeo2.hpp"
#include "advanceMass.hpp"
#include "advanceNdf.hpp"

using namespace std;

int main()
{
    cout << "Soot Stochastic" << endl << endl;
    
    //user inputs
    string pathProject("/Users/bouaniche/Xcode_projects/soot_stochastic");
    
    double lp0 = 0.1;             // nascent particles size
    int itTot = 10000;                 // number of iteration
    double pdfGrid(0.1);    // distance between two c bins for graphic representation of P(c)
    double LpdfGrid(0.02);      // distance between two l bins for graphic representation of P(l)
    double maxValC(1);       // maximum value of c considered for graphic representation of P(c)
    
    // Parameters for geometric grid construction
    int nBins = 30;
    double geoQ = 2.0;                     // must be >= 2 for geo2mesh
    
    double maxValL=lp0*pow(geoQ, double(nBins-1));     // maximum value of l considered for graphic representation of P(l).
                             // BE CAREFUL: if geo2mesh is used, please write a compatible maxvalL e.g lp0*geoQ^(Nbins-1)
    double minValC(0);      // minimum value of c for graphic representation of P(c)
    //double deltaL(0.02);        // spacing between two intervals Il*
    
    
    // model parameters
    double h = 0.0;              // constant used for source term of nucleation
    double a = 1.0;                 // constant used for source term of agglomeration
    double nT0 = 2.9297;             // initial total soot number density
    
    double dotAt(0);
    
    /*
    //mass parameters
    //double mT0 = 0.8978;
    double mT0 = 0.9606;
    double mT = mT0;
    double mTminusOne = mT;
    double dotAtm(0);
     */
    
    double uniformG = 0.0;
    //double linearG = 0.02;
    //double linearOxi = -0.0005;
    //double ageFactor = 20;
    
    // time and mixing parameters
    double time(0);                 // time
    double timePerIt = 1e-2;        // time per iteration
    //double tau(2);                // characteristic mixing time as a function of iterations.
    
    
    
    // construction of the mesh for the bins. Regular or geometric
    //vector<double> lVector = liVector(lp0, deltaL, maxValL);  // vector with all the li from lp0 to maxValL.
    //vector<double> lVector = initGeoMesh(lp0, maxValL, nBins, geoQ);
    vector<double> lVector = initGeo2Mesh(lp0, nBins, geoQ);
    
    
    int i(0);
    cout << endl;
    for(i=0; i<lVector.size();i++)
    {
        cout << "l["<<i<<"] = "<<lVector[i] << "   ";
    }
    cout << endl;
    i=0;
    
    // initiate particles
    vector<vector<double> > allParticles;  // col0: ci; col1: li
    
    
    //Customized Init
    
    //allParticles = initAllParticles(initVector);
    allParticles = initCustomAggloMass(lVector);
    
    
    
    
    
    // Initial state. write in output files
    double it(0);
    double nT = nT0;
    double nTtminusOne = nT;
    
    
    // write analytical Ref customized
    //writeCustomNv(pathProject, "/outputs/Nv_t/NvRef_0.0025_", 0, allParticles, 0.0025, 0.02, maxValL, 1, nT);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 0.0);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 0.05);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 0.2);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 1.0);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 2.0);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 3.0);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 5.0);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 20.0);
    writeCustomAggloCase(pathProject, "/outputs/Nd_t/Ref_t", lVector, 100.0);
    
    // advancing t, mixing (Cpdf), source terms, advancing nT and Lpdf
    vector<vector<double> > ndft;
    int j;
    for(j=0; j<itTot; j++ )
    {
        //printParticles(allParticles, it);
        
        
        //lAndNpL = liNpliNvli(allParticles, lVector, deltaL, nT);      // col0: li; col1: npli; col2: nvli. calculated BEFORE ADVANCING nT to nT(t+deltat) ! doesn't "see" particles out of bounds (lp0 and maxValL)
        //lAndNpL = geo2lNplNv(allParticles, lVector, nT, lp0, maxValL);
        
        
        // storing nT(t-deltat)
        nTtminusOne = nT;
        
        // values it0 and then itMinus1
        ndft = ndf(allParticles, lVector, nT);
        
        
        double totalMassBins(0);
        
        
        // printing for it-1
        cout << endl << endl;
        cout << "it = " << it << "   time = " << time << endl;
        totalMassBins = totalMassNdf(ndft);
        
        cout << "nT = " << nT << endl;
        cout << "dotAt = " << dotAt << endl;
        cout << "total mass bins = " << totalMassBins << endl;
        
        //cout << "dotH = " << dotH << "   dotAt = " << dotAt << endl;
        //cout << "dotG = " << dotG << endl;

        
        //writing for it-1
        writePdft(pathProject, "/outputs/Cpdf_t/Cpdf", it, allParticles, pdfGrid, minValC, maxValC, 0, ndft);
        writePdft(pathProject, "/outputs/Lpdf_t/Lpdf", it, allParticles, LpdfGrid, lp0, maxValL, 1, ndft);
        writeGeoNdt(pathProject, "/outputs/Nd_t/Nd", it, allParticles, LpdfGrid, lp0, maxValL, 1, nT,ndft);
        
        
        
        it = it + 1;                                         // advancing t
        time = it*timePerIt;
        //mix(allParticles, deltaT, tau, t);                    // advancing Cpdf = mixing
        uniformGrowth(allParticles, uniformG);
        //linerarSurfGrowth(allParticles, linearG, lp0, maxValL, deltaL);          // growth proportional to the surface
        //linerarSurfOxi(allParticles, linearOxi, lp0, deltaL);          // oxidation proportional to the surface
        //surfGrowthAging(allParticles, linearG, lp0, maxValL, deltaL, ageFactor);
        
        
        //double dotH = nuclSourceCustomized(it, deltaL);          // calculation of nucleation source term
        //double dotH = nuclSource(allParticles, h);         // calculation of nucleation source term
        double dotH =0;
        dotAt = aggloTotSource(allParticles, ndft, a, timePerIt);   // calculation of total agglomeration source term
        
        
        
        //double dotG = outOfBoundSource(allParticles, nT, maxValL, lp0, deltaL);
        double dotG = 0;
        
        //vector<vector<double> > lAndNpLg;            // updated vector lAndNpLg after growth before nT = nT(t+dt) and advancing pdf
        //lAndNpLg = liNpliNvli(allParticles, lVector, deltaL, nT);
        
        //advanceGrowthPdf(allParticles, nT, maxValL, lp0, deltaL, lAndNpLg);
        writeNvt(pathProject, "/outputs/Nv_tg/NvG", it, allParticles, LpdfGrid, lp0, maxValL, 1, nT, ndft);
        
        
        nT = nT + dotH +dotAt + dotG;                                 // advancing nT
        
        
        vector<double> alphaVector = allAlphaCoefNdf(allParticles, a, nT, h, ndft, timePerIt, lVector);  // coefs used for advancePdf
        
        advanceNdf(alphaVector, allParticles, ndft, h, nT, a, it, timePerIt, lVector);
        
    }
    
    // count of particles one more time for the final step
    //lAndNpL = liNpliNvli(allParticles, lVector, deltaL, nT);
    //lAndNpL = geo2lNplNv(allParticles, lVector, nT, lp0, maxValL);
    ndft = ndf(allParticles, lVector, nT);
    
    //printParticles(allParticles, it);
    
    writePdft(pathProject, "/outputs/Cpdf_t/Cpdf", it, allParticles, pdfGrid, minValC, maxValC, 0, ndft);
    writePdft(pathProject, "/outputs/Lpdf_t/Lpdf", it, allParticles, LpdfGrid, lp0, maxValL, 1, ndft);
    writeGeoNdt(pathProject, "/outputs/Nd_t/Nd", it, allParticles, LpdfGrid, lp0, maxValL, 1, nT, ndft);
    
    
    return 0;
}
