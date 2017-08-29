//
//  Psd.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef Psd_hpp
#define Psd_hpp

#include <stdio.h>
#include <vector>

#include "Bin.hpp"
#include "Grid.hpp"
#include "GridCell.hpp"
#include "AerosolPhase.hpp"

class Psd
{
public:
    Psd(Grid initialGrid, AerosolPhase* aeroSolPtr);
    Psd(Grid initialGrid, AerosolPhase* aeroSolPtr, double initNvT);
    Psd(Grid initialGrid, AerosolPhase* aeroSolPtr, double initNvT, std::vector<int> initNpVector);  // should not be needed if np is based on counting in AerosolPhase.
    
    int getSize();
    
    double getNvT();
    int getNpT();
    
    Bin getBin(int i);
    Bin* getBinPtr(int i);      // returns a pointer to Bin i
    
    double getmAvg(int i);      // [mass/particle]  average mass of a particle in Bin i
    double getmInf(int i);      // [mass/particle]  min mass of a particle in Bin i
    double getmSup(int i);      // [mass/particle]  max mass of a particle in Bin i
    
    int getNp(int i);        // [part] number of stochastic particles in Bin i
    double getNv(int i);        // [part/volume]
    double calculateNd(int i);  // [part/volume/xcoordinate in size space] -> Number density function. Calculated from the Nv in the Psd bins. These need to be up to date.
    
    double calculateMv(int i);  // [mass/volume]
    double calculateMvT();      // [mass/volume]
    
    
    // for advancing the PSD in time with s2solver
    void setNvt(double nvT);         // will be set by S2solver at each time step
    
    void countAndSetAllNp();           // counts Np in AerosolPhase and sets it in all Bins. (S2solver will reallocate particles directly modifying the AerosolPhase object)
    
    void calcAndSetAllNv();          // will be calculated and set from the npi and nvT at each time step
    
    void showPsd();
    
private:
    std::vector<Bin> m_binVector;
    double m_nvT;                      //  [part/m3]  sum of Nv for all Bins
    AerosolPhase* m_aeroSolPtr;
};



#endif /* Psd_hpp */
