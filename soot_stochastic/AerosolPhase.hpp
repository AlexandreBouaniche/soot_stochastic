//
//  AerosolPhase.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef AerosolPhase_hpp
#define AerosolPhase_hpp

#include <stdio.h>
#include <vector>

#include "SootStochasticParticle.hpp"
#include "GasParticle.hpp"

class AerosolPhase
{
public:
    AerosolPhase();
    AerosolPhase(int npTot);
    AerosolPhase(int npTot, std::vector<double> massVector);
    AerosolPhase(int npTot, std::vector<double> massVector, std::vector<double> cVector);
    
    int getNpTot();
    
    // access GasParticles
    GasParticle* getGasPartPtr(int i);
    
    double getCi(int i);
    void setCi(int i, double c);
    
    double getMass(int i);
    void setMass(int i, double mass);
    void addMass(int i, double deltaMass);
    
    int countNpl(double mInf, double mSup);
    
protected:
    std::vector<GasParticle> m_partVector;
    
};

#endif /* AerosolPhase_hpp */
