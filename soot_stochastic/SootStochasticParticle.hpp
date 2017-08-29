//
//  SootStochasticParticle.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef SootStochasticParticle_hpp
#define SootStochasticParticle_hpp

#include <stdio.h>
#include <vector>

class SootStochasticParticle
{
public:
    SootStochasticParticle();
    SootStochasticParticle(double mass);
    
    double getMass();
    void setMass(double mass);
    void addMass(double deltaMass);
    
protected:
    double m_mass;
};

#endif /* SootStochasticParticle_hpp */
