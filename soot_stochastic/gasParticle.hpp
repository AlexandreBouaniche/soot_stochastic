//
//  GasParticle.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef GasParticle_hpp
#define GasParticle_hpp

#include <stdio.h>
#include <vector>

#include "SootStochasticParticle.hpp"

class GasParticle : public SootStochasticParticle
{
public:
    GasParticle();
    GasParticle(double massSoot);
    GasParticle(double massSoot, double c);
    
    double getC();
    void setC(double c);
    
protected:
    double m_c;
    // cantera reactor / idealGasMix / vector<double> Yi ...
};

#endif /* GasParticle_hpp */
