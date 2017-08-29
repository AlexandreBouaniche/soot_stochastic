//
//  GasParticle.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "GasParticle.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

GasParticle::GasParticle() : SootStochasticParticle() , m_c(0.0)
{
    
}

GasParticle::GasParticle(double massSoot) : SootStochasticParticle(massSoot) , m_c(0.0)
{
    
}

GasParticle::GasParticle(double massSoot, double c) : SootStochasticParticle(massSoot), m_c(c)
{
    
}


double GasParticle::getC()
{
    return m_c;
}


void GasParticle::setC(double c)
{
    m_c = c;
}
