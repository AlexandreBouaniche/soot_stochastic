//
//  SootStochasticParticle.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "SootStochasticParticle.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

SootStochasticParticle::SootStochasticParticle()
{
    
}

SootStochasticParticle::SootStochasticParticle(double mass) : m_mass(mass)
{
    
}

double SootStochasticParticle::getMass()
{
    return m_mass;
}

void SootStochasticParticle::setMass(double mass)
{
    m_mass = mass;
}


void SootStochasticParticle::addMass(double mass)
{
    m_mass+=mass;
}
