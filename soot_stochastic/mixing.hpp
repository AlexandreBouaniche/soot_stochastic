//
//  mixing.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef mixing_hpp
#define mixing_hpp

#include <vector>


std::vector<std::vector<double> > initAllParticles(std::vector<std::vector<double> > initVector);

void mix(std::vector<std::vector<double> >& allParticles, double deltaT, double tau, double t);


void printParticles(std::vector<std::vector<double> > const& allParticles, double t);

#include <stdio.h>

#endif /* mixing_hpp */
