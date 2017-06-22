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


std::vector<double> initParticles(int Np0, int Np1, double c0, double c1);

void mix(std::vector<double>& allParticles, double deltaT, double tau, double t);


void printParticles(std::vector<double> const& allParticles, double t);

#include <stdio.h>

#endif /* mixing_hpp */
