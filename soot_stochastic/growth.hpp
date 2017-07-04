//
//  growth.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 03/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef growth_hpp
#define growth_hpp

#include <vector>

void uniformGrowth(std::vector<std::vector<double> > &allParticles, double deltaG);

void linerarSurfGrowth(std::vector<std::vector<double> > &allParticles, double deltaM0, double lp0);

void linerarSurfOxi(std::vector<std::vector<double> > &allParticles, double deltaM0, double lp0);

double dotOxi(std::vector<std::vector<double> > const& allParticles, double lp0, double deltaL);

#include <stdio.h>

#endif /* growth_hpp */
