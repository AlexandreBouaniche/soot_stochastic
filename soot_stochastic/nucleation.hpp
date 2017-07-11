//
//  nucleation.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 23/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef nucleation_hpp
#define nucleation_hpp

#include <vector>

double nuclSource(std::vector<std::vector<double> > const& allParticles, double h);

void LpdfAlphaH(std::vector<std::vector< double> >& allParticles, double nT, double dotH, double lp0, double t);

double nuclSourceCustomized(double it);

#include <stdio.h>

#endif /* nucleation_hpp */
