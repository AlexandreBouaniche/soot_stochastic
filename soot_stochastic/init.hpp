//
//  init.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 21/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef init_hpp
#define init_hpp

#include <vector>

std::vector<std::vector<double> > initAllParticles(std::vector<std::vector<double> > initVector);

std::vector<double> initGeoMesh(double lp0, double maxValL, int nBins, double geoQ);

std::vector<std::vector<double> > initCustomNuclGrowth();

std::vector<std::vector<double> > initCustomGrowth();

std::vector<std::vector<double> > initCustomAgglo(std::vector<double> liVector, double maxValL);


#include <stdio.h>

#endif /* init_hpp */
