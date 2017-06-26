//
//  agglomeration.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 26/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef agglomeration_hpp
#define agglomeration_hpp

#include <vector>

int geometricNpLstar(double lStar, std::vector<std::vector<double> > allParticles);

double geometricNLstar(double lStar, std::vector<std::vector<double> > allParticles, double nT);

double beta(double l1, double l2, double c);   // should depend on T and Knudsen not on c

double geoAggloTotSource(std::vector<std::vector< double> > const& allParticles, double maxValL, double lp0, double a, double nT);


#include <stdio.h>

#endif /* agglomeration_hpp */
