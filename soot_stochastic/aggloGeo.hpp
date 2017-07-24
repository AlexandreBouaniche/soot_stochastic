//
//  aggloGeo.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 24/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef aggloGeo_hpp
#define aggloGeo_hpp

#include <vector>



double nvLstarGeo(double lStar, std::vector<std::vector<double> > lNplNvl, double maxValL);


double dotAlStarGeo(double lStar, std::vector<std::vector< double> > const& allParticles, std::vector<std::vector< double> > const& lNplNvl, double a, double deltaL, double nT, double timePerIt, double maxValL);

std::vector<double> allAlphaCoefGeo(std::vector<std::vector< double> > const& allParticles, double lp0, double a, double nT, double nTtminusOne, double h, double deltaL, std::vector<std::vector< double> > const& lNplNvl, double t, double timePerIt, double maxValL);

void advancePdfGeo(std::vector<double>const& alphaVector, std::vector<std::vector< double> >& allParticles, std::vector<std::vector< double> > & lNplNvl, double h, double nT, double a, double deltaL, double t, double maxValL, double lp0, double nTtminusOne, double timePerIt);

#include <stdio.h>

#endif /* aggloGeo_hpp */
