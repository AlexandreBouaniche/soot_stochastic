//
//  aggloGeo2.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 26/07/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef aggloGeo2_hpp
#define aggloGeo2_hpp

#include <vector>

double aCoef(double lk, double lj);

double nvLstarGeo2(double lStar, std::vector<std::vector<double> > lNplNvl, double maxvalL);

double dotAlStarGeo2(double lStar, std::vector<std::vector< double> > const& allParticles, std::vector<std::vector< double> > const& lNplNvl, double a, double deltaL, double nT, double timePerIt, double maxValL);

std::vector<double> allAlphaCoefGeo2(std::vector<std::vector< double> > const& allParticles, double lp0, double a, double nT, double nTtminusOne, double h, double deltaL, std::vector<std::vector< double> > const& lNplNvl, double t, double timePerIt, double maxValL);

void advancePdfGeo2(std::vector<double>const& alphaVector, std::vector<std::vector< double> >& allParticles, std::vector<std::vector< double> > & lNplNvl, double h, double nT, double a, double deltaL, double it, double maxValL, double lp0, double nTtminusOne, double timePerIt);

std::vector<std::vector<double> > geo2lNplNv(std::vector<std::vector<double> > allParticles, std::vector<double> liVector, double nT, double lp0, double maxvalL);

#include <stdio.h>

#endif /* aggloGeo2_hpp */
