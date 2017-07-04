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

std::vector<double> liVector(double lp0, double deltaL, double maxValL);

int npLstar(double lStar, std::vector<std::vector<double> > allParticles, double deltaL);

std::vector<std::vector<double> > liNpliNvli(std::vector<std::vector<double> > allParticles, std::vector<double> liVector, double deltaL, double nT);

double nvLstar(double lStar, std::vector<std::vector<double> > allParticles, double nT, double deltaL);

double beta(double l1, double l2);   // should depend on T and Knudsen not on c

double aggloTotSource(std::vector<std::vector<double> > const& allParticles, std::vector<std::vector< double> > const& lNplNvl, double a);

double dotAlStar(double lStar, std::vector<std::vector< double> > const& allParticles, std::vector<std::vector< double> > const& lNplNvl, double a, double deltaL, double nT);

std::vector<double> allAlphaCoef(std::vector<std::vector< double> > const& allParticles, double lp0, double a, double nT, double h, double deltaL, std::vector<std::vector< double> > const& lNplNvl);

void advancePdf(std::vector<double>const& alphaVector, std::vector<std::vector< double> >& allParticles, std::vector<std::vector< double> > & lNplNvl, double h, double nT, double a, double deltaL, double t);


#include <stdio.h>

#endif /* agglomeration_hpp */
