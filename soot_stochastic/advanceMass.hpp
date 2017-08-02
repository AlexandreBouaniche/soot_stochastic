//
//  advanceMass.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 01/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef advanceMass_hpp
#define advanceMass_hpp

#include <vector>

std::vector<double> allAlphaCoefMass(std::vector<std::vector< double> > const& allParticles, double a, double mT, double h, std::vector<std::vector< double> > const& lNplNvl, double timePerIt, std::vector<double> lVector);

void advancePdfMass(std::vector<double>const& alphaVector, std::vector<std::vector< double> >& allParticles, std::vector<std::vector< double> > & lNplNvl, double h, double mT, double a,double it, double mTminusOne, double timePerIt, std::vector<double> lVector);

double aggloTotMassSource(std::vector<std::vector<double> > const& lNplNvl, double timePerIt, std::vector <double> const& lVector, double a);

double totalMassBins(std::vector<std::vector<double> > const& lNpNv);

std::vector<std::vector<double> > lNplNvMass(std::vector<std::vector<double> > allParticles, std::vector<double> liVector, double mT);

double nTfromMass(std::vector<std::vector<double> > lNpNvMass);

#include <stdio.h>

#endif /* advanceMass_hpp */
