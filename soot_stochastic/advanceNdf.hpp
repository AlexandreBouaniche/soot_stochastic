//
//  advanceNdf.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 03/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef advanceNdf_hpp
#define advanceNdf_hpp

#include <vector>

std::vector<std::vector<double> > ndf(std::vector<std::vector<double> > allParticles, std::vector<double> liVector, double nT);

std::vector<double> allAlphaCoefNdf(std::vector<std::vector< double> > const& allParticles, double a, double nT, double h, std::vector<std::vector< double> > const& ndft, double timePerIt, std::vector<double> lVector);

void advanceNdf(std::vector<double>const& alphaVector, std::vector<std::vector< double> >& allParticles, std::vector<std::vector< double> > & ndft, double h, double nT, double a,double it, double timePerIt, std::vector<double> lVector);

double totalMassNdf(std::vector<std::vector<double> > ndft);

std::vector<double> alphaCustomNuclAllSizes(std::vector<double> const& dotHvector, std::vector<std::vector< double> > const& allParticles, double a, double nT, double h, std::vector<std::vector< double> > const& ndft, double timePerIt, std::vector<double> lVector);

void advanceNdfNuclAllSizes(std::vector<double> const& dotHvector, std::vector<double>const& alphaVector, std::vector<std::vector< double> >& allParticles, std::vector<std::vector< double> > & ndft, double h, double nT, double a,double it, double timePerIt, std::vector<double> lVector);



#include <stdio.h>

#endif /* advanceNdf_hpp */
