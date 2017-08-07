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

void linerarSurfGrowth(std::vector<std::vector<double> > &allParticles, double deltaM0, double lp0, double maxValL, double deltaL);

void linerarSurfOxi(std::vector<std::vector<double> > &allParticles, double deltaM0, double lp0, double deltaL);

void surfGrowthAging(std::vector<std::vector<double> > &allParticles, double deltaM0, double lp0, double maxValL, double deltaL, double ageFactor);

double outOfBoundSource(std::vector<std::vector<double> >const& allParticles, double nT, double maxValL, double lp0, double deltaL);

void advanceGrowthPdf(std::vector<std::vector<double> >& allParticles, double nT, double maxValL, double lp0, double deltaL, std::vector<std::vector<double> >const& lAndNpL);

void linearGrowth(std::vector<std::vector<double> > &allParticles, double timePerIt);


#include <stdio.h>

#endif /* growth_hpp */
