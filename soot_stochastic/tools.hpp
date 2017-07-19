//
//  tools.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef tools_hpp
#define tools_hpp

#include <vector>

std::vector<int> randomList(double t, int nbPicked, int maxVal);


std::vector<int> randomListWithoutDuplicate(double t, int nbPicked, int maxVal);

double maxColi(std::vector<std::vector<double> > const& matrix, int col);


int rounding(double d);

bool isIntInlist(int rank, std::vector<int> list);

double frand_a_b(double a, double b);


#include <stdio.h>

#endif /* tools_hpp */
