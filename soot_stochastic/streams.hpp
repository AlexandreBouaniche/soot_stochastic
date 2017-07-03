//
//  streams.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 22/06/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef streams_hpp
#define streams_hpp

#include <vector>
#include <string>


void writePdft(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column);


#include <stdio.h>

#endif /* streams_hpp */
