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

void writeFile(std::string pathProject, std::string pathTarget, std::vector<std::vector<double> > allParticles);

void updateFile(std::string pathProject, std::string pathTarget, double t, std::vector<std::vector<double> > allParticles);

void writeCpdf(std::string pathProject, std::string pathTarget, std::vector<std::vector<double> > allParticles, double cSensibility);

void updateCpdf(std::string pathProject, std::string pathTarget, double t, std::vector<std::vector<double> > allParticles, double cSensibility);

void writeCpdft(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > allParticles, double cSensibility);


#include <stdio.h>

#endif /* streams_hpp */
