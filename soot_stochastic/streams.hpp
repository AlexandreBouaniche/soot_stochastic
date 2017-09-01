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


void writePdft(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, std::vector<std::vector<double> > lAndNpl);

void writeNvt(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, double nT, std::vector<std::vector<double> > lAndNpl);

void writeCustomNv(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, double nT);

void writeGeoNdt(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > allParticles, double pdfGrid,double minVal, double maxVal, int column, double nT, std::vector<std::vector<double> > lAndNpl);

void writeGeoNdtDi(std::string pathProject, std::string pathTarget, int it, std::vector<std::vector<double> > lAndNpl);

void writeCustomAggloCase(std::string pathProject, std::string pathTarget, std::vector<double> lVector, double time);

void writeCustomAggloGrowthCase(std::string pathProject, std::string pathTarget, std::vector<double> lVector, double time);

void writeCustomNuclGrowthCase(std::string pathProject, std::string pathTarget, std::vector<double> lVector, double time);



std::vector<double> readDataVectors(std::string pathProject, std::string pathTarget, std::string dataFilename);

std::vector<std::string> readLabels(std::string pathProject, std::string pathTarget ,std::string labelsFilename);

std::vector<std::vector<double> > readDataArray(std::string pathProject, std::string pathTarget, std::string dataFilename,std::string labelsFilename);

int dataIndex(std::vector<std::string> labelsVector, std::string label);

#include <stdio.h>

#endif /* streams_hpp */
