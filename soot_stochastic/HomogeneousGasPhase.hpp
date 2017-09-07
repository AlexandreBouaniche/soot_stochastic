//
//  HomogeneousGasPhase.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 30/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef HomogeneousGasPhase_hpp
#define HomogeneousGasPhase_hpp

#include <stdio.h>
#include <vector>
#include <string>

#include "streams.hpp"

// array with Y values of relevant species as a function of time. known in advance and homogeneous within the "reactor"

// could be replaced by a cantera object like a reactor or flame

class HomogeneousGasPhase
{
public:
    HomogeneousGasPhase();
    HomogeneousGasPhase(std::string pathProject, std::string pathTarget, std::string dataFilename,std::string labelsFilename);
    
    
    // returns the rank in m_t of the closest (superior) value to argument t
    int getRankt(double t);
    
    double getT(double t);
    double getYPAHnucl(double t);
    double getYA1(double t);
    double getYA2(double t);
    double getYA3(double t);
    double getYA4(double t);
    double getYC2H2(double t);
    double getYO2(double t);
    double getYOH(double t);
    
protected:
    
    std::vector<double> m_t;
    std::vector<double> m_T;
    std::vector<double> m_YPAHnucl;
    std::vector<double> m_YA1;
    std::vector<double> m_YA2;
    std::vector<double> m_YA3;
    std::vector<double> m_YA4;
    std::vector<double> m_YC2H2;
    std::vector<double> m_YO2;
    std::vector<double> m_YOH;
    
    // possible development: integrate here a Cantera object reactor/flame...
};

#endif /* HomogeneousGasPhase_hpp */
