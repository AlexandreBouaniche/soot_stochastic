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

// array with Y values of relevant species as a function of time. known in advance and homogeneous within the "reactor"

// could be replaced by a cantera object like a reactor or flame

class HomogeneousGasPhase
{
public:
    HomogeneousGasPhase();
    HomogeneousGasPhase(std::vector<double> t);
    
    void setVector(std::vector<double>, std::string label);
    double getValuet(double t, std::string label);
    
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
    
    // Cantera object reactor/flame...
};

#endif /* HomogeneousGasPhase_hpp */
