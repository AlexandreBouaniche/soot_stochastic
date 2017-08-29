//
//  Bin.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 28/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef Bin_hpp
#define Bin_hpp

#include <stdio.h>
#include <vector>

#include "GridCell.hpp"

class Bin : public GridCell
{
public:
    Bin(double mConst);
    Bin(double mConst, double mAvg, double mInf, double mSup);
    
    Bin(GridCell gridCell, int np);
    Bin(GridCell gridCell, int np, double nv);
    
    
    int getNp();
    double getNv();
    
    void setNp(int np);
    void setNv(double nv);
    
    void showBin();
    
    //TO IMPLEMENT  void countAndSetNp. // counts stochastic particles of AerosolPhase which are within the Bin interval. Then setNp
    
protected:
    
    int m_np;  // [stochastic particles]
    double m_nv;  // [real particles / m3]
    
};



#endif /* Bin_hpp */
