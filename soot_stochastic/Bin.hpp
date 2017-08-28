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
    
    Bin(GridCell gridCell, double np);
    Bin(GridCell gridCell, double np, double nv);
    
    
    double getNp();
    double getNv();
    
    void setNp(double np);
    void setNv(double nv);
    
    void showBin();
    
protected:
    
    double m_np;
    double m_nv;
    
};



#endif /* Bin_hpp */
