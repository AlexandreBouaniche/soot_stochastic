//
//  GridCell.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 28/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef GridCell_hpp
#define GridCell_hpp

#include <stdio.h>
#include <vector>


// Object that represents one "empty" cell of the grid in the space of sizes (mass).
// different attributes representing the mass of the cell "by construction" the corresponding
// average mass of the cell and the inf and sup born masses

class GridCell
{
public:
    
    //Creator
    GridCell(double mConst);
    GridCell(double mConst, double mAvg, double mInf, double mSup);
    
    
    double getmConst();
    double getmAvg();
    double getmInf();
    double getmSup();
    
    void setmConst(double mConst);
    void setmAvg(double mAvg);
    void setmInf(double mInf);
    void setmSup(double mSup);
    
    
protected:
    double m_mConst;  // mass "by construction" of the grid. For example corresponds to a geometric grid series pattern
    double m_mAvg;
    double m_mInf;
    double m_mSup;
    
};

#endif /* GridCell_hpp */
