//
//  Grid.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 28/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef Grid_hpp
#define Grid_hpp

#include <stdio.h>
#include <vector>
#include <string>

#include "GridCell.hpp"


// Object that represents the grid in the space of sizes (mass). It is constituted of GridCells

class Grid
{
public:
    Grid();  //  geometric2 without specified parameters
    
    Grid(std::string gridType, int numberTotCells, double firstCellmConst); // for type gemetric2 with specified parameters
    
    Grid(std::string gridType, int numberTotCells, double firstCellmConst, double step); // for regular type
    
    Grid(std::string gridType, int numberTotCells, double minInf, double maxSup, double qCoef);  //for geometricStandard
    
    Grid(std::string gridType, std::vector<double> mConstVector); // For customized
    
    void showGrid();
    //void refineGrid         Idea
    
    GridCell getGridCell(int i);
    
    int getSize();
    
protected:
    std::vector<GridCell> m_gridCellVector;
    std::string m_gridType;     // available options: "regular"; "geometric2"; "geometricStandard"; "custom".
    int m_numberTotCells;
    double m_firstCellmConst;
};


#endif /* Grid_hpp */
