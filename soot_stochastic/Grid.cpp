//
//  Grid.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 28/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "Grid.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

// by default a grid of type geometric2 is built, with numberTotCells = 30 and firstCellmConst = 1.0

Grid::Grid()
{
    m_gridType = "geometric2";
    m_numberTotCells = 30;
    m_firstCellmConst = 1.0;
    
    std::vector<GridCell> gridCellVector;
    
    int i = 0;
    for(i=0; i<m_numberTotCells; i++)
    {
        double power = double(i);
        double mConst = m_firstCellmConst * pow(2.0, power);
        double mAvg = 1.125 * mConst;
        double mInf = 0.75 * mConst;
        double mSup = 1.5 * mConst;
        GridCell gridCell = GridCell(mConst, mAvg, mInf, mSup);
        gridCellVector.push_back(gridCell);
    }
    
    m_gridCellVector = gridCellVector;
}


Grid::Grid(std::string gridType, int numberTotCells, double firstCellmConst) : m_gridType(gridType), m_numberTotCells(numberTotCells), m_firstCellmConst(firstCellmConst)
{
    if(m_gridType == "geometric2")
    {
        std::vector<GridCell> gridCellVector;
        
        int i = 0;
        for(i=0; i<m_numberTotCells; i++)
        {
            double power = double(i);
            double mConst = m_firstCellmConst * pow(2.0, power);
            double mAvg = 1.125 * mConst;
            double mInf = 0.75 * mConst;
            double mSup = 1.5 * mConst;
            GridCell gridCell = GridCell(mConst, mAvg, mInf, mSup);
            gridCellVector.push_back(gridCell);
        }
        
        m_gridCellVector = gridCellVector;
    }
    else
    {
        cout << "error Grid constructor gridType" << endl;
        
    }
    
}



Grid::Grid(std::string gridType, int numberTotCells, double firstCellmConst, double step) : m_gridType(gridType), m_numberTotCells(numberTotCells), m_firstCellmConst(firstCellmConst)
{
    if(gridType == "regular")
    {
        std::vector<GridCell> gridCellVector;
        
        int i = 0;
        for(i=0; i<m_numberTotCells; i++)
        {
            double mConst = m_firstCellmConst + double(i*step);
            double mAvg = mConst;
            double mInf = mConst - step/2.0;
            double mSup = mConst + step/2.0;
            GridCell gridCell = GridCell(mConst, mAvg, mInf, mSup);
            gridCellVector.push_back(gridCell);
        }
        
        m_gridCellVector = gridCellVector;
    }
    else
    {
        cout << "error Grid constructor gridType" << endl;
    }
}



Grid::Grid(std::string gridType, int numberTotCells, double minInf, double maxSup, double qCoef) : m_gridType(gridType), m_numberTotCells(numberTotCells)
{
    if(gridType == "geometricStandard")
    {
        std::vector<GridCell> gridCellVector;
        
        double mInf = minInf;
        
        // expression from thesis Molenaar p33. expression for L(i+1/2) which corresponds to mSup(i)
        int i=0;
        for(i=0; i<m_numberTotCells; i++)
        {
            int N = m_numberTotCells - 1;
            double power = double (i) - double(N);
            double mSup = minInf + pow(2.0,double(power/qCoef)) * (maxSup - minInf);
            double mAvg = mInf + (mSup - mInf)/2.0;
            double mConst = mAvg;
            
            GridCell gridCell = GridCell(mConst, mAvg, mInf, mSup);
            gridCellVector.push_back(gridCell);
            
            mInf = mSup;
        }
        
        m_gridCellVector = gridCellVector;
    }
    else
    {
        cout << "error Grid constructor gridType" << endl;
    }
}


//input vector with the first mInf, all the mConst and the last mSup

Grid::Grid(std::string gridType, std::vector<double> mConstVector) : m_gridType(gridType)
{
    if(m_gridType == "custom")
    {
        std::vector<GridCell> gridCellVector;
        
        double mInf = mConstVector[0];
        double mSup;
        
        int i=1;
        for(i=1; i<(mConstVector.size()-1); i++)
        {
            double mConst = mConstVector[i];
            
            if(i>1)
            {
                mInf = mConstVector[i-1] + (mConstVector[i] - mConstVector[i-1]) / 2.0;
            }
            else
            {
                mInf = mConstVector[0];
                
            }
            
            if(i<(mConstVector.size()-2))
            {
                mSup = mConstVector[i] + (mConstVector[i+1] - mConstVector[i]) /2.0;
            }
            else
            {
                mSup = mConstVector[i+1];
            }
            
            double mAvg = mInf + (mSup - mInf) / 2.0;
            
            GridCell gridCell = GridCell(mConst, mAvg, mInf, mSup);
            gridCellVector.push_back(gridCell);
            
        }
        
        m_gridCellVector = gridCellVector;
    }
    else
    {
        cout << "error Grid constructor gridType" << endl;
    }
}



void Grid::showGrid()
{
    int i=0;
    
    cout << endl << endl;
    cout << "mConst   mAvg   mInf   mSup" << endl;
    
    for(i=0; i<m_gridCellVector.size();i++)
    {
        GridCell gridCell = m_gridCellVector[i];
        cout << gridCell.getmConst() << "   " << gridCell.getmAvg() << "   " << gridCell.getmInf() << "   " << gridCell.getmSup() << endl;
    }
}
