//
//  GridCell.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 28/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "GridCell.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

GridCell::GridCell(double mConst) : m_mConst(mConst)
{
    
}

GridCell::GridCell(double mConst, double mAvg, double mInf, double mSup) : m_mConst(mConst), m_mAvg(mAvg), m_mInf(mInf), m_mSup(mSup)
{
    
}

double GridCell::getmConst()
{
    return m_mConst;
}

double GridCell::getmAvg()
{
    return m_mAvg;
}

double GridCell::getmInf()
{
    return m_mInf;
}

double GridCell::getmSup()
{
    return m_mSup;
}

void GridCell::setmConst(double mConst)
{
    m_mConst = mConst;
}

void GridCell::setmAvg(double mAvg)
{
    m_mAvg = mAvg;
}

void GridCell::setmInf(double mInf)
{
    m_mInf = mInf;
}

void GridCell::setmSup(double mSup)
{
    m_mSup = mSup;
}
