//
//  Bin.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 28/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "Bin.hpp"
#include "GridCell.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

Bin::Bin(double mConst) : GridCell(mConst), m_np(0.0), m_nv(0.0)
{
    
}

Bin::Bin(double mConst, double mAvg, double mInf, double mSup) : GridCell(mConst, mAvg, mInf, mSup) , m_np(0.0), m_nv(0.0)
{
    
}


Bin::Bin(GridCell gridCell, double np) : GridCell(gridCell.getmConst(), gridCell.getmAvg(), gridCell.getmInf(), gridCell.getmSup()), m_np(np), m_nv(0.0)
{
    
}


Bin::Bin(GridCell gridCell, double np, double nv) : GridCell(gridCell.getmConst(), gridCell.getmAvg(), gridCell.getmInf(), gridCell.getmSup()), m_np(np), m_nv(nv)
{
    
}


double Bin::getNp()
{
    return m_np;
}

double Bin::getNv()
{
    return m_nv;
}

void Bin::setNp(double np)
{
    m_np = np;
}

void Bin::setNv(double nv)
{
    m_nv = nv;
}

void Bin::showBin()
{
    cout << endl;
    cout << "mAvg   np   nv" << endl;
    cout << m_mAvg << "   " << m_np << "   " << m_nv << endl;
}

