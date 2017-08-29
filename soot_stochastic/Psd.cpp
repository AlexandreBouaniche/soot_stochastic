//
//  Psd.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "Psd.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

Psd::Psd(Grid initialGrid)
{
    double initNvT(1.0e10);
    
    std::vector<Bin> psd;
    double npInit = 0.0;
    double nvInit = 0.0;
    
    int i =0;
    for(i=0;i<initialGrid.getSize();i++)
    {
        GridCell gridCell = initialGrid.getGridCell(i);
        Bin bin = Bin(gridCell, npInit, nvInit);
        psd.push_back(bin);
    }
    
    m_binVector = psd;
    m_nvT = initNvT;
}


Psd::Psd(Grid initialGrid, double initNvT)
{
    std::vector<Bin> psd;
    double npInit = 0.0;
    double nvInit = 0.0;
    
    int i =0;
    for(i=0;i<initialGrid.getSize();i++)
    {
        GridCell gridCell = initialGrid.getGridCell(i);
        Bin bin = Bin(gridCell, npInit, nvInit);
        psd.push_back(bin);
    }
    
    m_binVector = psd;
    m_nvT = initNvT;
}



Psd::Psd(Grid initialGrid, double initNvT, std::vector<int> initNpVector)
{
    std::vector<Bin> psd;
    double nvInit = 0.0;
    
    int i =0;
    for(i=0;i<initialGrid.getSize();i++)
    {
        GridCell gridCell = initialGrid.getGridCell(i);
        Bin bin = Bin(gridCell, initNpVector[i], nvInit);
        psd.push_back(bin);
    }
    
    m_binVector = psd;
    m_nvT = initNvT;
}



int Psd::getSize()
{
    int size(0);
    int i(0);
    for(i=0; i<m_binVector.size();i++)
    {
        size++;
    }
    return size;
}


double Psd::getNvT()
{
    return m_nvT;
}


Bin Psd::getBin(int i)
{
    return m_binVector[i];
}


Bin* Psd::getBinPtr(int i)
{
    Bin* ptr = &m_binVector[i];
    return ptr;
}


int Psd::getNpT()
{
    int npT(0);
    
    int i(0);
    for(i=0; i<m_binVector.size();i++)
    {
        npT+= getNp(i);
    }
    
    return npT;
}


double Psd::getmAvg(int i)
{
    Bin* ptr = getBinPtr(i);
    return ptr->getmAvg();
}


double Psd::getmInf(int i)
{
    Bin* ptr = getBinPtr(i);
    return ptr->getmInf();
}

double Psd::getmSup(int i)
{
    Bin* ptr = getBinPtr(i);
    return ptr->getmSup();
}

int Psd::getNp(int i)
{
    Bin* ptr = getBinPtr(i);
    return ptr->getNp();
}

double Psd::getNv(int i)
{
    Bin* ptr = getBinPtr(i);
    return ptr->getNv();
}


double Psd::calculateNd(int i)
{
    double mInf = getmInf(i);
    double mSup = getmSup(i);
    double nv = getNv(i);
    double nd = nv / (mSup - mInf);
    
    return nd;
}


double Psd::calculateMv(int i)
{
    double mAvg = getmAvg(i);
    double nv = getNv(i);
    double mv = mAvg * nv;
    return mv;
}

double Psd::calculateMvT()
{
    double mvT(0);
    int i(0);
    for(i=0; i<m_binVector.size(); i++)
    {
        mvT+= calculateMv(i);
    }
    return mvT;
}


void Psd::setNvt(double nvt)
{
    m_nvT = nvt;
}



void Psd::calcAndSetAllNv()
{
    int npT = getNpT();
    int np(0);
    double nv(0.0);
    
    int i(0);
    for(i=0;i<m_binVector.size();i++)
    {
        np = getNp(i);
        nv = double(np) / double(npT) * m_nvT;
        Bin* ptr = getBinPtr(i);
        ptr->setNv(nv);
    }
}


void Psd::showPsd()
{
    cout << "nvT = " << m_nvT << "   npTot = " << getNpT() << endl;
    cout << "mAvg   np   nv   nd"<<endl;
    int i(0);
    for(i=0; i<m_binVector.size();i++)
    {
        cout << getmAvg(i) << "   " << getNp(i) << "   " << getNv(i) << "   " << calculateNd(i)<< endl;
    }
}

