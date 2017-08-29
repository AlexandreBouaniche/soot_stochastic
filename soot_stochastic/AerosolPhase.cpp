//
//  AerosolPhase.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 29/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "AerosolPhase.hpp"
#include "GasParticle.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

AerosolPhase::AerosolPhase()
{
    vector<GasParticle> aerosolPhase;
    
    int npTot = 10000;
    
    double m(1.0);
    double c0(0.5);
    
    int i(0);
    for(i=0; i<npTot;i++)
    {
        m = m*1.001;
        GasParticle gasPart(m,c0);
        aerosolPhase.push_back(gasPart);
    }
    
    m_partVector = aerosolPhase;
}


AerosolPhase::AerosolPhase(int npTot)
{
    vector<GasParticle> aerosolPhase;
    
    double m(1.0);
    double c0(0.5);
    
    int i(0);
    for(i=0; i<npTot;i++)
    {
        m = m*1.001;
        GasParticle gasPart(m,c0);
        aerosolPhase.push_back(gasPart);
    }
    
    m_partVector = aerosolPhase;
}


AerosolPhase::AerosolPhase(int npTot, std::vector<double> massVector)
{
    vector<GasParticle> aerosolPhase;
    
    double c0(0.5);
    double m(0);
    
    int i(0);
    for(i=0; i<npTot;i++)
    {
        m = massVector[i];
        GasParticle gasPart(m,c0);
        aerosolPhase.push_back(gasPart);
    }
    
    m_partVector = aerosolPhase;
}


AerosolPhase::AerosolPhase(int npTot, std::vector<double> massVector, std::vector<double> cVector)
{
    vector<GasParticle> aerosolPhase;
    
    double c(0.0);
    double m(0.0);
    
    int i(0);
    for(i=0; i<npTot;i++)
    {
        m = massVector[i];
        c = cVector[i];
        GasParticle gasPart(m,c);
        aerosolPhase.push_back(gasPart);
    }
    
    m_partVector = aerosolPhase;
}


int AerosolPhase::getNpTot()
{
    int npTot(0);
    int i(0);
    for(i=0; i<m_partVector.size();i++)
    {
        npTot++;
    }
    
    return npTot;
}


GasParticle* AerosolPhase::getGasPartPtr(int i)
{
    GasParticle* ptr = &m_partVector[i];
    return ptr;
}


double AerosolPhase::getCi(int i)
{
    GasParticle* ptr = &m_partVector[i];
    double c = ptr->getC();
    return c;
}


void AerosolPhase::setCi(int i, double c)
{
    GasParticle* ptr = &m_partVector[i];
    ptr->setC(c);
}

double AerosolPhase::getMass(int i)
{
    GasParticle* ptr = &m_partVector[i];
    double mass = ptr->getMass();
    return mass;
}

void AerosolPhase::setMass(int i, double mass)
{
    GasParticle* ptr = &m_partVector[i];
    ptr->setMass(mass);
}

void AerosolPhase::addMass(int i, double deltaMass)
{
    GasParticle* ptr = &m_partVector[i];
    ptr->addMass(deltaMass);
}


int AerosolPhase::countNpl(double mInf, double mSup)
{
    double npl(0);
    int i(0);
    for(i=0; i<m_partVector.size();i++)
    {
        double m = getMass(i);
        if(m >= mInf  &  m < mSup)
        {
            npl++;
        }
    }
    return npl;
}
