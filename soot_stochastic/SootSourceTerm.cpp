//
//  SootSourceTerm.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 04/09/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "SootSourceTerm.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

//constant values used for source terms calculation
const double SootSourceTerm::avogadro = 6.024e23;
const double SootSourceTerm::cmass = 12.0 * 1.6605e-24;  // [g]
const double SootSourceTerm::boltzmann = 1.3807e-16;
const double SootSourceTerm::rhoSoot = 1.8; // [g/cm^3]
const double SootSourceTerm::pi = 3.14159265;
const double SootSourceTerm::epsilon = 2.2;
const double SootSourceTerm::epsilonCond = 1.3;
//const double C1 = 1.4558e-3; // conversion to g [ g.m^-1.s^-1.K^(-1/2) ]
//const double C2 = 110.4;   // [K]

SootSourceTerm::SootSourceTerm(HomogeneousGasPhase* gasPhasePtr, Psd* psdPtr) : m_gasPhasePtr(gasPhasePtr), m_psdPtr(psdPtr), m_aggloType("standard"), m_nuclType("standard"), m_hacaType("standard"), m_condType("standard"), m_oxi02Type("standard"), m_oxiOHType("standard"), m_betaConst(0.0)
{
    
}

double SootSourceTerm::partNuclSource(double t)
{
    double dotH(0);
    double nuclMassSource(0);
    double T = m_gasPhasePtr->getT(t);
    double concentrationPAH = m_gasPhasePtr->getYPAHnucl(t); // [mol/cm^3]
    double NPAH = concentrationPAH * avogadro;   // [part/cm^3] cf Nv
    
    if(m_nuclType == "standard")
    {
        double mPAH = 16.0 * cmass;   // pyrene. Only C atoms mass considered as a simplification.
        double mPrimaryParticle = 2.0 * mPAH;
        double beta = betaNucl(mPAH, T);
        
        nuclMassSource = 2 * mPAH * beta * pow(NPAH,2.0);
        
        dotH = nuclMassSource / mPrimaryParticle ;
        
    }
    else
    {
    
        dotH = 0;
        cout << "error m_nuclType" << endl;
    }
    
    return dotH;
}

vector<double> SootSourceTerm::partAggloSourceVector(double t)
{
    vector<double> dotAlVector;
    
    double T = m_gasPhasePtr->getT(t);
    double dotAl(0);
    int j(0);
    
    for(j=0;j<m_psdPtr->getSize();j++)
    {
        dotAl = wmTotj(j, T) / m_psdPtr->getmAvg(j);
        dotAlVector.push_back(dotAl);
    }
    
    return dotAlVector;
}




// test function
void SootSourceTerm::showBoltzmann()
{
    double t = 1.0e-3;
    
    cout << endl;
    cout << endl;
    cout << "kb = " << boltzmann;
    cout << endl;
    cout << 1.0/3.0 << endl;
    cout << "T(t) = " << m_gasPhasePtr->getT(t);
    cout << endl;
    cout << endl;
}


double SootSourceTerm::betaFree(double m1, double m2, double T)
{
    double beta(0);
    double v1 = m1/rhoSoot;    // particle volume in [cm^3]
    double v2 = m2/rhoSoot;
    
    beta = epsilon * pow(  3.0 / (4.0*pi), 1.0/6.0 ) * pow(  6.0*boltzmann*T / rhoSoot  , 1.0/2.0  ) * pow(  1.0/v1 + 1.0/v2  , 1.0/2.0  ) * pow(  pow(v1, 1.0/3.0) + pow(v2, 1.0/3.0)  , 2.0  );
    
    return beta;
}


double SootSourceTerm::betaNucl(double mPAH, double T)
{
    double beta(0);
    double vPAH = mPAH / rhoSoot;  // approx PAH = sphere with density rhoSoot
    
    beta = epsilon * pow(  3.0 / (4.0*pi), 1.0/6.0 ) * pow(  6.0*boltzmann*T / rhoSoot  , 1.0/2.0  ) * pow(  1.0/vPAH + 1.0/vPAH  , 1.0/2.0  ) * pow(  pow(vPAH, 1.0/3.0) + pow(vPAH, 1.0/3.0)  , 2.0  );
    
    return beta;
}


double SootSourceTerm::betaCond(double mPAH, double m1, double T)
{
    double beta(0);
    double vPAH = mPAH / rhoSoot;
    double v1 = m1 / rhoSoot;
    
    beta = epsilonCond * pow(  3.0 / (4.0*pi), 1.0/6.0 ) * pow(  6.0*boltzmann*T / rhoSoot  , 1.0/2.0  ) * pow(  1.0/v1 + 1.0/vPAH  , 1.0/2.0  ) * pow(  pow(v1, 1.0/3.0) + pow(vPAH, 1.0/3.0)  , 2.0  );
    
    return beta;
    
}



// "reaction rate" collision BINk + BINj -> nu_j*BINj + nu_(j+1)*BINj+1
double SootSourceTerm::w_kj(int k, int j, double T)
{
    double wkj(0);
    double mAvgj = m_psdPtr->getmAvg(j);
    double mAvgk = m_psdPtr->getmAvg(k);
    double nvj = m_psdPtr->getNv(j);
    double nvk = m_psdPtr->getNv(k);
    
    double betakj(0.0);
    
    if(m_aggloType == "standard")
    {
        betakj = betaFree(mAvgj, mAvgk, T);
    }
    else
    {
        cout << "problem m_aggloType" << endl;
        betakj = 0.0;
        
    }
    
    
    // !!! Nv not Nd is used for collision rates [part/volume]  !!!
    if(k==j)
    {
        wkj = 0.5 * betakj * nvj * nvk;  // "factor 1/2"
    }
    else
    {
        wkj = betakj * nvj * nvk;
    }
    return wkj;
}


double SootSourceTerm::avgMbinkTojPlusOne(int k, int j)
{
    double mInfk = m_psdPtr->getmInf(k);
    double mSupk = m_psdPtr->getmSup(k);
    
    double avgM_k_jTojplusOne = 2.0/3.0*(pow(mSupk, 3.0)-pow(mInfk, 3.0))/(pow(mSupk, 2.0)-pow(mInfk, 2.0));
    
    return avgM_k_jTojplusOne;
}


double SootSourceTerm::avgMbinjTojPlusOne(int k, int j)
{
    double mSupj = m_psdPtr->getmSup(j);
    
    double avgM_j_jTojplusOne = mSupj -0.5 * avgMbinkTojPlusOne(k, j);
    
    return avgM_j_jTojplusOne;
}



double SootSourceTerm::avgMbinkToj(int k, int j)
{
    double Mkmin = m_psdPtr->getmInf(k);  // Mkmin is mInfk
    double Mkmax = m_psdPtr->getmSup(k);   // Mkmax is mSupk
    double Mjmin = m_psdPtr->getmInf(j);
    double Mjmax = m_psdPtr->getmSup(j);
    double avgM_k_jToj = ((Mjmax-Mjmin)*(pow(Mkmax, 2.0)/2.0-pow(Mkmin, 2.0)/2.0)-(pow(Mkmax, 3.0)/3.0-pow(Mkmin, 3.0)/3.0))/((Mjmax-Mjmin)*(Mkmax-Mkmin)-(pow(Mkmax, 2.0)/2.0-pow(Mkmin, 2.0)/2.0));
    
    return avgM_k_jToj;
}



double SootSourceTerm::aCoef(int k, int j)
{
    double mAvgk = m_psdPtr->getmAvg(k);
    double mInfj = m_psdPtr->getmInf(j);
    double mSupj = m_psdPtr->getmSup(j);
    
    double aCoefficient(1);
    if(k<j)
    {
        aCoefficient = 1 - mAvgk/(mSupj-mInfj);
    }
    else
    {
        cout << "error aCoef";
    }
    
    
    return aCoefficient;
}




double SootSourceTerm::nuj(int k, int j)
{
    
    double mAvgj = m_psdPtr->getmAvg(j);
    
    double nu_j = (aCoef(k, j)*avgMbinkToj(k, j) - (1-aCoef(k, j))*avgMbinjTojPlusOne(k, j))/(mAvgj) + 1;
    
    return nu_j;
}



double SootSourceTerm::nujPlusOne(int k, int j)
{
    
    double nu_jPlusOne(0);
    int jPlusOne = j+1;
    double mAvgjPlusOne = m_psdPtr->getmAvg(j+1);
    
    if(jPlusOne<m_psdPtr->getSize())
    {
        nu_jPlusOne = (1-aCoef(k, j)) * ( avgMbinkTojPlusOne(k, j) + avgMbinjTojPlusOne(k, j) ) / (mAvgjPlusOne);
    }
    else
    {
        nu_jPlusOne = 0;
        cout << "error nujPlusOne j+1 > lVectorsize " << endl;
    }
    
    return nu_jPlusOne;
}



double SootSourceTerm::wmNegRj(int k, int j, double T)
{
    
    double wmNeg_Rj(0);
    double mAvgj = m_psdPtr->getmAvg(j);
    
    if(k<j)
    {
        wmNeg_Rj = -mAvgj* w_kj(k, j, T);
        
    }
    if(k==j)
    {
        wmNeg_Rj = -mAvgj* w_kj(k, j, T);  // for k=j -1 but summing with wmNegRk will give -2. but factor 1/2 for w_kj. compensate
    }
    if(k>j)
    {
        cout << "error k>j " << endl;
        wmNeg_Rj = 0;
    }
    
    return wmNeg_Rj;
}



double SootSourceTerm::wmNegRk(int k, int j, double T)
{
    double wmNeg_Rk(0);
    double mAvgk = m_psdPtr->getmAvg(k);
    
    if(k<j)
    {
        wmNeg_Rk = -mAvgk*w_kj(k, j, T);
        
    }
    if(k==j)
    {
        wmNeg_Rk = -mAvgk*w_kj(k, j, T);  // for k=j -2 but factor 1/2 for w_kj. compensate
    }
    if(k>j)
    {
        cout << "error k>j " << endl;
        wmNeg_Rk = 0;
    }
    
    return wmNeg_Rk;
}


double SootSourceTerm::wmPosRj(int k, int j, double T)
{
    double mAvgj = m_psdPtr->getmAvg(j);
    double wmPosR_j(0);
    if(k<j)
    {
        wmPosR_j = nuj(k, j)* mAvgj * w_kj(k, j, T);
        
    }
    if(k==j)
    {
        wmPosR_j = 0;
    }
    if(k>j)
    {
        cout << "error k>j " << endl;
        wmPosR_j = 0;
    }
    
    return wmPosR_j;
}


double SootSourceTerm::wmPosPlusOneRj(int k, int jPlusOne, double T)
{
    int j = jPlusOne -1;
    
    double wmPosR_jPlusOne(0);
    double mAvgjPlusOne = m_psdPtr->getmAvg(jPlusOne);
    
    if(jPlusOne>0)
    {
        if(k<j)
        {
            
            wmPosR_jPlusOne += nujPlusOne(k, j) * mAvgjPlusOne * w_kj(k, j, T);   // case k < j
            
            
        }
        if(k==j)
        {
            wmPosR_jPlusOne += 1 * mAvgjPlusOne * w_kj(k, j, T);  // case k = j
        }
        if(k>j)
        {
            cout << "error k>j " << endl;
            wmPosR_jPlusOne = 0;
        }
    }
    else
    {
        cout << "error jPlusOne = 0" << endl;
        wmPosR_jPlusOne = 0;
    }
    
    return wmPosR_jPlusOne;
}



double SootSourceTerm::wmNegJ(int j, double T)
{
    double wmNegJ(0);
    int k=0;
    for(k=0; k<j; k++)      // case collision of j with a smaller particle
    {
        wmNegJ += wmNegRj(k, j, T);
    }
    
    wmNegJ += wmNegRj(k, j, T) + wmNegRk(k, j, T); //k=j
    
    for(k=(j+1);(k+1)<m_psdPtr->getSize(); k++) // case collision of j with a bigger particle. only possible if smaller than biggest size  -> for until (k+1)<psd.size
    {
        wmNegJ += wmNegRk(j, k, T);
    }
    
    return wmNegJ;
}



double SootSourceTerm::wmPosJ(int j, double T)
{
    double wmPosj(0);
    int k(0);
    for(k=0; k<j; k++)
    {
        wmPosj += wmPosRj(k, j, T);
    }
    
    //k=j
    wmPosj += wmPosRj(k, j, T); // = 0 in this case
    
    return wmPosj;
}


double SootSourceTerm::wmPosJPlusOne(int jPlusOne, double T)
{
    double wmPosjPlus(0);
    int k(0);
    int j = jPlusOne - 1;
    
    if(jPlusOne>0)
    {
        for(k=0; k<j;k++)
        {
            wmPosjPlus += wmPosPlusOneRj(k, (j+1), T);
        }
        wmPosjPlus += wmPosPlusOneRj(k, (j+1), T);  // k=j
    }
    else
    {
        wmPosjPlus = 0;
    }
    
    return wmPosjPlus;
}


double SootSourceTerm::wmTotj(int j, double T)
{
    double wmTotj(0);
    
    int jPlusOne = j+1;
    if(jPlusOne<m_psdPtr->getSize())
    {
        wmTotj = wmNegJ(j, T) + wmPosJ(j, T) + wmPosJPlusOne(j, T);
    }
    else
    {
        wmTotj = wmPosJPlusOne(j, T);
        
    }
    
    return wmTotj;
}
