//
//  HomogeneousGasPhase.cpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 30/08/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#include "HomogeneousGasPhase.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

HomogeneousGasPhase::HomogeneousGasPhase()
{
    
}

HomogeneousGasPhase::HomogeneousGasPhase(string pathProject, string pathTarget, string dataFilename,string labelsFilename)
{
    vector<string> labels = readLabels(pathProject, pathTarget , labelsFilename);
    
    vector<vector<double> > dataArray = readDataArray(pathProject, pathTarget, dataFilename, labelsFilename);
    
    int indext = dataIndex(labels, "t");
    int indexT = dataIndex(labels, "T");
    int indexYPAHnucl = dataIndex(labels, "YPAHnucl");
    int indexYA1 = dataIndex(labels, "YA1");
    int indexYA2 = dataIndex(labels, "YA2");
    int indexYA3 = dataIndex(labels, "YA3");
    int indexYA4 = dataIndex(labels, "YA4");
    int indexYC2H2 = dataIndex(labels, "YC2H2");
    int indexYO2 = dataIndex(labels, "YO2");
    int indexYOH = dataIndex(labels, "YOH");
    
    int i(0);
    for(i=0; i<dataArray.size();i++)
    {
        if(indext>0)
            m_t.push_back(dataArray[i][indext]);
        
        if(indexT>0)
            m_T.push_back(dataArray[i][indexT]);
        
        if(indexYPAHnucl>0)
            m_YPAHnucl.push_back(dataArray[i][indexYPAHnucl]);
        
        if(indexYA1>0)
            m_YA1.push_back(dataArray[i][indexYA1]);
        
        if(indexYA2>0)
            m_YA2.push_back(dataArray[i][indexYA2]);
        
        if(indexYA3>0)
            m_YA3.push_back(dataArray[i][indexYA3]);
        
        if(indexYA4>0)
            m_YA4.push_back(dataArray[i][indexYA4]);
        
        if(indexYC2H2>0)
            m_YC2H2.push_back(dataArray[i][indexYC2H2]);
        
        if(indexYO2>0)
            m_YO2.push_back(dataArray[i][indexYO2]);
        
        if(indexYOH>0)
            m_YOH.push_back(dataArray[i][indexYOH]);
        
    }
}



// returns the rank in m_t of the closest (superior) value to argument t

int HomogeneousGasPhase::getRankt(double t)
{
    double tCloser(0);
    int rankCloser(0);
    
    int rank(0);
    
    int size(0);
    int i(0);
    for(i=0; i<m_t.size();i++)
    {
        size++;
    }
    int rankMax = size - 1;
    
    double tmax = m_t[rankMax];
    double tmin = m_t[0];
    
    if(t>tmax)
    {
        cout << endl;
        cout << "t > tmax" << endl;
        tCloser = tmax;
        rankCloser = rankMax;
        
    }
    else if(t<tmin)
    {
        cout << endl;
        cout << "t < tmin" << endl;
        tCloser = tmin;
        rankCloser = 0;
    }
    else
    {
        
        while(m_t[rank] < t)
        {
            rank++;
        }
        
        rankCloser = rank;
    }
    
    return rankCloser;
}



double HomogeneousGasPhase::getT(double t)
{
    int rank = getRankt(t);
    
    
    // for interpolation between two adjascent data lines

    double tSup = m_t[rank];
    double coefSup(1.0);
    double coefInf(0.0);
    double interpolT(0);
    
    if(m_T.size()>rank && m_T.size()>0)
    {
        if(t>=tSup || t<=m_t[0])
        {
            return m_T[rank];
        }
        else  // we go from equations coefSup * m_t[rank] + coefInf * m_t[rank - 1] = t
            // and  coefInf + coefSup = 1
        {
            coefSup = (t - m_t[rank-1]) / (m_t[rank] - m_t[rank - 1]);
            coefInf = 1.0 - coefSup;
            interpolT = coefInf * m_T[rank -1] + coefSup * m_T[rank];
            return interpolT;
        }
    }
    else
    {
        cout << "insufficient information in m_T" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYPAHnucl(double t)
{
    int rank = getRankt(t);
    
    if(m_YPAHnucl.size()>rank && m_YPAHnucl.size()>0)
    {
        return m_YPAHnucl[rank];
    }
    
    else
    {
        cout << "insufficient information in m_YPAHnucl" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYA1(double t)
{
    int rank = getRankt(t);
    
    if(m_YA1.size()>rank && m_YA1.size()>0)
    {
        return m_YA1[rank];
    }
    
    else
    {
        cout << "insufficient information in m_YA1" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYA2(double t)
{
    int rank = getRankt(t);
    
    if(m_YA2.size()>rank && m_YA2.size()>0)
    {
        return m_YA2[rank];
    }
    
    else
    {
        cout << "insufficient information in m_YA2" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYA3(double t)
{
    int rank = getRankt(t);
    
    if(m_YA3.size()> rank  &&  m_YA3.size()>0)
    {
        return m_YA3[rank];
    }
    
    else
    {
        cout << "insufficient information in m_YA3" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYA4(double t)
{
    int rank = getRankt(t);
    
    if(m_YA4.size()>rank  &&  m_YA4.size()>0)
    {
        return m_YA4[rank];
    }
    
    else
    {
        cout << "insufficient information in m_YA4" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYC2H2(double t)
{
    int rank = getRankt(t);
    
    if(m_YC2H2.size()>rank  &&  m_YC2H2.size()>0)
    {
        return m_YC2H2[rank];
    }
    
    else
    {
        cout << "insufficient information in m_YC2H2" << endl;
        return 0.0;
    }
}


double HomogeneousGasPhase::getYO2(double t)
{
    int rank = getRankt(t);
    
    if(m_YO2.size()>rank  &&  m_YO2.size()>0)
    {
        return m_YO2[rank];
    }
    else
    {
        cout << "insufficient information in m_YO2" << endl;
        return 0.0;
    }
    
}

double HomogeneousGasPhase::getYOH(double t)
{
    int rank = getRankt(t);
    
    if(m_YOH.size()>rank  &&  m_YOH.size()>0)
    {
        return m_YOH[rank];
    }
    else
    {
        cout << "insufficient information in m_YOH" << endl;
        return 0.0;
    }
}

