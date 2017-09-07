//
//  SootSourceTerm.hpp
//  soot_stochastic
//
//  Created by Alexandre Bouaniche on 04/09/2017.
//  Copyright © 2017 Alexandre Bouaniche. All rights reserved.
//

#ifndef SootSourceTerm_hpp
#define SootSourceTerm_hpp

#include <stdio.h>

#include <vector>
#include <string>

#include "HomogeneousGasPhase.hpp"
#include "Psd.hpp"

class SootSourceTerm
{
public:
    SootSourceTerm(HomogeneousGasPhase* gasPhasePtr, Psd* psdPtr);
    
    // here we did not put the timePerIt value in the source terms like in previous version (there was *timePerIt in beta for example). Instead, all source terms will be calculated in [s^-1] and will need to be multiplied by timePerIt in S2solver

    // UNITS: [ g, mol or particles (Navogadro linking), cm, s, K ]  !!
    
        // The chosen mass unit is g. the mass of 1 carbon atom Cmass = 12.0 AMU. 1 AMU = 1.66e-24 g. Then 1 Cmass = 1.99e-23 g  and 1 mol of Cmass = 12.00 g
    
        // The particles and PAH diameters are calculated (in cm) from the mass considering spheres with rhoSoot = 1.8 g/cm^3
    
        // gas phase species concentrations are given in [mol/cm^3], particles Nv in [part/cm^3]. concentration * Navogadro = Nv   :    [mol/cm^3]  *  [part/mol] = [part/cm^3]
    
    
    
    
    
    //main public functions (information needed for advancing the solver)
    
    double partNuclSource (double t);      // dotH. PAH pyrene C16H10
    std::vector<double> partAggloSourceVector(double t);  // dotAlVector
    
    double deltamtotGrowth(double t, double particleMass);
    
    //The source terms for nucleation and agglomeration are in particles above   (alphaVector calculated in npSolver because need to have nT(t) and nT(t-dt) )
    // the source terms for surface growth are in mass as they are directly added on each particle.
    
    
    
    // model configuration
    void setAggloType(std::string aggloType);
    void setNuclType(std::string nuclType);
    void setGrowthType(std::string growthType);
    
    void setHacaType(std::string hacaType);
    void setCondType(std::string condType);
    void setOxiO2Type(std::string oxiO2Type);
    void setOxiOHType(std::string oxiOHType);
    
    void setbetaConst(double betaConst);
    
    void showBoltzmann();
    
    // sub-functions of deltamtotGrowth
    double deltamHaca(double t, double particleMass);
    double deltamCond(double t, double particleMass);
    double deltamOxiO2(double t, double particleMass);
    double deltamOxiOH(double t, double particleMass);
    
protected:
    HomogeneousGasPhase* m_gasPhasePtr;
    Psd* m_psdPtr;
    
    std::string m_aggloType; // "standard", "custom", "constantBeta", "fractal"...?
    std::string m_nuclType;  // "standard", "custom", "YA4", "dimers"?...
    std::string m_growthType;  // "standard", "custom"
    
    std::string m_hacaType;
    std::string m_condType;
    std::string m_oxi02Type;
    std::string m_oxiOHType;
    
    double m_betaConst;
    
    //constant values used for source terms calculation
    
    static const double avogadro;
    static const double cmass;
    static const double boltzmann;
    static const double rhoSoot;
    static const double pi;
    static const double epsilon;      //Van der Waals enhancement factor for nucl and agglo
    static const double epsilonCond;  //Van der Waals enhancement factor for condensation
    //static const double C1;           // constants C1 and C2 used for Sutherland's law for continuum regime
    //static const double C2;
    
    
    
    //internal use object subfunctions
    
        // beta: sub-functions of massAggloSource and massNuclSource 1
    
    double betaFree(double m1, double m2, double T);
    //double betaContinuum(double m1, double m2, double T, double P, double sigma);
    //double betaTransition(double m1, double m2, double T, double P, double sigma);
    double betaNucl(double mPAH, double T);
    double betaCond(double mPAH, double m1, double T);
    
        // agglo: sub-functions of massAggloSource and massNuclSource 2
    double w_kj(int k, int j, double T);  // "reaction rate" collision BINk + BINj -> nu_j*BINj + nu_(j+1)*BINj+1
    
    double avgMbinkTojPlusOne(int k, int j);  // used to calculate the nu coefficients for the transfer from one bin to the superior one.
    double avgMbinjTojPlusOne(int k, int j);
    double avgMbinkToj(int k, int j);
    
    double aCoef(int k, int j);
    double nuj(int k, int j);
    double nujPlusOne(int k, int j);
    
    double wmNegRj(int k, int j, double T);
    double wmNegRk(int k, int j, double T);
    double wmPosRj(int k, int j, double T);
    double wmPosPlusOneRj(int k, int jPlusOne, double T);
    double wmNegJ(int j, double T);
    double wmPosJ(int j, double T);
    double wmPosJPlusOne(int j, double T);
    double wmTotj(int j, double T);
    
};

#endif /* SootSourceTerm_hpp */
