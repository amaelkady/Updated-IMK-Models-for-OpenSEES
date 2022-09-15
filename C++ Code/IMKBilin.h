/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

//Modified Ibarra-Medina-Krawinkler with Bilin Hysteretic Response

//**********************************************************************                                                                     
// Code Developed by: Ahmed Elkady
// Lecturer,	University of Southampton
// Last Updated: October 25th 2021
//**********************************************************************

#ifndef IMKBilin_h
#define IMKBilin_h

#include <UniaxialMaterial.h>

class IMKBilin : public UniaxialMaterial
{
public:
    IMKBilin(int tag, double Ke,
        double posUp_0,  double posUpc_0, double posUu_0,  double posFy_0, double posFcapFy_0, double posResF_0,
        double negUp_0,  double negUpc_0, double negUu_0,  double negFy_0, double negFcapFy_0, double negResF_0,
        double LAMBDA_S, double LAMBDA_C, double LAMBDA_K, double c_S, double c_C, double c_K, double D_pos, double D_neg);
    IMKBilin();
    ~IMKBilin();
    const char *getClassType(void) const { return "IMKBilin"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double  getStrain(void);
    double  getStress(void);
    double  getTangent(void);
    double  getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);


protected:

private:
    //my functions

    //Fixed input material parameters 
    double  Ke;
    double  posUp_0;
    double  posUpc_0;
    double  posUu_0;
    double  posFy_0;
    double  posFcapFy_0;
    double  posResF_0;
    double  negUp_0;
    double  negUpc_0;
    double  negUu_0;
    double  negFy_0;
    double  negFcapFy_0;
    double  negResF_0;
    double  LAMBDA_S;
    double  LAMBDA_C;
    double  LAMBDA_K;
    double  c_S;
    double  c_C;
    double  c_K;
    double  D_pos;
    double  D_neg;

    //State variables 
    double  U,              cU;

    //History variables
    // U and F
    double  posUy_0,        negUy_0;
    double  posUcap_0,      negUcap_0;
    double  posFcap_0,      negFcap_0;
    double  posKp_0,        negKp_0;
    double  posKpc_0,       negKpc_0;
    double  posFr_0,        cPosFr_0;
    double  negFr_0,        cNegFr_0;
    double  posFyProj_0,    negFyProj_0;
    double  posFcapProj_0,  negFcapProj_0;

    double  posUy_1,        cPosUy_1;
    double  posUcap_1,      cPosUcap_1;
    double  posFy_1,        cPosMy_1;
    double  posFcap_1,      cPosMmax_1;
    double  posFyProj_1,    cPosMyProj_1;
    double  posFcapProj_1,  cPosMmaxProj_1;
    double  posKp_1,        cPosKp_1;
    double  posKpc_1,       cPosKpc_1;

    double  negUy_1,		cNegUy_1;
    double  negUcap_1,		cNegUcap_1;
    double  negFy_1,		cNegMy_1;
    double  negFcap_1,		cNegMmax_1;
    double  negFyProj_1,	cNegMyProj_1;
    double  negFcapProj_1,	cNegMmaxProj_1;
    double  negKp_1,		cNegKp_1;
    double  negKpc_1,		cNegKpc_1;
    // Basic
    double  Ui,				cUi;
    double  Fi,				cFi;
    double  Ui_1,			cUi_1;
    double  Fi_1,			cFi_1;
    // int     Di,				cDi;
    int     Di_1,			cDi_1;
    // Energy

    // double  betaS,		cBetaS;
    // double  betaC,		cBetaC;
    // double  betaK,		cBetaK;

    double  refEnergyS;
    double  refEnergyC;
    double  refEnergyK;

    double  engExcr_1,		cEngExcr_1;
    double  engExcr,		cEngExcr;
    double  engRvrs,		cEngRvrs;
    double  engTotl,		cEngTotl;
    // Stiffness
    double  TangentK,       cTangentK;
    double  K_j,          cK_j;
    // Flag
    bool    Excursion_Flag, cExcursion_Flag;
    bool    Reversal_Flag,  cReversal_Flag;
    bool    Yield_Flag,     cYield_Flag;
    bool    Fail_FlagPos,   cFail_FlagPos;
    bool    Fail_FlagNeg,   cFail_FlagNeg;
    bool    Mrpos_Flag,     cMrpos_Flag;
    bool    Mrneg_Flag,     cMrneg_Flag;
    bool    Energy_Flag,    cEnergy_Flag;
    // Local
    double  Ulocal,         cUlocal;
    double  Flocal,	        cFlocal;
    // Other

};

#endif

//////////////////////////// Variables Definitions /////////////////////////////
/*
Ke 					Initial elastic stiffness
posUp_0 		Initial pre-capping plastic rotation in the +ve loading direction
posUpc_0   	Initial post-capping plastic rotation in the +ve loading direction
posUu_0 		Ultimate rotation in the +ve loading direction
posFy_0 			Initial effective plastic moment in the +ve loading direction
posFcapFy_0 		Initial maximum-to-effective plastic moment ratio in the +ve loading direction
posResF_0 			Residual moment in the +ve loading direction
negUp_0 		Initial pre-capping plastic rotation in the -ve loading direction
negUpc_0   	Initial post-capping plastic rotation in the -ve loading direction
negUu_0    	Ultimate rotation in the -ve loading direction
negFy_0        	Initial effective plastic moment in the -ve loading direction
negFcapFy_0    	Initial maximum-to-effective plastic moment ratio in the -ve loading direction
negResF_0       	Residual moment in the -ve loading direction
LAMBDA_S 			Cyclic deterioration parameter for strength deterioration 
LAMBDA_C 			Cyclic deterioration parameter for post-capping strength deterioration
LAMBDA_K 			Cyclic deterioration parameter for unloading stiffness deterioration 
c_S 				Rate of strength deterioration.
c_C 				Rate of post-capping strength deterioration.
c_K 				Rate of unloading stiffness deterioration.
D_pos 				Rate of cyclic deterioration in the +ve loading direction
D_neg 				Rate of cyclic deterioration in the -ve loading direction
n					Parameter identifying the offset rotation on the unloading side
Roffset				Offset rotation identifying the Smooth Transition region
LAMBDA_F			Cyclic deterioration parameter for Smooth Transition deterioration 
c_F					Cyclic deterioration parameter for Smooth Transition deterioration
Ui 					Rotation at current step
Fi 					Moment at current step 
Di 					Rotation Direction at current step 
Ui_1  				Rotation at previous step
Fi_1            	Moment at previous step 
Di_1            	Rotation Direction at previous step
Ulocal 			Rotation at direction reversal points
Flocal 			Moment at direction reversal points
TangentK 			Tangent stiffness
K_j 				Unloading stiffness in the previous excursion
posUy_1 	Yielding rotation in the previous +ve excursion
posUcap_1 	Capping point rotation in the previous +ve excursion
posKp_1 	Pre-capping slope in the previous +ve excursion
posKpc_1 	Post-capping slope in the previous +ve excursion
posFy_1 		Effective plastic moment in the previous +ve excursion
posFyProj_1  Projed effective plastic moment in the previous +ve excursion
posFcap_1        Maximum moment in the previous +ve excursion
posFcapProj_1 Projed maximum  moment in the previous +ve excursion
negUy_1     Yielding rotation in the previous -ve excursion
negUcap_1   Capping point rotation in the previous -ve excursion
negKp_1     Pre-capping slope in the previous -ve excursion
negKpc_1    Post-capping slope in the previous -ve excursion
negFy_1         Effective plastic moment in the previous -ve excursion
negFyProj_1  Projed effective plastic moment in the previous -ve 
negFcap_1        Maximum moment in the previous -ve excursion
negFcapProj_1 Projed maximum  moment in the previous -ve excursion
betaS    
betaC    
betaK    
betaF 	  
refEnergyS    	Refernence energy for strength deterioration
refEnergyC    	Refernence energy for post-capping strength deterioration
refEnergyK    	Refernence energy for unloading stiffness deterioration
refEnergyF    	Refernence energy for Smooth Transition deterioration
Excursion_Flag 		Flag for Excursion occurrence (i.e.,	crossing the x-axis)
Reversal_Flag 		Flag for Loading direction reversal occurrence
Yield_Flag 			Flag for Yielding occurrence
Fail_FlagPos 		Flag for reaching the ultimate rotation in the +ve loading direction
Fail_FlagNeg 		Flag for reaching the ultimate rotation in the -ve loading direction
Mrpos_Flag 			Flag for reaching the residual moment in the +ve loading direction
Mrneg_Flag 			Flag for reaching the residual moment in the -ve loading direction
Energy_Flag 		Flag for reaching the reference energy
engExcr_1	Dissipated energy in previous excursion
engExcr 		Dissipated energy in current excursion
engRvrs 			Total dissipated energy till previous load reversal point
engTotl 		Total dissipated energy till current step
*/
