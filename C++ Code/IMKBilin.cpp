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
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

#include <math.h>
#include <IMKBilin.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <algorithm>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

static int numIMKBilinMaterials	= 0;

void *
OPS_IMKBilin(void)
{
    if (numIMKBilinMaterials == 0) {
        numIMKBilinMaterials++;
        OPS_Error("Mod. IMK Bilinear Model - AE-Aug22\n", 1);
    }

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial	= 0;

    int    iData[1];
    double	dData[21];
    int numData	= 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial IMKBilin tag" << endln;
        return 0;
    }

    numData	= 21;

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid Args want: uniaxialMaterial IMKBilin tag? Ke? ";
        opserr << "Up_pos? Upc_pos? Uu_pos? My_pos? MmaxMy_pos? ResM_pos? ";
        opserr << "Up_neg? Upc_neg? Uu_neg? My_neg? MmaxMy_neg? ResM_neg? ";
        opserr << "LamdaS?  LamdaC? LamdaK? Cs? Cc? Ck? D_pos? D_neg? ";
        return 0;
    }


    // Parsing was successful, allocate the material
    theMaterial	= new IMKBilin(iData[0],
        dData[0],
        dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
        dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type IMKBilin Material\n";
        return 0;
    }

    return theMaterial;
}

IMKBilin::IMKBilin(int tag, double	p_Ke,
    double	p_posUp_0, double	p_posUpc_0, double	p_posUu_0, double	p_posFy_0, double	p_posFcapFy_0, double	p_posResF_0,
    double	p_negUp_0, double	p_negUpc_0, double	p_negUu_0, double	p_negFy_0, double	p_negFcapFy_0, double	p_negResF_0,
    double	p_LAMBDA_S, double	p_LAMBDA_C, double	p_LAMBDA_K, double	p_c_S, double	p_c_C, double	p_c_K, double	p_D_pos, double	p_D_neg)
    :UniaxialMaterial(tag, 0), Ke(p_Ke),
    posUp_0(p_posUp_0), posUpc_0(p_posUpc_0), posUu_0(p_posUu_0), posFy_0(p_posFy_0), posFcapFy_0(p_posFcapFy_0), posResF_0(p_posResF_0),
    negUp_0(p_negUp_0), negUpc_0(p_negUpc_0), negUu_0(p_negUu_0), negFy_0(p_negFy_0), negFcapFy_0(p_negFcapFy_0), negResF_0(p_negResF_0),
    LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
    this->revertToStart();
}

IMKBilin::IMKBilin()
    :UniaxialMaterial(0, 0), Ke(0),
    posUp_0(0), posUpc_0(0), posUu_0(0), posFy_0(0), posFcapFy_0(0), posResF_0(0),
    negUp_0(0), negUpc_0(0), negUu_0(0), negFy_0(0), negFcapFy_0(0), negResF_0(0),
    LAMBDA_S(0), LAMBDA_C(0), LAMBDA_K(0), c_S(0), c_C(0), c_K(0), D_pos(0), D_neg(0)
{
    this->revertToStart();
}

IMKBilin::~IMKBilin()
{
    // does nothing
}

int IMKBilin::setTrialStrain(double	strain, double	strainRate)
{
    //all variables to the last commit
    this->revertToLastCommit();

    //state determination algorithm: defines the current force and tangent stiffness
    U       = strain; //set trial displacement
    Ui_1    = Ui;
    Fi_1    = Fi;
    // Di_1    = Di;
    Ui      = U;

    double dU   = Ui - Ui_1;
    double dEi;
    if (dU == 0) {
        Fi  = Fi_1;
        dEi = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%% INITIALIZE CURRENT BACKBONE VALUES AS PREVIOUS %%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    double	posFcap     = posFcap_1;
    double	posFy	    = posFy_1;
    double	posFyProj	= posFyProj_1;
    double	posFcapProj = posFcapProj_1;
    double	posKp       = posKp_1;
    double	posKpc	    = posKpc_1;
    double	posUy       = posUy_1;
    double	posUcap     = posUcap_1;

    double	negFcap     = negFcap_1;
    double	negFy	    = negFy_1;
    double	negFyProj	= negFyProj_1;
    double	negFcapProj = negFcapProj_1;
    double	negKp       = negKp_1;
    double	negKpc	    = negKpc_1;
    double	negUy       = negUy_1;
    double	negUcap     = negUcap_1;

    double	Mi_boundary = 0.0;

    int	    QuarterFlag = 0, Di;
    double  betaS=0, betaC=0, betaK=0;
    double  Rintrsct_K, DISP_Rev;
    double	Ki, Kpi, Kpci, Ucapi, FyProji, FcapProji;

    double	Mi_temp;

    Reversal_Flag   = false;

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ////////////////// BRANCH DETERMINATION AND FLAG RAISE ////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
    //  Branch
    //      0:  Elastic
    //      1:  Unloading Branch

    //      5:  Towards Capping Point   +
    //      6:  Towards Residual Point  +
    //      7:  Residual Branch         +

    //      15: Towards Capping Point   -
    //      16: Towards Residual Point  -
    //      17: Residual Branch         -

    //  Flag
    //      Yield_Flag:     Preserved.      When the deformation exceeds yield capacity for the first time.
    //      Excursion_Flag: Not preserved.  When crossing X-axis. Evokes re-considering of the deteriorations and which peak to go for.
    //      Reversal_Flag:  Not preserved.  When unloading starts. Evokes re-condiersing of the stiffness deterioration and peak point registration.

    // Find the direction of the current increment "Di": 1:if Ui is moving right    and -1: if Ui is moving left
    // 方向を判定する
    if (Ui > Ui_1) {
        Di	= 1;
    }
    else {
        Di	= -1;
    }

    //  Simple Notation for current step parameters
    Fi	= Fi_1 + K_j * dU;

    //  Get Information before first Yield
    // CHECK FOR YIELDING
    //  Check if previous point was a reversal point
    // CHECK FOR REVERSAL
    if (Di_1 * Di < 0) {
        Reversal_Flag   = true;
        Ulocal	= Ui_1;
        Flocal	= Fi_1;
    }

    // Update loading / unloading stiffness at load reversals
    // UPDATE PEAK POINTS and DETERIORATIONS PARAMETERS
    if (Reversal_Flag) {
        Rintrsct_K  = Ulocal - Flocal / K_j;      // 除荷終了点
        DISP_Rev    = engTotl - engExcr_1 - 0.5*Flocal *(Rintrsct_K - Ulocal);  // 割引吸収エネルギー
        betaK       = pow((DISP_Rev / (2 * refEnergyK - engTotl + 0.5*Flocal * (Rintrsct_K - Ulocal))), c_K);

        K_j         = K_j * (1 - betaK);
        if (Mrpos_Flag || Mrneg_Flag) {
            K_j     = 0.5*Ke;
        }
    }

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        /////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Update Positive Backbone and Target Peak Point

    //  Calculate Backbone parameters at current excursion based on Energy Dissipated in the previous Excursion
    if (Excursion_Flag) {
        betaS   = pow((engExcr / (refEnergyS - engTotl)), c_S);
        betaC   = pow((engExcr / (refEnergyC - engTotl)), c_C);

        if (Ui > Ui_1) {
            //  Update My, Mmax Projion, Kp, and Kpc for Current Step
            posFy       = posFy_1       * (1.0 - betaS * D_pos);
            posKp       = posKp_1       * (1.0 - betaS * D_pos);
            posFcapProj = posFcapProj_1 * (1.0 - betaC * D_pos);
            if (posFr_0 == 0.0) {
                posKpc  = posKpc_0 * (posFcapProj - posFr_0) / posFcapProj;
            }
            else {
                posKpc  = posKpc_0 * (posFy - posFr_0) / (posFy_0 - posFr_0);
            }

            //  Calculate Rotation at Capping Point (Intersection of the two Slopes)
            posUy           = posFy / K_j;
            posFyProj       = posFy - posKp * posUy;
            posUcap         = fabs((posFcapProj - posFyProj) / (posKpc + posKp));
            posFcap         = posFyProj + posUcap * posKp;

            if ((posFcap - posFr_0) / (posUcap + fabs(Ui)- posFr_0/K_j) < posKp) {
                posKp       = (posFcap - posFr_0) / (posUcap + fabs(Ui) - posFr_0 / K_j);
                posFyProj   = posFy - posKp * posUy;
                posUcap     = fabs((posFcapProj - posFyProj) / (posKpc + posKp));
                posFcap     = posFyProj + posUcap * posKp;
            }
        }
        else {
            //  Update My, Mmax Projion, Kp, and Kpc for Current Step
            negFy	    = negFy_1 * (1.0 - betaS * D_neg);
            negFcapProj	= negFcapProj_1 * (1.0 - betaC * D_neg);
            negKp	= negKp_1 * (1.0 - betaS * D_neg);
            if (negFr_0 == 0.0) {
                negKpc	= negKpc_0 * (negFcapProj - negFr_0) / negFcapProj;
            }
            else {
                negKpc	= negKpc_0 * (negFy - negFr_0) / (negFy_0 - negFr_0);
            }

            //  Calculate Rotation at Capping Point (Intersection of the two Slopes)
            negUy       = negFy / K_j;
            negFyProj  = negFy - negKp * negUy;
            negUcap     = fabs((negFcapProj - negFyProj) / (negKpc + negKp));
            negFcap     = negFyProj + negUcap * negKp;

            if ((negFcap - negFr_0) / (negUcap + fabs(Ui) - negFr_0 / K_j) < negKp) {
                negKp       = (negFcap - negFr_0) / (negUcap + fabs(Ui) - negFr_0 / K_j);
                negFyProj  = negFy - negKp * negUy;
                negUcap     = fabs((negFcapProj - negFyProj) / (negKpc + negKp));
                negFcap     = negFyProj + negUcap * negKp;
            }
        }
    }
    // If the residual moment is reached in a given direction, Override the values of Mmax, Ucap, Kp and Kpc
    // Residual Curveに入れる
    if (Di == 1) {
        if (posFcap < posFr_0) {
            posFcap = posFr_0;
            posKpc  = pow(10., -6);
            posKp   = pow(10., -6);
            posUcap = pow(10., -6);
        }
        posUy_1         = posUy;
        posUcap_1       = posUcap;
        posKp_1         = posKp;
        posKpc_1        = posKpc;
        posFy_1         = posFy;
        posFyProj_1     = posFyProj;
        posFcap_1       = posFcap;
        posFcapProj_1   = posFcapProj;
    }
    else {
        if (negFcap < negFr_0) {
            negFcap = negFr_0;
            negKpc  = pow(10., -6);
            negKp   = pow(10., -6);
            negUcap = pow(10., -6);
        }
        negUy_1         = negUy;
        negUcap_1       = negUcap;
        negKp_1         = negKp;
        negKpc_1        = negKpc;
        negFy_1         = negFy;
        negFyProj_1     = negFyProj;
        negFcap_1       = negFcap;
        negFcapProj_1   = negFcapProj;
    }

    // %%%%%%%%%% PREPARE RETURN VALUES %%%%%%%%%%%%%
    //  Simple and unified notation for current bacbone parameters
    // x軸を越えたら更新する
    if (Fi_1 + K_j * dU > 0.0) {
        Ki          = K_j;
        Kpi         = posKp;
        Kpci        = posKpc;
        Ucapi       = posUcap;
        FyProji     = posFyProj;
        FcapProji   = posFcapProj;
    }
    else {
        Ki          = K_j;
        Kpi         = negKp;
        Kpci        = negKpc;
        Ucapi       = negUcap;
        FyProji     = negFyProj;
        FcapProji   = negFcapProj;
    }

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////// COMPUTE FORCE INCREMENT /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////

    //  Moment Calculation Based on unloading/reloading stiffeness
    Fi              = Fi_1 + Ki * dU;
// Capを越えたかどうかの判定

    //  Location Flags
    if (Ui > 0){
        if (Fi > 0){
            QuarterFlag = 1;
        }
        else{
            QuarterFlag = 2;
        }
    }
    else{
        if (Fi > 0){
            QuarterFlag = 4;
        }
        else{
            QuarterFlag = 3;
        }
    }

    // if (Ui > 0){
    //     if (Ui < Ucapi){
    //         Mi_boundary = FyProji + Kpi*Ui;
    //     }else{
    //         Mi_boundary = FcapProji - Kpci*Ui;
    //     }
    //     if (Mi_boundary < posFr_0){
    //         Mi_boundary = posFr_0;
    //         Mrpos_Flag  = 1;
    //     }
    // }else{
    //     if (Ui > -Ucapi){
    //         Mi_boundary = -FyProji + Kpi*Ui;
    //     }else{
    //         Mi_boundary = -FcapProji - Kpci*Ui;
    //     }
    //     if (Mi_boundary > -negFr_0){
    //         Mi_boundary = -negFr_0;
    //         Mrneg_Flag  = 1;
    //     }
    // }


    // Get Boundary Moment at Current Step Based on Current BackBone Curve
    if (QuarterFlag == 1) {         // Ui > 0
        if (fabs(Ui) < Ucapi) {
            Mi_boundary =               FyProji     + Kpi * Ui;
        }
        else{
            Mi_boundary = max(posFr_0,  FcapProji   - Kpci * Ui);
        }
        if (Mi_boundary < posFr_0) {
            Mrpos_Flag  = true;
        }
    }
    else if (QuarterFlag == 3) {    // Ui < 0
        if (fabs(Ui) < Ucapi) {
            Mi_boundary = -FyProji + Kpi * Ui;
            //cout << "        FyProji=" << FyProji << " Kpi=" << Kpi << " Mbound=" << Mi_boundary << endln;

        }
        else if (fabs(Ui) > Ucapi) {
            Mi_boundary = min(-negFr_0, -FcapProji - Kpci * Ui);
        }
        if (Mi_boundary > -negFr_0) {
            Mrneg_Flag  = true;
        }
    }
    else if (QuarterFlag == 2) {
        Mi_boundary     = min(-negFr_0, -FyProji + Kpi * fabs(Ui));
        if (Mi_boundary == -negFr_0 && TangentK==1.e-6) {
            Mrneg_Flag  = true;
        }
    }
    else if (QuarterFlag == 4) {
        Mi_boundary     = max(posFr_0, FyProji - Kpi * fabs(Ui));
        if (Mi_boundary == posFr_0 && TangentK == 1.e-6) {
            Mrneg_Flag  = true;
        }
    }

    //cout << "                Fi_1=" << Fi_1 << " Fi=" << Fi << " TangentK=" << TangentK << " Mbound=" << Mi_boundary << " Q=" << QuarterFlag << endln;

    //  Check for Fail Flag
    // 限界変形を超えた時
    if ((Ui > posUu_0)) {
        Fail_FlagPos    = true;
    }
    if ((Ui < -negUu_0)) {
        Fail_FlagNeg    = true;
    }

    // If Failure took place in a given direction (Fail_Flag_dir=1), Set the Boundary Moment in the opposite direction to Mr
    if ((Ui < 0) && (Di == 1) && (Fail_FlagNeg)) {
        Mi_boundary     = posFr_0;
    }
    else if ((Ui > 0) && (Di == -1) && (Fail_FlagPos)) {
        Mi_boundary     = -negFr_0;
    }


    // %%%%%%% Current Step Moment Calculation %%%%%%%
    // If current moment based on unloading/reloading Ki is larger than the boundary moment, set it equal to the boundary moment
    if (QuarterFlag == 1 && Di == 1 && Fi > Mi_boundary) {
        Fi  = Mi_boundary;
    }
    else if (QuarterFlag == 3 && Di == -1 && Fi < Mi_boundary) {
        Fi  = Mi_boundary;
    }
    else if (QuarterFlag == 2 && Fi < Mi_boundary) {
        Fi  = Mi_boundary;
    }
    else if (QuarterFlag == 4 && Fi > Mi_boundary) {
        Fi  = Mi_boundary;
    }
// Residual Curveの判定
    if (Mrneg_Flag || Mrpos_Flag) {
        if (QuarterFlag == 1 && Di == 1 && Fi_1 == posFr_0)  {
            Fi  = posFr_0;
        }
        if  (QuarterFlag == 3 && Di == -1 && Fi_1 == -negFr_0) {
            Fi  = -negFr_0;
        }
    }
    // if fail flag is reached in any loading direction, set current moment equal to zero
    if ( Fail_FlagPos || Fail_FlagNeg || Energy_Flag ) {
        Fi  = 0.0;
    }
// Yieldの判定
    if (!Yield_Flag) {
        if (Ui > posUy_0) {
            Fi  = posFy_0 + posKp_0 * (Ui - posUy_0);
            Yield_Flag  = true;
        }
        // else {
        //     Fi  = Ke * (Ui);
        // }

        else if (Ui < -negUy_0) {
            Fi  = -negFy_0 - negKp_0 * fabs(Ui - negUy_0);
        }
        else {
            Fi  = Ke * (Ui);
        }
    }
// エネルギーを計算して、耐力低下の判定
    //cout << "                Fi_1=" << Fi_1 << " Fi=" << Fi << " TangentK=" << TangentK << " Mbound=" << Mi_boundary << " Q=" << QuarterFlag << endln;

    // %%%%%%%%%%%%% Energy Calculation %%%%%%%%%%%%%
    dEi     = (Fi + Fi_1) * 0.5 * dU;
    engTotl += dEi; //  total energy dissipated till current incremental step

    // Energy calculation at each new excursion
    if (Fi * Fi_1 <= 0.0) {
        engExcr       = max(0.,engTotl - engExcr_1);  // total energy dissipated in current excursion
        engExcr_1     = engTotl;                // total energy dissipated in previous excursion
        Excursion_Flag  = true;
    }
    else {
        Excursion_Flag  = false;
    }

    // Check if the Component inherit Reference Energy is Consumed
    if (Excursion_Flag) {
        if ( engTotl > refEnergyS || engTotl > refEnergyC || betaS > 1 || betaC > 1 ) {
            Energy_Flag = true;
        }
    }
    else if (Reversal_Flag) {
        if ( engTotl > refEnergyK || betaK > 1 ) {
            Energy_Flag = true;
        }
    }

    // Tangent Stiffeness Calculation
    if (Fi == posFr_0 || Fi == -negFr_0) {
        TangentK    = pow(10., -6);
    }

    if (Ui == Ui_1) {
        TangentK    = Ke;
        Fi          = Fi_1;
    }
    else {
        TangentK    = (Fi - Fi_1) / dU;
        if (TangentK == 0) {
            TangentK    = pow(10., -6);
        }
    }
    Di_1    = Di;
    //cout << "                Fi_1=" << Fi_1 << " Fi=" << Fi << " Ke=" << Ke << " TangentK=" << TangentK << " Mbound=" << Mi_boundary << endln;

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return 0;
}

double IMKBilin::getStress(void)
{
    //cout << " getStress" << endln;
    return (Fi);
}

double IMKBilin::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (TangentK);
}

double IMKBilin::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (Ke);
}

double IMKBilin::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (U);
}

int IMKBilin::commitState(void)
{
    //cout << " commitState" << endln;

    //commit trial  variables

    cU              = U;

    cUi             = Ui;
    cFi             = Fi;
    // cDi             = Di;
    cUi_1           = Ui_1;
    cFi_1           = Fi_1;
    cDi_1           = Di_1;

    cUlocal         = Ulocal;
    cFlocal         = Flocal;
    cTangentK       = TangentK;

    cK_j          = K_j;
    cPosUy_1        = posUy_1;
    cPosUcap_1      = posUcap_1;
    cPosKp_1        = posKp_1;
    cPosKpc_1       = posKpc_1;
    cPosMy_1        = posFy_1;
    cPosMyProj_1    = posFyProj_1;
    cPosMmax_1      = posFcap_1;
    cPosMmaxProj_1  = posFcapProj_1;

    cNegUy_1	    = negUy_1;
    cNegUcap_1	    = negUcap_1;
    cNegKp_1	= negKp_1;
    cNegKpc_1	= negKpc_1;
    cNegMy_1	    = negFy_1;
    cNegMyProj_1	= negFyProj_1;
    cNegMmax_1	    = negFcap_1;
    cNegMmaxProj_1	= negFcapProj_1;

    // cBetaS	    = betaS;
    // cBetaC	    = betaC;
    // cBetaK	    = betaK;

    cExcursion_Flag	= Excursion_Flag;
    cReversal_Flag	= Reversal_Flag;
    cYield_Flag	    = Yield_Flag;
    cFail_FlagPos	= Fail_FlagPos;
    cFail_FlagNeg	= Fail_FlagNeg;
    cMrpos_Flag	    = Mrpos_Flag;
    cMrneg_Flag	    = Mrneg_Flag;
    cEnergy_Flag	= Energy_Flag;

    cEngExcr_1	= engExcr_1;
    cEngExcr	    = engExcr;
    cEngRvrs	    = engRvrs;
    cEngTotl	    = engTotl;

    return 0;
}

int IMKBilin::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;

    //the opposite of commit trial history variables
    U	            = cU;

    Ui	            = cUi;
    Fi	            = cFi;
    // Di	            = cDi;
    Ui_1	        = cUi_1;
    Fi_1	        = cFi_1;
    Di_1	        = cDi_1;

    Ulocal	    = cUlocal;
    Flocal	    = cFlocal;
    TangentK	    = cTangentK;

    K_j	        = cK_j;
    posUy_1	        = cPosUy_1;
    posUcap_1	    = cPosUcap_1;
    posKp_1	= cPosKp_1;
    posKpc_1	= cPosKpc_1;
    posFy_1	    = cPosMy_1;
    posFyProj_1	= cPosMyProj_1;
    posFcap_1	    = cPosMmax_1;
    posFcapProj_1	= cPosMmaxProj_1;

    negUy_1	        = cNegUy_1;
    negUcap_1	    = cNegUcap_1;
    negKp_1	= cNegKp_1;
    negKpc_1	= cNegKpc_1;
    negFy_1	    = cNegMy_1;
    negFyProj_1	= cNegMyProj_1;
    negFcap_1	    = cNegMmax_1;
    negFcapProj_1	= cNegMmaxProj_1;

    // betaS	    = cBetaS;
    // betaC	    = cBetaC;
    // betaK	    = cBetaK;

    Excursion_Flag	= cExcursion_Flag;
    Reversal_Flag	= cReversal_Flag;
    Yield_Flag	    = cYield_Flag;
    Fail_FlagPos	= cFail_FlagPos;
    Fail_FlagNeg	= cFail_FlagNeg;
    Mrpos_Flag	    = cMrpos_Flag;
    Mrneg_Flag	    = cMrneg_Flag;
    Energy_Flag	    = cEnergy_Flag;

    engExcr_1	    = cEngExcr_1;
    engExcr	    = cEngExcr;
    engRvrs	    = cEngRvrs;
    engTotl	    = cEngTotl;

    return 0;
}

int IMKBilin::revertToStart(void)
{
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE TIME CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    if (posResF_0 == 0.0) {
        posResF_0	= 0.01;
    }
    if (negResF_0 == 0.0) {
        negResF_0	= 0.01;
    }

    posUy_0	                            = posFy_0 / Ke;
    posUcap_0	                        = posUy_0 + posUp_0;
    posKp_0	                    = (posFcap_0 - posFy_0) / (posUp_0);
    posKpc_0	                    = posFcap_0 / (posUpc_0);
    // posFy_0	                        = posFy_0;
    posFcap_0	                        = posFcapFy_0 * posFy_0;
    posFyProj_0	                    = posFcap_0 - posKp_0  * posUcap_0;
    posFcapProj_0	                    = posFcap_0 + posKpc_0 * posUcap_0;

    negUy_0	                            = negFy_0 / Ke;
    negUcap_0	                        = negUy_0 + negUp_0;
    negKp_0	                    = (negFcap_0 - negFy_0) / (negUp_0);
    negKpc_0	                    = negFcap_0 / (negUpc_0);
    // negFy_0	                        = negFy_0;
    negFcap_0	                        = negFcapFy_0 * negFy_0;
    negFyProj_0	                    = negFcap_0 - negKp_0  * negUcap_0;
    negFcapProj_0	                    = negFcap_0 + negKpc_0 * negUcap_0;

    posFr_0	                            = posResF_0*posFy_0;
    negFr_0	                            = negResF_0*negFy_0;

    refEnergyS	                    = LAMBDA_S*posFy_0;
    refEnergyC	                    = LAMBDA_C*posFy_0;
    refEnergyK	                    = LAMBDA_K*posFy_0;

    K_j	        = cK_j	        = Ke;
    posUy_1	        = cPosUy_1	        = posFy_0 / Ke;
    posUcap_1	    = cPosUcap_1	    = posUy_0 + posUp_0;
    posKp_1	= cPosKp_1	    = (posFcap_0 - posFy_0) / (posUp_0);
    posKpc_1	= cPosKpc_1	= posFcap_0 / (posUpc_0);
    posFy_1	    = cPosMy_1	        = posFy_0;
    posFcap_1	    = cPosMmax_1	    = posFcapFy_0*posFy_0;
    posFyProj_1	= cPosMyProj_1	    = posFcap_1 - posKp_1 *posUcap_1;
    posFcapProj_1	= cPosMmaxProj_1	= posFcap_1 + posKpc_1*posUcap_1;

    negUy_1	        = cNegUy_1	        = negFy_0 / Ke;
    negUcap_1	    = cNegUcap_1	    = negUy_0 + negUp_0;
    negKp_1	= cNegKp_1	    = (negFcap_0 - negFy_0) / (negUp_0);
    negKpc_1	= cNegKpc_1	= negFcap_0 / (negUpc_0);
    negFy_1	    = cNegMy_1	        = negFy_0;
    negFcap_1	    = cNegMmax_1	    = negFcapFy_0*negFy_0;
    negFyProj_1	= cNegMyProj_1	    = negFcap_1 - negKp_1 *negUcap_1;
    negFcapProj_1	= cNegMmaxProj_1	= negFcap_1 + negKpc_1*negUcap_1;

    //initially I zero everything   
    U	            = cU	            = 0;

    Ui	            = cUi	            = 0;
    Fi	            = cFi	            = 0;
    // Di	            = cDi	            = 0;
    Ui_1	        = cUi_1	            = 0;
    Fi_1	        = cFi_1	            = 0;
    Di_1	        = cDi_1	            = 0;
    TangentK	    = cTangentK	        = Ke;

    Ulocal	    = cUlocal	    = 0;
    Flocal	    = cFlocal	    = 0;

    // betaS	        = cBetaS	        = 0;
    // betaC	        = cBetaC	        = 0;
    // betaK	        = cBetaK	        = 0;

    Excursion_Flag	= cExcursion_Flag	= false;
    Reversal_Flag	= cReversal_Flag	= false;
    Yield_Flag	    = cYield_Flag	    = false;
    Fail_FlagPos	= cFail_FlagPos	    = false;
    Fail_FlagNeg	= cFail_FlagNeg	    = false;
    Mrpos_Flag	    = cMrpos_Flag	    = false;
    Mrneg_Flag	    = cMrneg_Flag	    = false;
    Energy_Flag	    = cEnergy_Flag	    = false;

    engExcr_1	    = cEngExcr_1	    = 0;
    engExcr	    = cEngExcr	    = 0;
    engRvrs	    = cEngRvrs	    = 0;
    engTotl	    = cEngTotl	    = 0;

    return 0;
}

UniaxialMaterial *
IMKBilin::getCopy(void)
{
    IMKBilin *theCopy	= new IMKBilin(this->getTag(), Ke,
        posUp_0, posUpc_0, posUu_0, posFy_0, posFcapFy_0, posResF_0,
        negUp_0, negUpc_0, negUu_0, negFy_0, negFcapFy_0, negResF_0,
        LAMBDA_S, LAMBDA_C, LAMBDA_K, c_S, c_C, c_K, D_pos, D_neg);

    //cout << " getCopy" << endln;

    theCopy->posFr_0	        = posFr_0;
    theCopy->negFr_0	        = negFr_0;

    theCopy->U	                = U;
    theCopy->cU	                = cU;

    theCopy->Ui	                = Ui;
    theCopy->Fi	                = Fi;
    // theCopy->Di	                = Di;
    theCopy->Ui_1	            = Ui_1;
    theCopy->Fi_1	            = Fi_1;
    theCopy->Di_1	            = Di_1;

    theCopy->Ulocal	        = Ulocal;
    theCopy->Flocal	        = Flocal;
    theCopy->TangentK	        = TangentK;

    theCopy->K_j	            = K_j;
    theCopy->posUy_1	        = posUy_1;
    theCopy->posUcap_1	        = posUcap_1;
    theCopy->posKp_1	    = posKp_1;
    theCopy->posKpc_1	    = posKpc_1;
    theCopy->posFy_1	        = posFy_1;
    theCopy->posFyProj_1	    = posFyProj_1;
    theCopy->posFcap_1	        = posFcap_1;
    theCopy->posFcapProj_1	    = posFcapProj_1;

    theCopy->negUy_1	        = negUy_1;
    theCopy->negUcap_1	        = negUcap_1;
    theCopy->negKp_1	    = negKp_1;
    theCopy->negKpc_1	    = negKpc_1;
    theCopy->negFy_1	        = negFy_1;
    theCopy->negFyProj_1	    = negFyProj_1;
    theCopy->negFcap_1	        = negFcap_1;
    theCopy->negFcapProj_1	    = negFcapProj_1;

    // theCopy->betaS	        = betaS;
    // theCopy->betaC	        = betaC;
    // theCopy->betaK	        = betaK;

    theCopy->Excursion_Flag	    = Excursion_Flag;
    theCopy->Reversal_Flag	    = Reversal_Flag;
    theCopy->Yield_Flag	        = Yield_Flag;
    theCopy->Fail_FlagPos	    = Fail_FlagPos;
    theCopy->Fail_FlagNeg	    = Fail_FlagNeg;
    theCopy->Mrpos_Flag	        = Mrpos_Flag;
    theCopy->Mrneg_Flag	        = Mrneg_Flag;
    theCopy->Energy_Flag	    = Energy_Flag;

    theCopy->engExcr_1	    = engExcr_1;
    theCopy->engExcr	        = engExcr;
    theCopy->engRvrs	        = engRvrs;
    theCopy->engTotl	        = engTotl;


    theCopy->cPosFr_0	        = cPosFr_0;
    theCopy->cNegFr_0	        = cNegFr_0;

    theCopy->cUi	            = cUi;
    theCopy->cFi	            = cFi;
    // theCopy->cDi	            = cDi;
    theCopy->cUi_1	            = cUi_1;
    theCopy->cFi_1	            = cFi_1;
    theCopy->cDi_1	            = cDi_1;

    theCopy->cUlocal	        = cUlocal;
    theCopy->cFlocal	        = cFlocal;
    theCopy->cTangentK	        = cTangentK;

    theCopy->cK_j	            = cK_j;
    theCopy->cPosUy_1	        = cPosUy_1;
    theCopy->cPosUcap_1	        = cPosUcap_1;
    theCopy->cPosKp_1	    = cPosKp_1;
    theCopy->cPosKpc_1	    = cPosKpc_1;
    theCopy->cPosMy_1	        = cPosMy_1;
    theCopy->cPosMyProj_1	    = cPosMyProj_1;
    theCopy->cPosMmax_1	        = cPosMmax_1;
    theCopy->cPosMmaxProj_1	    = cPosMmaxProj_1;

    theCopy->cNegUy_1	        = cNegUy_1;
    theCopy->cNegUcap_1	        = cNegUcap_1;
    theCopy->cNegKp_1	    = cNegKp_1;
    theCopy->cNegKpc_1	    = cNegKpc_1;
    theCopy->cNegMy_1	        = cNegMy_1;
    theCopy->cNegMyProj_1	    = cNegMyProj_1;
    theCopy->cNegMmax_1	        = cNegMmax_1;
    theCopy->cNegMmaxProj_1	    = cNegMmaxProj_1;

    // theCopy->cBetaS	        = cBetaS;
    // theCopy->cBetaC	        = cBetaC;
    // theCopy->cBetaK	        = cBetaK;

    theCopy->cExcursion_Flag	= cExcursion_Flag;
    theCopy->cReversal_Flag	    = cReversal_Flag;
    theCopy->cYield_Flag	    = cYield_Flag;
    theCopy->cFail_FlagPos	    = cFail_FlagPos;
    theCopy->cFail_FlagNeg	    = cFail_FlagNeg;
    theCopy->cMrpos_Flag	    = cMrpos_Flag;
    theCopy->cMrneg_Flag	    = cMrneg_Flag;
    theCopy->cEnergy_Flag	    = cEnergy_Flag;

    theCopy->cEngExcr_1	    = cEngExcr_1;
    theCopy->cEngExcr	        = cEngExcr;
    theCopy->cEngRvrs	        = cEngRvrs;
    theCopy->cEngTotl	        = cEngTotl;

    return theCopy;
}

int IMKBilin::sendSelf(int cTag, Channel &theChannel)
{
    int res	= 0;
    cout << " sendSelf" << endln;

    static Vector data(113);
    data(0)	    = this->getTag();
    data(1)	    = Ke;
    data(2)	    = posUp_0;
    data(3)	    = posUpc_0;
    data(4)	    = posUu_0;
    data(5)	    = posFy_0;
    data(6)	    = posFcapFy_0;
    data(7)	    = posResF_0;
    data(8)	    = negUp_0;
    data(9)	    = negUpc_0;
    data(10)	= negUu_0;
    data(11)	= negFy_0;
    data(12)	= negFcapFy_0;
    data(13)	= negResF_0;
    data(14)	= LAMBDA_S;
    data(15)	= LAMBDA_C;
    data(16)	= LAMBDA_K;
    data(17)	= c_S;
    data(18)	= c_C;
    data(19)	= c_K;
    data(20)	= D_pos;
    data(21)	= D_neg;

    data(22)	= Ui;
    data(23)	= Fi;
    // data(24)	= Di;
    data(25)	= Ui_1;
    data(26)	= Fi_1;
    data(27)	= Di_1;
    data(28)	= Ulocal;
    data(29)	= Flocal;
    data(30)	= TangentK;

    data(31)	= K_j;
    data(32)	= posUy_1;
    data(33)	= posUcap_1;
    data(34)	= posKp_1;
    data(35)	= posKpc_1;
    data(36)	= posFy_1;
    data(37)	= posFyProj_1;
    data(38)	= posFcap_1;
    data(39)	= posFcapProj_1;

    data(40)	= negUy_1;
    data(41)	= negUcap_1;
    data(42)	= negKp_1;
    data(43)	= negKpc_1;
    data(44)	= negFy_1;
    data(45)	= negFyProj_1;
    data(46)	= negFcap_1;
    data(47)	= negFcapProj_1;

    // data(48)	= betaS;
    // data(49)	= betaC;
    // data(50)	= betaK;

    data(51)	= refEnergyS;
    data(52)	= refEnergyC;
    data(53)	= refEnergyK;

    data(54)	= Excursion_Flag;
    data(55)	= Reversal_Flag;
    data(56)	= Yield_Flag;
    data(57)	= Fail_FlagPos;
    data(58)	= Fail_FlagNeg;
    data(59)	= Mrpos_Flag;
    data(60)	= Mrneg_Flag;
    data(61)	= Energy_Flag;

    data(62)	= engExcr_1;
    data(63)	= engExcr;
    data(64)	= engRvrs;
    data(65)	= engTotl;

    data(66)	= cUi;
    data(67)	= cFi;
    // data(68)	= cDi;
    data(69)	= cUi_1;
    data(70)	= cFi_1;
    data(71)	= cDi_1;
    data(72)	= cUlocal;
    data(73)	= cFlocal;
    data(74)	= cTangentK;

    data(75)	= cK_j;
    data(76)	= cPosUy_1;
    data(77)	= cPosUcap_1;
    data(78)	= cPosKp_1;
    data(79)	= cPosKpc_1;
    data(80)	= cPosMy_1;
    data(81)	= cPosMyProj_1;
    data(82)	= cPosMmax_1;
    data(83)	= cPosMmaxProj_1;

    data(84)	= cNegUy_1;
    data(85)	= cNegUcap_1;
    data(86)	= cNegKp_1;
    data(87)	= cNegKpc_1;
    data(88)	= cNegMy_1;
    data(89)	= cNegMyProj_1;
    data(90)	= cNegMmax_1;
    data(91)	= cNegMmaxProj_1;

    // data(92)	= cBetaS;
    // data(93)	= cBetaC;
    // data(94)	= cBetaK;

    data(95)	= cExcursion_Flag;
    data(96)	= cReversal_Flag;
    data(97)	= cYield_Flag;
    data(98)	= cFail_FlagPos;
    data(99)	= cFail_FlagNeg;
    data(100)	= cMrpos_Flag;
    data(101)	= cMrneg_Flag;
    data(102)	= cEnergy_Flag;

    data(103)	= cEngExcr_1;
    data(104)	= cEngExcr;
    data(105)	= cEngRvrs;
    data(106)	= cEngTotl;

    data(107)	= posFr_0;
    data(108)	= negFr_0;
    data(109)	= cPosFr_0;
    data(110)	= cNegFr_0;

    data(111)	= U;
    data(112)	= cU;

    res	= theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "IMKBilin::sendSelf() - failed to send data\n";

    return res;
}

int IMKBilin::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res	= 0;
    static Vector data(113);
    res	= theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "IMKBilin::recvSelf() - failed to receive data\n";
        this->setTag(0);
    }
    else {
        cout << " recvSelf" << endln;
        this->setTag((int)data(0));
        Ke	            = data(1);
        posUp_0	        = data(2);
        posUpc_0	    = data(3);
        posUu_0	        = data(4);
        posFy_0	    = data(5);
        posFcapFy_0	= data(6);
        posResF_0	    = data(7);
        negUp_0	        = data(8);
        negUpc_0	    = data(9);
        negUu_0	        = data(10);
        negFy_0	    = data(11);
        negFcapFy_0	= data(12);
        negResF_0	    = data(13);
        LAMBDA_S	    = data(14);
        LAMBDA_C	    = data(15);
        LAMBDA_K	    = data(16);
        c_S	            = data(17);
        c_C	            = data(18);
        c_K	            = data(19);
        D_pos	        = data(20);
        D_neg	        = data(21);

        Ui	            = data(22);
        Fi	            = data(23);
        // Di	            = data(24);
        Ui_1	        = data(25);
        Fi_1	        = data(26);
        Di_1	        = data(27);
        Ulocal	    = data(28);
        Flocal	    = data(29);
        TangentK	    = data(30);

        K_j	        = data(31);
        posUy_1	        = data(32);
        posUcap_1	    = data(33);
        posKp_1	= data(34);
        posKpc_1	= data(35);
        posFy_1	    = data(36);
        posFyProj_1	= data(37);
        posFcap_1	    = data(38);
        posFcapProj_1	= data(39);

        negUy_1	        = data(40);
        negUcap_1	    = data(41);
        negKp_1	= data(42);
        negKpc_1	= data(43);
        negFy_1	    = data(44);
        negFyProj_1	= data(45);
        negFcap_1	    = data(46);
        negFcapProj_1	= data(47);

        // betaS	        = data(48);
        // betaC	        = data(49);
        // betaK	        = data(50);

        refEnergyS	= data(51);
        refEnergyC	= data(52);
        refEnergyK	= data(53);

        Excursion_Flag	= data(54);
        Reversal_Flag	= data(55);
        Yield_Flag	    = data(56);
        Fail_FlagPos	= data(57);
        Fail_FlagNeg	= data(58);
        Mrpos_Flag	    = data(59);
        Mrneg_Flag	    = data(60);
        Energy_Flag	    = data(61);

        engExcr_1	    = data(62);
        engExcr	    = data(63);
        engRvrs	    = data(64);
        engTotl	    = data(65);

        cUi	            = data(66);
        cFi	            = data(67);
        // cDi	            = data(68);
        cUi_1	        = data(69);
        cFi_1	        = data(70);
        cDi_1	        = data(71);

        cUlocal	    = data(72);
        cFlocal	    = data(73);
        cTangentK	    = data(74);

        cK_j	        = data(75);
        cPosUy_1	    = data(76);
        cPosUcap_1	    = data(77);
        cPosKp_1	= data(78);
        cPosKpc_1	= data(79);
        cPosMy_1	    = data(80);
        cPosMyProj_1	= data(81);
        cPosMmax_1	    = data(82);
        cPosMmaxProj_1	= data(83);

        cNegUy_1	    = data(84);
        cNegUcap_1	    = data(85);
        cNegKp_1	= data(86);
        cNegKpc_1	= data(87);
        cNegMy_1	    = data(88);
        cNegMyProj_1	= data(89);
        cNegMmax_1	    = data(90);
        cNegMmaxProj_1	= data(91);

        // cBetaS	    = data(92);
        // cBetaC	    = data(93);
        // cBetaK	    = data(94);

        cExcursion_Flag	= data(95);
        cReversal_Flag	= data(96);
        cYield_Flag	    = data(97);
        cFail_FlagPos	= data(98);
        cFail_FlagNeg	= data(99);
        cMrpos_Flag	    = data(100);
        cMrneg_Flag	    = data(101);
        cEnergy_Flag	= data(102);

        cEngExcr_1	= data(103);
        cEngExcr	    = data(104);
        cEngRvrs	    = data(105);
        cEngTotl	    = data(106);

        posFr_0	        = data(107);
        negFr_0	        = data(108);
        cPosFr_0	    = data(109);
        cNegFr_0	    = data(110);

        U	            = data(111);
        cU	            = data(112);

    }

    return res;
}

void IMKBilin::Print(OPS_Stream &s, int flag)
{
    cout << "IMKBilin tag: " << this->getTag() << endln;
}
