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
#include <IMKPeakOriented_n.h>
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
static int numIMKPeakOrientedMaterials	= 0;
void *
OPS_IMKPeakOriented()
{
	if (numIMKPeakOrientedMaterials == 0) {
		numIMKPeakOrientedMaterials++;
		OPS_Error("IMK with Peak-Oriented Response - Code by Elkady & Eljisr (July22)\n", 1);
	}
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial	= 0;
	int    iData[1];
	double dData[23];
	int numData	= 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial IMKPeakOriented tag" << endln;
		return 0;
	}
	numData	= 23;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid Args want: uniaxialMaterial IMKPeakOriented tag? Ke? ";
		opserr << "Up_pos? Upc_pos? Uu_pos? Fy_pos? FcapFy_pos? ResF_pos? ";
		opserr << "Up_neg? Upc_neg? Uu_neg? Fy_neg? FcapFy_neg? ResF_neg? ";
		opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? ";
		return 0;
	}
	// Parsing was successful, allocate the material
	theMaterial	= new IMKPeakOriented(iData[0],
		dData[0],
		dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
		dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20],
		dData[21], dData[22]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type IMKPeakOriented Material\n";
		return 0;
	}
	return theMaterial;
}
IMKPeakOriented::IMKPeakOriented(int tag, double p_Ke,
	double p_Up_pos, double p_Upc_pos, double p_Uu_pos, double p_Fy_pos, double p_FcapFy_pos, double p_ResF_pos,
	double p_Up_neg, double p_Upc_neg, double p_Uu_neg, double p_Fy_neg, double p_FcapFy_neg, double p_ResF_neg,
	double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_A, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_A, double p_c_K, double p_D_pos, double p_D_neg)
	: UniaxialMaterial(tag, 0), Ke(p_Ke),
	Up_pos(p_Up_pos), Upc_pos(p_Upc_pos), Uu_pos(p_Uu_pos), Fy_pos(p_Fy_pos), FcapFy_pos(p_FcapFy_pos), ResF_pos(p_ResF_pos),
	Up_neg(p_Up_neg), Upc_neg(p_Upc_neg), Uu_neg(p_Uu_neg), Fy_neg(p_Fy_neg), FcapFy_neg(p_FcapFy_neg), ResF_neg(p_ResF_neg),
	LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_A(p_LAMBDA_A), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_A(p_c_A), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
	this->revertToStart();
}
IMKPeakOriented::IMKPeakOriented()
	:UniaxialMaterial(0, 0), Ke(0),
	Up_pos(0), Upc_pos(0), Uu_pos(0), Fy_pos(0), FcapFy_pos(0), ResF_pos(0),
	Up_neg(0), Upc_neg(0), Uu_neg(0), Fy_neg(0), FcapFy_neg(0), ResF_neg(0),
	LAMBDA_S(0), LAMBDA_C(0), LAMBDA_A(0), LAMBDA_K(0), c_S(0), c_C(0), c_A(0), c_K(0), D_pos(0), D_neg(0)
{
	this->revertToStart();
}
IMKPeakOriented::~IMKPeakOriented()
{
	// does nothing
}
int IMKPeakOriented::setTrialStrain(double strain, double strainRate)
{
	//all variables to the last commit
	this->revertToLastCommit();
	//state determination algorithm: defines the current force and tangent stiffness
	U		= strain; //set trial displacement
	ui_1	= ui;
	fi_1	= fi;
	ui		= U;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	du		= ui - ui_1;	// Incremental Deformation at Current Step
	if (Failure_Flag) {		// When a failure has already occured
		fi	= 0;
		dEi	= 0;
	} else if (du == 0) {	// When deformation doesn't change from the last
		fi	= fi_1;
		dEi	= 0;
	} else {
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		////////////////// BRANCH DETERMINATION AND FLAG RAISE ////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
	// 	Branch 
	// 		0:	Elastic
	// 		1:	Unloading Branch
	// 		3:	Towards Local Peak		+
	// 		4:	Towards Global Peak		+
	// 		5:	Towards Capping Point	+
	// 		6:	Towards Residual Point	+
	// 		7:	Residual Branch			+
	// 		13:	Towards Local Peak		-
	// 		14:	Towards Global Peak		-
	// 		15:	Towards Capping Point	-
	// 		16:	Towards Residual Point	-
	// 		17:	Residual Branch			-
	// 	Flag
	// 		Yield_Flag:		Preserved.		When the deformation exceeds yield capacity for the first time.
	// 		Excursion_Flag:	Not preserved.	When crossing X-axis. Evokes re-considering of the deteriorations and which peak to go for.
	// 		Reversal_Flag:	Not preserved.	When unloading starts. Evokes re-condiersing of the stiffness deterioration and peak point registration.
		exBranch		= Branch;
		Excursion_Flag	= false;
		Reversal_Flag	= false;
		if (Branch == 0) {
			// CHECK FOR YIELDING
			if (ui > posUy_1) {
				Yield_Flag	= true;
				Branch	= 5;
			} else if (ui < negUy_1) {
				Yield_Flag	= true;
				Branch	= 15;
			}
		} else if (Branch == 1) {
			if (fi_1*(fi_1+du*K_unload) <= 0) {
			// CHECK FOR NEW EXCURSION
				Excursion_Flag	= true;
			} else if (ui > posULocal_1) {
				Branch	= 4;
			} else if (ui < negULocal_1) {
				Branch	= 14;
			}
		} else if (fi_1*du < 0) {
			Reversal_Flag	= true;
			Branch	= 1;
		}
	// Branch shifting from 3 -> 4 -> 5 -> 6 -> 7 can be considered.
		if (Branch == 3 && ui > posULocal_1) {
			Branch	= 4;
		}
		if (Branch == 4 && ui > posUGlobal_1) {
			Branch	= 5;
		}
		if (Branch == 5 && ui > posUcap_1) {
			Branch	= 6;
		}
		if (Branch == 6 && ui > posUres_1) {
			Branch	= 7;
		}
		if (Branch == 13 && ui < negULocal_1) {
			Branch	= 14;
		}
		if (Branch == 14 && ui < negUGlobal_1) {
			Branch	= 15;
		}
		if (Branch == 15 && ui < negUcap_1) {
			Branch	= 16;
		}
		if (Branch == 16 && ui < negUres_1) {
			Branch	= 17;
		}
	// UPDATE PEAK POINTS
		if (Reversal_Flag) {
			if ( fi_1 > 0 ){
				posULocal_1	= ui_1;				// UPDATE LOCAL
				posFLocal_1	= fi_1;
				if ( ui_1 > posUGlobal_1 ) {	// UPDATE GLOBAL
					posUGlobal_1	= ui_1;
					posFGlobal_1	= fi_1;
				}			
			} else {
				negULocal_1	= ui_1;				// UPDATE LOCAL
				negFLocal_1	= fi_1;
				if ( ui_1 < negUGlobal_1 ) {	// UPDATE GLOBAL
					negUGlobal_1	= ui_1;
					negFGlobal_1	= fi_1;
				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// UPDATE DETERIORATION PARAMETERS ///////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		if (Reversal_Flag) {
			EpjK	= Energy_Acc				- 0.5*(fi_1 / K_unload)*fi_1;
			EiK		= Energy_Acc - Energy_Diss	- 0.5*(fi_1 / K_unload)*fi_1;
			betaK	= pow( (EiK / (EtK - EpjK)), c_K );
			K_unload	= K_unload * (1 - betaK);
		// Detect unloading completed in a step.
			if (fi_1*(fi_1+du*K_unload) <= 0) {
				Excursion_Flag	= true;
				Reversal_Flag	= false;
			}
		} else {
			betaK	= 0;
		}
		if (Excursion_Flag) {
			Ei		= fmax( 0, Energy_Acc - Energy_Diss);
			betaS	= pow( (Ei / (EtS - Energy_Acc)), c_S);
			betaC	= pow( (Ei / (EtC - Energy_Acc)), c_C);
			betaA	= pow( (Ei / (EtA - Energy_Acc)), c_A);
			Energy_Diss	= Energy_Acc;
		} else {
			betaS	= 0;
			betaC	= 0;
			betaA	= 0;
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		if ( Excursion_Flag && Yield_Flag ) {
			// Positive
			if (fi_1 < 0) {
				// Basic strength deterioration: Yield point
				// Basic strength deterioration: Post-yield Stiffness
				posFy_1	= posFy_1	* (1 - betaS * D_pos);
				posKp_1	= posKp_1	* (1 - betaS * D_pos);
				if (posFy_1 < posFres_1) {
					posFy_1	= posFres_1;
					posKp_1	= 0;
				}
				posUy_1	= posFy_1 / Ke;
				// Basic strength deterioration: Capping Point
				sPCsp		= (posFy_1 - posUy_1 * posKp_1 - posFcap_1 + posKpc_1 * posUcap_1) / (posKpc_1 - posKp_1);
				posFcap_1	= posFcap_1 + (sPCsp - posUcap_1)*posKpc_1;
				posUcap_1	= sPCsp;
				// Post-capping strength deterioration: Capping point
				sPCpcp		= max(posUcap_1 + betaC * D_pos*(posFcap_1 - posKpc_1 * posUcap_1) / (posKpc_1 - posKp_1), posUy_1);
				posFcap_1	= posFcap_1 + (sPCpcp - posUcap_1)*posKp_1;
				posUcap_1	= sPCpcp;
				// Accelerated reloading stiffness deterioration: Target peak deformation point
				posUGlobal_1	= (1 + betaA * D_pos)*posUGlobal_1;
				if (posUGlobal_1 < posUy_1) {
					posFGlobal_1	= Ke * posUGlobal_1;
					// Target peak deformation in post-yield branch of the updated backbone
				} else if (posUGlobal_1 < posUcap_1) {
					posFGlobal_1	= posKp_1 * (posUGlobal_1 - posUy_1) + posFy_1;
					// Target peak deformation in post-capping branch of the updated backbone
				} else {
					posFGlobal_1	= max(posKpc_1*(posUGlobal_1 - posUcap_1) + posFcap_1, posFres_1);
				}
				posUres_1	= (posFres_1 - posFcap_1 + posKpc_1 * posUcap_1) / posKpc_1;
			// Negative
			} else {
				// Basic strength deterioration: Yield point
				// Basic strength deterioration: Post-yield stiffness
				negFy_1	= negFy_1	* (1 - betaS * D_neg);
				negKp_1	= negKp_1	* (1 - betaS * D_neg);
				if (negFy_1 > negFres_1) {
					negFy_1	= negFres_1;
					negKp_1	= 0;
				}
				negUy_1	= negFy_1 / Ke;
				// Basic strength deterioration: Capping point
				sPCsn		= (negFy_1 - negUy_1 * negKp_1 - negFcap_1 + negKpc_1 * negUcap_1) / (negKpc_1 - negKp_1);
				negFcap_1	= negFcap_1 + (sPCsn - negUcap_1)*negKpc_1;
				negUcap_1	= sPCsn;
				// Post-capping strength deterioration: Capping point
				sPCpcn		= min(negUcap_1 + betaC * D_neg*(negFcap_1 - negKpc_1 * negUcap_1) / (negKpc_1 - negKp_1), negUy_1);
				negFcap_1	= negFcap_1 + (sPCpcn - negUcap_1)*negKp_1;
				negUcap_1	= sPCpcn;
				// Accelerated reloading stiffness deterioration: Target peak deformation point
				negUGlobal_1	= (1 + betaA * D_neg)*negUGlobal_1;
				// Target peak deformation in reloading branch of the updated backbone
				if (negUGlobal_1 > negUy_1) {
					negFGlobal_1	= Ke * negUGlobal_1;
					// Target peak deformation in post-yield branch of the updated backbone
				} else if (negUGlobal_1 > negUcap_1) {
					negFGlobal_1	= negKp_1 * (negUGlobal_1 - negUy_1) + negFy_1;
					// Target peak deformation in post-capping branch of the updated backbone
				} else {
					negFGlobal_1	= min(negKpc_1*(negUGlobal_1 - negUcap_1) + negFcap_1, negFres_1);
				}
				negUres_1	= (negFres_1 - negFcap_1 + negKpc_1 * negUcap_1) / negKpc_1;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////// COMPUTE FORCE INCREMENT /////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		// CASE 1: EACH NEW EXCURSION
		if (Excursion_Flag) {
			// Detection of reloading completed in a step might be needed, while it's not as severe as a one step unloading.
			u0	= ui_1 - (fi_1 / K_unload);
			if (du > 0) {
				K_Local		= posFLocal_1	/ (posULocal_1	- u0);
				K_Global	= posFGlobal_1	/ (posUGlobal_1	- u0);
				if ( (posFLocal_1 < posFGlobal_1) && (K_Local > K_Global)) {
					Branch		= 3;
					K_reload	= K_Local;
				} else {
					Branch		= 4;
					K_reload	= K_Global;
				}
			} else {
				K_Local		= negFLocal_1	/ (negULocal_1	- u0);
				K_Global	= negFGlobal_1	/ (negUGlobal_1	- u0);
				if ( (negFLocal_1 > negFGlobal_1) && (K_Local > K_Global)) {
					Branch		= 13;
					K_reload	= K_Local;
				} else {
					Branch		= 14;
					K_reload	= K_Global;
				}
			}
			df	= 0				- fi_1 + K_reload*	(ui - u0);

// With Branch Change
	// Positive Force
		} else if (Branch == 4 && exBranch != 4) {
			K_reload	= (posFGlobal_1 - posFLocal_1) / (posUGlobal_1 - posULocal_1);
			df	= posFLocal_1	- fi_1 + K_reload*	(ui - posULocal_1);
		} else if (Branch == 5 && exBranch == 0) {
			df	= posFy_1		- fi_1 + posKp_1*	(ui - posUy_1);
		} else if (Branch == 5 && exBranch != 5) {
			df	= posFGlobal_1	- fi_1 + posKp_1*	(ui - posUGlobal_1);
		} else if (Branch == 6 && exBranch == 5) {
			df	= posFcap_1		- fi_1 + posKpc_1*	(ui - posUcap_1);
		} else if (Branch == 6 && exBranch != 6) {
			df	= posFGlobal_1	- fi_1 + posKpc_1*	(ui - posUGlobal_1);
		} else if (Branch == 7 && exBranch != 7) {
			df	= posFres_1		- fi_1;
	// Negative Force
		} else if (Branch == 14 && exBranch != 14) {
			K_reload	= (negFGlobal_1 - negFLocal_1) / (negUGlobal_1 - negULocal_1);
			df			= negFLocal_1 - fi_1 + K_reload*(ui - negULocal_1);
		} else if (Branch == 15 && exBranch == 0) {
			df	= negFy_1 - fi_1 + negKp_1*(ui - negUy_1);
		} else if (Branch == 15 && exBranch != 15) {
			df	= negFGlobal_1 - fi_1 + negKp_1*(ui - negUGlobal_1);
		} else if (Branch == 16 && exBranch == 15) {
			df	= negFcap_1 - fi_1 + negKpc_1*(ui - negUcap_1);
		} else if (Branch == 16 && exBranch != 16) {
			df	= negFGlobal_1 - fi_1 + negKpc_1*(ui - negUGlobal_1);
		} else if (Branch == 17 && exBranch != 17) {
			df	= negFres_1 - fi_1;
// Without Branch Change
	// Positive Force
		// CASE 0: AT THE ELASTIC SLOPE
		} else if (Branch == 0) {
			df	= du*Ke;
		// CASE 2: WHEN RELOADING
		// CASE 3: WHEN UNLOADING
		} else if (Branch == 1) {
			df	= du*K_unload;
		// CASE 4: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
		// CASE 5: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
		// CASE 6: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
		} else if (Branch == 3 || Branch == 4 || Branch == 13 || Branch == 14) {
			df	= du*K_reload;
		// CASE 7: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
		} else if (Branch == 5) {
			df	= du*posKp_1;
		// CASE 8: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
		} else if (Branch == 6) {
			df	= du*posKpc_1;
		// CASE 9: WHEN LOADING AND BEYOND THE RESIDUAL POINT
		} else if (Branch == 7) {
			df	= 0.0;
	// Negative Force
		// CASE 7: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
		} else if (Branch == 15) {
			df	= du*negKp_1;
		// CASE 8: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
		} else if (Branch == 16) {
			df	= du*negKpc_1;
		// CASE 9: WHEN LOADING AND BEYOND THE RESIDUAL POINT
		} else if (Branch == 17) {
			df	= 0.0;
		}
	// Branch Change check
		// if (Branch!=exBranch) {
		// 	std::cout << exBranch << " -> " << Branch << "\n";
		// }
		fi	= fi_1+df;
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		// CHECK FOR FAILURE
		///////////////////////////////////////////////////////////////////////////////////////////		
		///////////////////////////////////////////////////////////////////////////////////////////		
		///////////////////////////////////////////////////////////////////////////////////////////		

		// Failure criteria (Tolerance = 1//)
	// I have no idea about why it can' t be 0 nor 1.
		FailS	= ( betaS < -0.01 || betaS > 1.01	);
		FailC	= ( betaC < -0.01 || betaC > 1.01	);
		FailA	= ( betaA < -0.01 || betaA > 1.01	);
		FailK	= ( betaK < -0.01 || betaK > 1.01	);
		FailPp	= ( posFGlobal_1 == 0				);
		FailPn	= ( negFGlobal_1 == 0				);
		FailDp	= ( ui >  Uu_pos					);
		FailDn	= ( ui < -Uu_neg					);	
		FailRp	= ( Branch ==  7 && posFres_1 == 0	);
		FailRn	= ( Branch == 17 && negFres_1 == 0	);
		if (FailS||FailC||FailA||FailK||FailPp||FailPn||FailRp||FailRn||FailDp||FailDn) {
			Failure_Flag	= true;
		}
		if (Failure_Flag) {
			fi	= 0;
		}
		dEi	= 0.5*(fi + fi_1)*du; // Internal energy increment
	}
	Energy_Acc	= Energy_Acc + dEi; 	// Energy
	du_i_1		= du;					// Update Displacement
	if ( du == 0 ) {					// Stiffness Calculation
		ki			= Ke;
		TangentK	= Ke;
	} else {
		ki			= (fi - fi_1) / (du);
		TangentK	= (fi - fi_1) / (du);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////// END OF MAIN CODE ///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return 0;
}
double IMKPeakOriented::getStress(void)
{
	//cout << " getStress" << endln;
	return (fi);
}
double IMKPeakOriented::getTangent(void)
{
	//cout << " getTangent" << endln;
	return (TangentK);
}
double IMKPeakOriented::getInitialTangent(void)
{
	//cout << " getInitialTangent" << endln;
	return (Ke);
}
double IMKPeakOriented::getStrain(void)
{
	//cout << " getStrain" << endln;
	return (U);
}
int IMKPeakOriented::commitState(void)
{
	cU				= U;
	cui				= ui;
	cfi				= fi;
	cui_1			= ui_1;
	cfi_1			= fi_1;
	cTangentK		= TangentK;
	cdu_i_1			= du_i_1;
	cposUy_1		= posUy_1;
	cposUcap_1		= posUcap_1;
	cposFy_1		= posFy_1;
	cposFcap_1		= posFcap_1;
	cposUGlobal_1	= posUGlobal_1;
	cposFGlobal_1	= posFGlobal_1;
	cposUres_1		= posUres_1;
	cposFres_1		= posFres_1;
	cposKp_1		= posKp_1;
	cposKpc_1		= posKpc_1;
	cnegUy_1		= negUy_1;
	cnegUcap_1		= negUcap_1;
	cnegFy_1		= negFy_1;
	cnegFcap_1		= negFcap_1;
	cnegUGlobal_1	= negUGlobal_1;
	cnegFGlobal_1	= negFGlobal_1;
	cnegUres_1		= negUres_1;
	cnegFres_1		= negFres_1;
	cnegKp_1		= negKp_1;
	cnegKpc_1		= negKpc_1;
	cK_unload		= K_unload;
	cK_reload		= K_reload;
	cEnergy_Acc		= Energy_Acc;
	cEnergy_Diss	= Energy_Diss;
	// cu0				= u0;
	cposULocal_1	= posULocal_1;
	cposFLocal_1	= posFLocal_1;
	cnegULocal_1	= negULocal_1;
	cnegFLocal_1	= negFLocal_1;
	cFailure_Flag	= Failure_Flag;
	// cExcursion_Flag	= Excursion_Flag;
	cexBranch		= exBranch;
	cBranch			= Branch;
	// cTargetPeak_Flag= TargetPeak_Flag;
	cYield_Flag		= Yield_Flag;
	// cReversal_Flag	= Reversal_Flag;
	return 0;
}
int IMKPeakOriented::revertToLastCommit(void)
{
	//cout << " revertToLastCommit" << endln;
	//the opposite of commit trial history variables
	U				= cU;
	ui				= cui;
	fi				= cfi;
	ui_1			= cui_1;
	fi_1			= cfi_1;
	TangentK		= cTangentK;
	du_i_1			= cdu_i_1;
	posUy_1			= cposUy_1;
	posUcap_1		= cposUcap_1;
	posFy_1			= cposFy_1;
	posFcap_1		= cposFcap_1;
	posUGlobal_1	= cposUGlobal_1;
	posFGlobal_1	= cposFGlobal_1;
	posUres_1		= cposUres_1;
	posFres_1		= cposFres_1;
	posKp_1			= cposKp_1;
	posKpc_1		= cposKpc_1;
	negUy_1			= cnegUy_1;
	negUcap_1		= cnegUcap_1;
	negFy_1			= cnegFy_1;
	negFcap_1		= cnegFcap_1;
	negUGlobal_1	= cnegUGlobal_1;
	negFGlobal_1	= cnegFGlobal_1;
	negUres_1		= cnegUres_1;
	negFres_1		= cnegFres_1;
	negKp_1			= cnegKp_1;
	negKpc_1		= cnegKpc_1;
	K_unload		= cK_unload;
	K_reload		= cK_reload;
	Energy_Acc		= cEnergy_Acc;
	Energy_Diss		= cEnergy_Diss;
	posULocal_1		= cposULocal_1;
	posFLocal_1		= cposFLocal_1;
	negULocal_1		= cnegULocal_1;
	negFLocal_1		= cnegFLocal_1;
	Failure_Flag	= cFailure_Flag;
	// Excursion_Flag	= cExcursion_Flag;
	exBranch  		= cexBranch;
	Branch  			= cBranch;
	// TargetPeak_Flag = cTargetPeak_Flag;
	Yield_Flag  	= cYield_Flag;
	// Reversal_Flag	= cReversal_Flag;
	// u0				= cu0;
	return 0;
}
int IMKPeakOriented::revertToStart(void)
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
//////////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	betaS			= 0;
	betaC			= 0;
	betaK			= 0;
	betaA			= 0;
	u0				= 0;
	cdu_i_1			= 0;
	EtS				= LAMBDA_S *Fy_pos;
	EtC				= LAMBDA_C *Fy_pos;
	EtA				= LAMBDA_A *Fy_pos;
	EtK				= LAMBDA_K *Fy_pos;
	Uy_pos  		= Fy_pos / Ke;
	Ucap_pos		= Uy_pos + Up_pos;
	Fcap_pos		= FcapFy_pos*Fy_pos;
	Kp_pos 			= (Fcap_pos - Fy_pos) / Up_pos;
	Kpc_pos 		= Fcap_pos / Upc_pos;
	Uy_neg 			= Fy_neg / Ke;
	Ucap_neg		= Uy_neg + Up_neg;
	Fcap_neg		= FcapFy_neg*Fy_neg;
	Kp_neg 			= (Fcap_neg - Fy_neg) / Up_neg;
	Kpc_neg 		= Fcap_neg / Upc_neg;
	U				= cU			= 0;
	ui		 		= cui			= 0;
	fi		 		= cfi			= 0;
	ui_1 			= cui_1			= 0;
	fi_1 			= cfi_1			= 0;
	posULocal_1		= cposULocal_1	=  Uy_pos;
	posFLocal_1		= cposFLocal_1	=  Fy_pos;
	negULocal_1		= cnegULocal_1	= -Uy_neg;
	negFLocal_1		= cnegFLocal_1	= -Fy_neg;
	posUGlobal_1	= cposUGlobal_1	=  Uy_pos;
	posFGlobal_1	= cposFGlobal_1	=  Fy_pos;
	negUGlobal_1	= cnegUGlobal_1	= -Uy_neg;
	negFGlobal_1	= cnegFGlobal_1	= -Fy_neg;
	posUy_1 		= cposUy_1		=  Uy_pos;
	posFy_1 		= cposFy_1		=  Fy_pos;
	posKp_1 		= cposKp_1		=  Kp_pos;
	posKpc_1		= cposKpc_1		= -Kpc_pos;
	negUy_1 		= cnegUy_1		= -Uy_neg;
	negFy_1 		= cnegFy_1		= -Fy_neg;
	negKp_1 		= cnegKp_1		=  Kp_neg;
	negKpc_1		= cnegKpc_1		= -Kpc_neg;
	posUcap_1		= cposUcap_1	=  Ucap_pos;
	posFcap_1		= cposFcap_1	=  Fcap_pos;
	posFres_1		= cposFres_1	=  Fy_pos*ResF_pos;
	negUcap_1		= cnegUcap_1	= -Ucap_neg;
	negFcap_1		= cnegFcap_1	= -Fcap_neg;
	negFres_1		= cnegFres_1	= -Fy_neg*ResF_neg;
	posUres_1		= cposUres_1	=  (posFres_1 - posFcap_1) / posKpc_1 + posUcap_1;
	negUres_1		= cnegUres_1	=  (negFres_1 - negFcap_1) / negKpc_1 + negUcap_1;
	Energy_Acc 		= cEnergy_Acc	=  0;
	Energy_Diss		= cEnergy_Diss	=  0;
	Branch 			= cBranch 			= 0;
	Failure_Flag 	= cFailure_Flag	  	= false;
	Excursion_Flag 	= cExcursion_Flag 	= false;
	exBranch		= cexBranch 	= false;
	// TargetPeak_Flag	= cTargetPeak_Flag	= false;
	Yield_Flag		= cYield_Flag	  	= false;
	Reversal_Flag	= cReversal_Flag  	= false;
	K_unload		= cK_unload			= Ke;
	K_reload 	 	= cK_reload			= Ke;
	TangentK		= cTangentK			= Ke;
	return 0;
}
UniaxialMaterial *
IMKPeakOriented::getCopy(void)
{
	IMKPeakOriented *theCopy	= new IMKPeakOriented(this->getTag(), Ke,
		Uy_pos, Ucap_pos, Uu_pos, Fy_pos, FcapFy_pos, ResF_pos,
		Uy_neg, Ucap_neg, Uu_neg, Fy_neg, FcapFy_neg, ResF_neg,
		LAMBDA_S, LAMBDA_C, LAMBDA_A, LAMBDA_K, c_S, c_C, c_A, c_K, D_pos, D_neg);
	theCopy->U				= U;
	theCopy->cU				= cU;
	theCopy->TangentK		= TangentK;
	theCopy->ui				= ui;
	theCopy->fi				= fi;
	theCopy->ui_1			= ui_1;
	theCopy->fi_1			= fi_1;
	theCopy->du_i_1			= du_i_1;
	theCopy->posUy_1		= posUy_1;
	theCopy->posUcap_1		= posUcap_1;
	theCopy->posFy_1		= posFy_1;
	theCopy->posFcap_1		= posFcap_1;
	theCopy->posUGlobal_1	= posUGlobal_1;
	theCopy->posFGlobal_1	= posFGlobal_1;
	theCopy->posUres_1		= posUres_1;
	theCopy->posFres_1		= posFres_1;
	theCopy->posKp_1		= posKp_1;
	theCopy->posKpc_1		= posKpc_1;
	theCopy->negUy_1		= negUy_1;
	theCopy->negUcap_1		= negUcap_1;
	theCopy->negFy_1		= negFy_1;
	theCopy->negFcap_1		= negFcap_1;
	theCopy->negUGlobal_1	= negUGlobal_1;
	theCopy->negFGlobal_1	= negFGlobal_1;
	theCopy->negUres_1		= negUres_1;
	theCopy->negFres_1		= negFres_1;
	theCopy->negKp_1		= negKp_1;
	theCopy->negKpc_1		= negKpc_1;
	theCopy->K_unload		= K_unload;
	theCopy->Energy_Acc		= Energy_Acc;
	theCopy->Energy_Diss	= Energy_Diss;
	theCopy->u0				= u0;
	theCopy->posULocal_1	= posULocal_1;
	theCopy->posFLocal_1	= posFLocal_1;
	theCopy->negULocal_1	= negULocal_1;
	theCopy->negFLocal_1	= negFLocal_1;
	theCopy->Failure_Flag 	= Failure_Flag;
	theCopy->Excursion_Flag = Excursion_Flag;
	theCopy->exBranch 		= exBranch;
	// theCopy->TargetPeak_Flag= TargetPeak_Flag;
	theCopy->Yield_Flag 	= Yield_Flag;
	theCopy->Reversal_Flag	= Reversal_Flag;
	theCopy->K_reload		= K_reload;
	theCopy->cTangentK		= cTangentK;
	theCopy->cui			= cui;
	theCopy->cfi			= cfi;
	theCopy->cui_1			= cui_1;
	theCopy->cfi_1			= cfi_1;
	theCopy->cdu_i_1		= cdu_i_1;
	theCopy->cposUy_1		= cposUy_1;
	theCopy->cposUcap_1		= cposUcap_1;
	theCopy->cposFy_1		= cposFy_1;
	theCopy->cposFcap_1		= cposFcap_1;
	theCopy->cposUGlobal_1	= cposUGlobal_1;
	theCopy->cposFGlobal_1	= cposFGlobal_1;
	theCopy->cposUres_1		= cposUres_1;
	theCopy->cposFres_1		= cposFres_1;
	theCopy->cposKp_1		= cposKp_1;
	theCopy->cposKpc_1		= cposKpc_1;
	theCopy->cnegUy_1		= cnegUy_1;
	theCopy->cnegUcap_1		= cnegUcap_1;
	theCopy->cnegFy_1		= cnegFy_1;
	theCopy->cnegFcap_1		= cnegFcap_1;
	theCopy->cnegUGlobal_1	= cnegUGlobal_1;
	theCopy->cnegFGlobal_1	= cnegFGlobal_1;
	theCopy->cnegUres_1		= cnegUres_1;
	theCopy->cnegFres_1		= cnegFres_1;
	theCopy->cnegKp_1		= cnegKp_1;
	theCopy->cnegKpc_1		= cnegKpc_1;
	theCopy->cK_unload		= cK_unload;
	theCopy->cEnergy_Acc 	= cEnergy_Acc;
	theCopy->cEnergy_Diss	= cEnergy_Diss;
	theCopy->cu0			= cu0;
	theCopy->cposULocal_1	= cposULocal_1;
	theCopy->cposFLocal_1	= cposFLocal_1;
	theCopy->cnegULocal_1	= cnegULocal_1;
	theCopy->cnegFLocal_1	= cnegFLocal_1;
	theCopy->cFailure_Flag	= cFailure_Flag;
	theCopy->cExcursion_Flag= cExcursion_Flag;
	theCopy->cexBranch		= cexBranch;
	// theCopy->cTargetPeak_Flag= cTargetPeak_Flag;
	theCopy->cYield_Flag 	= cYield_Flag;
	theCopy->cReversal_Flag	= cReversal_Flag;
	theCopy->cK_reload		= cK_reload;
	return theCopy;
}
int IMKPeakOriented::sendSelf(int cTag, Channel &theChannel)
{
	int res		= 0;
	cout << " sendSelf" << endln;
	static Vector data(137);
	data(0)		= this->getTag();
	data(1)  	= Ke;
	data(2)  	= Uy_pos;
	data(3)  	= Ucap_pos;
	data(4)  	= Uu_pos;
	data(5)  	= Fy_pos;
	data(6)  	= FcapFy_pos;
	data(7)  	= ResF_pos;
	data(8)  	= Uy_neg;
	data(9)  	= Ucap_neg;
	data(10) 	= Uu_neg;
	data(11) 	= Fy_neg;
	data(12) 	= FcapFy_neg;
	data(13) 	= ResF_neg;
	data(14) 	= LAMBDA_S;
	data(15) 	= LAMBDA_C;
	data(16) 	= LAMBDA_A;
	data(17) 	= LAMBDA_K;
	data(18) 	= c_S;
	data(19) 	= c_C;
	data(20) 	= c_A;
	data(21) 	= c_K;
	data(22) 	= D_pos;
	data(23) 	= D_neg;
	data(24) 	= ui;
	data(25) 	= fi;
	data(26) 	= ui_1;
	data(27) 	= fi_1;
	data(28) 	= du_i_1;
	data(29) 	= posUy_1;
	data(30) 	= posUcap_1;
	data(31) 	= posFy_1;
	data(32) 	= posFcap_1;
	data(33) 	= posUGlobal_1;
	data(34) 	= posFGlobal_1;
	data(35) 	= posUres_1;
	data(36) 	= posFres_1;
	data(37) 	= posKp_1;
	data(38) 	= posKpc_1;
	data(39) 	= negUy_1;
	data(40) 	= negUcap_1;
	data(41) 	= negFy_1;
	data(42) 	= negFcap_1;
	data(43) 	= negUGlobal_1;
	data(44) 	= negFGlobal_1;
	data(45) 	= negUres_1;
	data(46) 	= negFres_1;
	data(47) 	= negKp_1;
	data(48) 	= negKpc_1;
	data(49) 	= K_unload;
	data(50) 	= Failure_Flag;
	data(51) 	= Excursion_Flag;
	data(52) 	= Branch;
	data(53) 	= exBranch;
	// data(54) 	= TargetPeak_Flag;
	data(55) 	= Yield_Flag;
	data(56) 	= Energy_Acc;
	data(57) 	= Energy_Diss;
	data(58) 	= u0;
	data(59) 	= du;
	data(60) 	= df;
	data(61) 	= FailS;
	data(62) 	= FailC;
	data(63) 	= FailA;
	data(64) 	= FailK;
	data(65)	= Ei;
	data(66)	= dEi;
	data(67)	= Epj;
	data(68)	= EpjK;
	data(69)	= EiK;
	data(70)	= c_S;
	data(71)	= c_C;
	data(72)	= c_A;
	data(73)	= c_K;
	data(74)	= EtS;
	data(75)	= EtC;
	data(76)	= EtA;
	data(77)	= EtK;
	data(78)	= betaS;
	data(79)	= betaC;
	data(80)	= betaA;
	data(81)	= betaK;
	data(82)	= sPCsp;
	data(83)	= sPCpcp;
	data(84)	= TangentK;
	data(85)	= Uy_pos;
	data(86)	= Ucap_pos;
	data(87)	= Fcap_pos;
	data(88)	= Kp_pos;
	data(89)	= Kpc_pos;
	data(90)	= Uy_neg;
	data(91)	= Ucap_neg;
	data(92)	= Fcap_neg;
	data(93)	= Kp_neg;
	data(94)	= Kpc_neg;
	data(95)	= cui;
	data(96)	= cfi;
	data(97)	= cui_1;
	data(98)	= cfi_1;
	data(99)	= cdu_i_1;
	data(100)	= cposUy_1;
	data(101)	= cposUcap_1;
	data(102)	= cposFy_1;
	data(103)	= cposFcap_1;
	data(104)	= cposUGlobal_1;
	data(105)	= cposFGlobal_1;
	data(106)	= cposUres_1;
	data(107)	= cposFres_1;
	data(108)	= cposKp_1;
	data(109)	= cposKpc_1;
	data(110)	= cnegUy_1;
	data(111)	= cnegUcap_1;
	data(112)	= cnegFy_1;
	data(113)	= cnegFcap_1;
	data(114)	= cnegUGlobal_1;
	data(115)	= cnegFGlobal_1;
	data(116)	= cnegUres_1;
	data(117)	= cnegFres_1;
	data(118)	= cnegKp_1;
	data(119)	= cnegKpc_1;
	data(120)	= cK_unload;
	data(121)	= cposULocal_1;
	data(122)	= cposFLocal_1;
	data(123)	= cnegULocal_1;
	data(124)	= cnegFLocal_1;
	data(125)	= cFailure_Flag;
	data(126)	= cExcursion_Flag;
	data(127)	= cexBranch;
	data(128)	= cBranch;
	// data(129)	= cTargetPeak_Flag;
	data(130)	= cYield_Flag;
	data(131)	= cK_reload;
	data(132)	= K_Local;
	data(133)	= K_Global;
	// data(134)	= K_check;
	data(135)	= cReversal_Flag;
	data(136)	= Reversal_Flag;
	res			= theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "IMKPeakOriented::sendSelf() - failed to send data\n";
	return res;
}
int IMKPeakOriented::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res	= 0;
	static Vector data(137);
	res	= theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		opserr << "IMKPeakOriented::recvSelf() - failed to receive data\n";
		this->setTag(0);
	} else {
		cout << " recvSelf" << endln;
		this->setTag( (int)data(0) );
		Ke				= data(1);
		Up_pos			= data(2);
		Upc_pos			= data(3);
		Uu_pos			= data(4);
		Fy_pos			= data(5);
		FcapFy_pos		= data(6);
		ResF_pos		= data(7);
		Up_neg			= data(8);
		Upc_neg			= data(9);
		Uu_neg			= data(10);
		Fy_neg			= data(11);
		FcapFy_neg		= data(12);
		ResF_neg		= data(13);
		LAMBDA_S		= data(14);
		LAMBDA_C		= data(15);
		LAMBDA_A		= data(16);
		LAMBDA_K		= data(17);
		c_S				= data(18);
		c_C				= data(19);
		c_A				= data(20);
		c_K				= data(21);
		D_pos			= data(22);
		D_neg			= data(23);
		ui				= data(24);
		fi				= data(25);
		ui_1			= data(26);
		fi_1			= data(27);
		du_i_1			= data(28);
		posUy_1			= data(29);
		posUcap_1		= data(30);
		posFy_1			= data(31);
		posFcap_1		= data(32);
		posUGlobal_1	= data(33);
		posFGlobal_1	= data(34);
		posUres_1		= data(35);
		posFres_1		= data(36);
		posKp_1			= data(37);
		posKpc_1		= data(38);
		negUy_1			= data(39);
		negUcap_1		= data(40);
		negFy_1			= data(41);
		negFcap_1		= data(42);
		negUGlobal_1	= data(43);
		negFGlobal_1	= data(44);
		negUres_1		= data(45);
		negFres_1		= data(46);
		negKp_1			= data(47);
		negKpc_1		= data(48);
		Failure_Flag	= data(49);
		Excursion_Flag	= data(50);
		exBranch		= data(51);
		Branch			= data(52);
		// TargetPeak_Flag	= data(53);
		Yield_Flag		= data(54);
		K_unload		= data(55);
		Energy_Acc		= data(56);
		Energy_Diss		= data(57);
		u0				= data(58);
		du				= data(59);
		df				= data(60);
		FailS			= data(61);
		FailC			= data(62);
		FailA			= data(63);
		FailK			= data(64);
		Ei				= data(65);
		dEi				= data(66);
		Epj				= data(67);
		EpjK			= data(68);
		EiK				= data(79);
		c_S				= data(70);
		c_C				= data(71);
		c_A				= data(72);
		c_K				= data(73);
		EtS				= data(74);
		EtC				= data(75);
		EtA				= data(76);
		EtK				= data(77);
		betaS			= data(78);
		betaC			= data(79);
		betaA			= data(80);
		betaK			= data(81);
		sPCsp			= data(82);
		sPCpcp			= data(83);
		TangentK		= data(84);
		Uy_pos			= data(85);
		Ucap_pos		= data(86);
		Fcap_pos		= data(87);
		Kp_pos			= data(88);
		Kpc_pos			= data(89);
		Uy_neg			= data(90);
		Ucap_neg		= data(91);
		Fcap_neg		= data(92);
		Kp_neg			= data(93);
		Kpc_neg			= data(94);
		cui				= data(95);
		cfi				= data(96);
		cui_1			= data(97);
		cfi_1			= data(98);
		cdu_i_1			= data(99);
		cposUy_1		= data(100);
		cposUcap_1		= data(101);
		cposFy_1		= data(102);
		cposFcap_1		= data(103);
		cposUGlobal_1	= data(104);
		cposFGlobal_1	= data(105);
		cposUres_1		= data(106);
		cposFres_1		= data(107);
		cposKp_1		= data(108);
		cposKpc_1		= data(109);
		cnegUy_1		= data(110);
		cnegUcap_1		= data(111);
		cnegFy_1		= data(112);
		cnegFcap_1		= data(113);
		cnegUGlobal_1	= data(114);
		cnegFGlobal_1	= data(115);
		cnegUres_1		= data(116);
		cnegFres_1		= data(117);
		cnegKp_1		= data(118);
		cnegKpc_1		= data(119);
		cK_unload		= data(120);
		cposULocal_1	= data(121);
		cposFLocal_1	= data(122);
		cnegULocal_1	= data(123);
		cnegFLocal_1	= data(124);
		cFailure_Flag	= data(125);
		cExcursion_Flag	= data(126);
		cexBranch		= data(127);
		cBranch			= data(128);
		// cTargetPeak_Flag= data(129);
		cYield_Flag   	= data(130);
		cK_reload      	= data(131);
		K_Local			= data(132);
		K_Global		= data(133);
		// K_check			= data(134);
		cReversal_Flag	= data(135);
		Reversal_Flag	= data(136);
	}
	return res;
}
void IMKPeakOriented::Print(OPS_Stream &s, int flag)
{
	cout << "IMKPeakOriented tag: " << this->getTag() << endln;
}
