#include "pnet_model.h"

void pnet_model::AllocateMo(veg_struct* veg, share_struct* share, int rstep, int CN_Mode)
{
	//
	// C allocation for the PnET ecosystem model.
	//


	double WoodMRespMo, GDDWoodEff, delGDDWoodEff, WoodProdCMo, WoodGRespMo, TMult, RootCAdd, RootAllocCMo, RootProdCMo, RootMRespMo, RootGRespMo;
	double rpotential, wpotential;
	double budc, budn, woodn, rootn, totaln;
	double folnconnew;



	share->PlantC = share->PlantC + share->NetPsnMo - share->FolGRespMo;
//	WoodMRespMo = share->CanopyGrossPsnActMo * veg->WoodMRespA;
	WoodMRespMo = share->GrsPsnMo * veg->WoodMRespA;
	share->WoodMRespYr = share->WoodMRespYr + WoodMRespMo;
	share->FolProdCYr = share->FolProdCYr + share->FolProdCMo;
	share->FolGRespYr = share->FolGRespYr + share->FolGRespMo;

	
	if (share->FolGDD == 1) {

		/*
		share->BudN = (share->BudC / veg->CFracBiomass) * veg->FLPctN * (1 / (1 - veg->FolNRetrans)) * share->NRatio;//ZZX


		
		if (share->BudN > share->PlantN)
		{
			if (share->PlantN < 0)
			{
				share->BudC = share->BudC * 0.1;
				share->BudN = share->BudN * 0.1;
			}
			else
			{
//				share->BudC = share->BudC * (share->PlantN / share->BudN);
				budc = share->BudN * veg->CFracBiomass / veg->FLPctN / (1 / (1 - veg->FolNRetrans))/ share->NRatio; 
				if (budc < share->BudC) { share->BudC = budc; }//LM 20240730 adjust budc based on available budn
				share->BudN = share->BudN * (share->PlantN / share->BudN);


			}
		}
		*/

		folnconnew = (share->FolMass * (veg->FolNCon / 100) + share->BudN) / (share->FolMass + (share->BudC / veg->CFracBiomass)) * 100;
		veg->FolNCon = folnconnew;

		share->FolGDD = 0; //LM 812024
	}
	

	//calculate root carbon addition
	TMult = (exp(0.1 * (share->Tave - 7.1)) * 0.68) * 1.0;
	RootCAdd = veg->RootAllocA * (share->dayspan / 365.0) + veg->RootAllocB * share->FolProdCMo;
	
	
	//if plantc is not enough, stop allocate carbon
	if (share->PlantC < RootCAdd) {
		if (share->PlantC < 0) {
			RootCAdd = 0;
		}
		else
		{
			RootCAdd = share->PlantC;
		}
	}
	
	share->RootC = share->RootC + RootCAdd;
	share->PlantC = share->PlantC - RootCAdd;
	//calculate root allocation
	RootAllocCMo = (share->dayspan / 365.0) * TMult;   // modified from matlab version
	if (RootAllocCMo > 1.0) RootAllocCMo = 1.0;
	RootAllocCMo = RootAllocCMo * share->RootC;
	RootProdCMo = RootAllocCMo / (1.0 + veg->RootMRespFrac + veg->GRespFrac);



	//calculate wood production
	if (share->GDDTot >= veg->GDDWoodStart)
	{
		GDDWoodEff = (share->GDDTot - veg->GDDWoodStart) / (veg->GDDWoodEnd - veg->GDDWoodStart);
		if (GDDWoodEff > 1.0)GDDWoodEff = 1;
		if (GDDWoodEff < 0)GDDWoodEff = 0;

		delGDDWoodEff = GDDWoodEff - share->OldGDDWoodEff;
		WoodProdCMo = share->WoodC * delGDDWoodEff;
		share->OldGDDWoodEff = GDDWoodEff;

	}
	else
	{
		WoodProdCMo = 0;
		WoodGRespMo = 0;
	}


	if (CN_Mode == 1)
	{
		//theoritical pontential root can growth with N limititation 
		rootn = RootProdCMo / veg->CFracBiomass* veg->RLPctN* share->NRatio;
		woodn = WoodProdCMo / veg->CFracBiomass* veg->WLPctN* share->NRatio;
		totaln = rootn + woodn;

		//plantN is less than N requirement to build up biomass, reduce biomass LM 822024 
		if (share->PlantN < totaln && totaln >0) {

			RootProdCMo = RootProdCMo * (share->PlantN / totaln);
			WoodProdCMo = WoodProdCMo * (share->PlantN / totaln);

			rootn = rootn * (share->PlantN / totaln);
			woodn = woodn * (share->PlantN / totaln);
		}

		//allocate wood and root
		WoodProdCMo = woodn * veg->CFracBiomass / veg->WLPctN / share->NRatio;
		RootProdCMo = rootn * veg->CFracBiomass / veg->WLPctN / share->NRatio;
	}



	share->WoodMass = share->WoodMass + (WoodProdCMo / veg->CFracBiomass);
	share->RootMass = share->RootMass + (RootProdCMo / veg->CFracBiomass);

	//substract cost from N pool
	share->WoodMassN = share->WoodMassN + ((WoodProdCMo / veg->CFracBiomass) * veg->WLPctN * share->NRatio);
	share->RootMassN = share->RootMassN + ((RootProdCMo / veg->CFracBiomass) * veg->RLPctN * share->NRatio);
	share->PlantN = share->PlantN - ((WoodProdCMo / veg->CFracBiomass) * veg->WLPctN * share->NRatio) - ((RootProdCMo / veg->CFracBiomass) * veg->RLPctN * share->NRatio);

//	share->WoodMassN = share->WoodMassN + woodn;
//	share->RootMassN = share->RootMassN + rootn;
//	share->PlantN = share->PlantN - woodn - rootn;


	if (share->PlantN < 0) { share->PlantN = 0; }

	//substract cost from C pool
	share->RootC = share->RootC - RootAllocCMo;
	RootMRespMo = RootProdCMo * veg->RootMRespFrac;
	RootGRespMo = RootProdCMo * veg->GRespFrac;
	WoodGRespMo = WoodProdCMo * veg->GRespFrac;

	share->WoodProdCYr = share->WoodProdCYr + WoodProdCMo;
	share->WoodGRespYr = share->WoodGRespYr + WoodGRespMo;


	share->WoodProdCMo = WoodProdCMo;
	share->PlantC = share->PlantC - WoodMRespMo - WoodGRespMo;
	share->WoodC = share->WoodC - WoodProdCMo;

	share->RootProdCYr = share->RootProdCYr + RootProdCMo;
	share->RootMRespYr = share->RootMRespYr + RootMRespMo;
	share->RootGRespYr = share->RootGRespYr + RootGRespMo;




	share->RootProdCMo = RootProdCMo;

	//PnET-CN Only -----------------------------------------------------------------
	if (CN_Mode == 1)
	{
		share->NetCBal = share->NetCBal + share->NetPsnMo - WoodMRespMo - WoodGRespMo - share->FolGRespMo - RootMRespMo - RootGRespMo;
		// needs -share->SoilDecResp - share->WoodDecResp, and will be updated in the respective routine.
		share->SppMaxFol = share->BudN + share->PlantN;
		share->EnvMaxFol = share->WoodMassN + share->RootMassN;
	}
	else
	{
		//------------------------------------------------------------------------------
		share->NetCBal = share->NetCBal + share->NetPsnMo - WoodMRespMo - WoodGRespMo - share->FolGRespMo - share->SoilRespMo;
	}

}

