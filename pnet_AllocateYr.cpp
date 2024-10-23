#include "pnet_model.h"

void pnet_model::AllocateYr(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode)
{
	//
	//Annual C allocation for the PnET ecosystem model.
	//

	int i;
	double EnvMaxFol, SppMaxFol, FolRegen, BiomLossFrac, TotalC, nr;
	double PotLightEff;  // potential light for growth
	double budc, woodc, rootc, budn; //param use to calculate potental plant pool
	double folnconnew;
	double woodmass, woodmassN;

	//Check for a disturbance year
	FolRegen = 100;
	BiomLossFrac = 0;

	for (i = 1; i <= site->distyrs; i++)
	{
		if (clim->year[rstep] == site->distyear[i])
		{
			BiomLossFrac = site->distintensity[i];
			veg->FolMassMax = veg->FolMassMax * (1 - BiomLossFrac);
			if (veg->FolMassMax < site->folregen[i])veg->FolMassMax = site->folregen[i];
			break;
		}
	}

	//Agriculture  Linghui Meng 05062020 move agriculture to year end calculation
	if (clim->year[rstep] >= site->agstart && clim->year[rstep] < site->agstop)
	{

		share->WoodMass = share->WoodMass * (1 - site->agrem * (1 - site->agrem));
		share->WoodMassN = share->WoodMassN * (1 - site->agrem * (1 - site->agrem));

		share->DeadWoodM = share->DeadWoodM * (1 - site->agrem * (1 - site->agrem)); //Linghui 03252020
		share->DeadWoodN = share->DeadWoodN * (1 - site->agrem * (1 - site->agrem)); //Linghui 03252020

		share->PlantC = share->PlantC * (1 - site->agrem * (1 - site->agrem));
		share->PlantN = share->PlantN * (1 - site->agrem * (1 - site->agrem));

		share->HON = share->HON + 1.5; //Linghui Meng 11262020 assume N fertilization
	}

	if (clim->year[rstep] >= site->disestart && clim->year[rstep] < site->disestop)
	{

		woodmass = share->WoodMass * site->distintens;
		woodmassN = share->WoodMassN * site->distintens;

		share->WoodMass = share->WoodMass - woodmass; 
		share->WoodMassN = share->WoodMassN - woodmassN;

		share->DeadWoodM = share->DeadWoodM + woodmass; //Linghui 03252020
		share->DeadWoodN = share->DeadWoodN + woodmassN; //Linghui 03252020

//		share->PlantC = share->PlantC * (1 - site->distintens);
//		share->PlantN = share->PlantN * (1 - site->distintens);

	}

	//if the decidous tree foliar mass greater than 0 at the end of year, make it to 0 and relocate C and N
	if ((veg->FolReten == 1)&&(share->FolMass>0))
	{
		double FolNLoss, Retrans, FolLitN;
		share->FolLitM = share->FolMass;
		FolNLoss = share->FolLitM * (veg->FolNCon / 100);
		Retrans = FolNLoss * veg->FolNRetrans;
		share->PlantN = share->PlantN + Retrans;
		FolLitN = FolNLoss - Retrans;

		share->TotalLitterM = share->TotalLitterM + share->FolLitM;
		share->TotalLitterN = share->TotalLitterN + FolLitN;

		share->HOM = share->HOM + share->TotalLitterM / share->dayspan;
		share->HON = share->HON + share->TotalLitterN / share->dayspan;


		share->TotalLitterM = share->TotalLitterM - share->TotalLitterM / share->dayspan;
		share->TotalLitterN = share->TotalLitterN - share->TotalLitterN / share->dayspan;

		share->FolMass = 0;
	}

		share->NPPFolYr = share->FolProdCYr / veg->CFracBiomass;
		share->NPPWoodYr = share->WoodProdCYr / veg->CFracBiomass;
		share->NPPRootYr = share->RootProdCYr / veg->CFracBiomass;

	//reallocate all C for next year allocation LM 20240802
	share->PlantC = share->PlantC + share->WoodC + share->RootC;
	share->WoodC = 0;
	share->RootC = 0;

	if (share->DwaterIx > 0)
	{
		share->AvgDWater = share->Dwatertot / share->DwaterIx;
		//	share->AvgDWater = 1.0;
	}
	else
	{
		share->AvgDWater = 1.0;
	}


	if (share->PosCBalMassIx > 0)
	{
		share->avgPCBM = (share->PosCBalMassTot / share->PosCBalMassIx);
	}
	else
	{
		share->avgPCBM = share->FolMass;
	}


	share->CanopyDO3Avg = share->CanopyDO3Tot / share->CanopyDO3TotDay;


	if (share->LightEffCBalIx > 0)
	{
		PotLightEff = (share->LightEffCBalTot / share->LightEffCBalIx);
	}
	else
	{
		PotLightEff = 1.0;
	}

	share->LightEffMin = PotLightEff;

	EnvMaxFol = (share->AvgDWater * share->avgPCBM) * (1.0 + (veg->FolRelGrowMax * share->LightEffMin)); //convert foliar mass to budc
	SppMaxFol = share->avgPCBM * (1.0 + (veg->FolRelGrowMax * share->LightEffMin));

//	share->EnvMaxFol = EnvMaxFol;
//	share->SppMaxFol = SppMaxFol;


	if (EnvMaxFol < SppMaxFol)
	{
		veg->FolMassMax = EnvMaxFol;
	}
	else
	{
		veg->FolMassMax = SppMaxFol;
	}

	veg->FolMassMin = (veg->FolMassMax - veg->FolMassMax * (1.0 / veg->FolReten));

	//balance budc, woodc, and plantC, replace old algorithm below LM 812024

	share->BudC = ((veg->FolMassMax - share->FolMass) * veg->CFracBiomass);
	share->WoodC = (1.0 - veg->PlantCReserveFrac) * share->PlantC;
	
	if (share->BudC < 0)
	{
		share->BudC = 0;
		if (veg->FolReten > 1) share->BudC = share->FolMass * 1 / (veg->FolReten - 1) * veg->CFracBiomass * 0.5;  // evergreen has half budc

	}

	TotalC = share->WoodC + share->BudC;

	if (share->WoodC < (veg->MinWoodFolRatio * share->BudC))
	{
		
		share->WoodC = TotalC * (veg->MinWoodFolRatio / (1.0 + veg->MinWoodFolRatio));
		share->BudC = TotalC - share->WoodC;
		veg->FolMassMax = share->FolMass + (share->BudC / veg->CFracBiomass);
		veg->FolMassMin = (veg->FolMassMax - veg->FolMassMax * (1.0 / veg->FolReten));
	}



	if (TotalC > share->PlantC) {
		share->BudC = share->BudC * share->PlantC / TotalC;
		share->WoodC = share->WoodC * share->PlantC / TotalC;
		
	}

	share->PlantC = share->PlantC - share->WoodC - share->BudC;



	/*
	share->BudC = ((veg->FolMassMax - share->FolMass) * veg->CFracBiomass);
	//	share->BudC = 308. *veg->CFracBiomass;

	if (share->BudC < 0)
	{
		share->BudC = 0;
		if (veg->FolReten>1) share->BudC = share->FolMass * 1 / (veg->FolReten - 1)*veg->CFracBiomass*0.5;  // evergreen has half budc

	}


	share->PlantC = share->PlantC - share->BudC;
	share->WoodC = (1.0 - veg->PlantCReserveFrac) * share->PlantC;




	if (share->WoodC < (veg->MinWoodFolRatio * share->BudC))
	{
		TotalC = share->WoodC + share->BudC;
		share->WoodC = TotalC * (veg->MinWoodFolRatio / (1.0 + veg->MinWoodFolRatio));
		share->BudC = TotalC - share->WoodC;
		veg->FolMassMax = share->FolMass + (share->BudC / veg->CFracBiomass);
		veg->FolMassMin = (veg->FolMassMax - veg->FolMassMax * (1.0 / veg->FolReten));
	}


	share->PlantC = share->PlantC - share->WoodC;
	*/

	//NEP calculation for PnET-II
	share->NEP = share->TotPsn - share->WoodMRespYr - share->WoodGRespYr - share->FolGRespYr - share->SoilRespYr;
	// save current foliar N for output
	share->FolNConOld = veg->FolNCon;

	// PnET-CN Only -----------------------------------------------------------------
	if (CN_Mode == 1)
	{
		
		if (share->PlantN > veg->MaxNStore)
		{
			share->NH4 = share->NH4 + (share->PlantN - veg->MaxNStore);  // ZZX revised
			share->PlantN = veg->MaxNStore;
		}


		share->NRatio = 1 + (share->PlantN / veg->MaxNStore) * veg->FolNConRange;

		if (share->NRatio < 1)
		{
			share->NRatio = 1;
		}

		if (share->NRatio >(1 + veg->FolNConRange))
		{
			share->NRatio = 1 + veg->FolNConRange;
		}

		if (share->NRatio < 1)
		{
			share->NRatioNit = 0;
		}
		else
		{

			nr = share->NRatio - 1 - (veg->FolNConRange / 3);
			if (nr < 0)nr = 0;

			share->NRatioNit = (nr / (0.6667 * veg->FolNConRange))*(nr / (0.6667 * veg->FolNConRange));
			if (share->NRatioNit > 1) share->NRatioNit = 1;

		}

		//adjust budN make sure there is enough N for boost LM 812024
		budn = (share->BudC / veg->CFracBiomass) * veg->FLPctN * (1 / (1 - veg->FolNRetrans)) * share->NRatio;

		if (budn > share->PlantN)
		{
			budn = share->PlantN;

			share->BudC = budn * veg->CFracBiomass / veg->FLPctN / (1 / (1 - veg->FolNRetrans)) / share->NRatio;

		}

		//set up foliar N concentration for the next year in PnET-CN LM 962024
		if (clim->timestep == 0) {
			folnconnew = (share->FolMass * (veg->FolNCon / 100) + share->BudN) / (share->FolMass + (share->BudC / veg->CFracBiomass)) * 100;
			veg->FolNCon = folnconnew;
		}

		share->BudN = budn;
		share->PlantN = share->PlantN - share->BudN;


//		share->RootNSinkEff = sqrt(1 - (share->PlantN / veg->MaxNStore));   // this is based on monthly rate

		//Annual total variables for PnET-CN
		share->NEP = share->TotPsn - share->SoilDecRespYr - share->WoodDecRespYr - share->WoodMRespYr
			- share->WoodGRespYr - share->FolGRespYr - share->RootMRespYr - share->RootGRespYr;
//		share->FolN = (share->FolMass * veg->FolNCon / 100); LM 07312024
		share->FolC = share->FolMass * veg->CFracBiomass;
		share->SppMaxFol = share->FolN + share->WoodMassN + share->RootMassN + share->HON
			+ share->NH4 + share->NO3 + share->BudN + share->DeadWoodN + share->PlantN+share->TotalLitterN;
		share->EnvMaxFol = share->FolN + share->WoodMassN + share->RootMassN + share->PlantN;
 		share->TotalM = (share->BudC / veg->CFracBiomass) + share->FolMass + (share->WoodMass + share->WoodC / veg->CFracBiomass)
			+ share->RootMass + share->DeadWoodM + share->HOM + (share->PlantC / veg->CFracBiomass) + (share->RootC / veg->CFracBiomass);
	}
	
	//calculate ISA in each year Linghui Meng 20220807
//	if (modelmode == 7 && clim->year[rstep] >= site->yrurstart) {share->In = 1;}
	if (clim->year[rstep] >= site->yrurstart&&share->ISA>0) { share->In = 1; } //turn on ISA function after the year yrurstart

	//calculate ISA for next year
	if (share->In>0&&clim->year[rstep] < site->yrurstart2) {

		share->ISA = share->ISA + (share->ISA1 / 120); //calculated ISA in during stage 1
	}
	else
	{
		share->ISA = share->ISA + (share->ISA2 / 30);

	}
	
	//Agriculture  Linghui Meng 05062020 move agriculture to year end calculation
	if (clim->year[rstep] >= site->agstart && clim->year[rstep] < site->agstop)
	{
		share->WoodMass = share->WoodMass * (1 - site->agrem * (1 - site->agrem));
		share->WoodMassN = share->WoodMassN * (1 - site->agrem * (1 - site->agrem));
		share->DeadWoodM = share->DeadWoodM * (1 - site->agrem * (1 - site->agrem)); //Linghui 03252020
		share->DeadWoodN = share->DeadWoodN * (1 - site->agrem * (1 - site->agrem)); //Linghui 03252020
		share->PlantC = share->PlantC * (1 - site->agrem * (1 - site->agrem));
		share->PlantN = share->PlantN * (1 - site->agrem * (1 - site->agrem));
		share->HON = share->HON + 1.5; //Linghui Meng 11262020 assume N addition
	}

}
