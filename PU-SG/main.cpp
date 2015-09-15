#include "Limits.h"
#include "Laplace.h"
#include "Solver.h"
#include "BaseLine.h"
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;

int main()
{
	//		city:		Houston, Chicago, LA, NYC
	int lunum[] =	{ 36199, 89796, 124040, 338144 };
	int levnum[] = { 16694, 36009, 54538, 108170 };
	int lgnum[] = { 3813, 6804, 13441, 19918};
	int luenum[] = { 97216, 209203, 327711, 685029 };
	int lugnum[] = { 60133, 97434, 214422, 427159 };

	// Houston
	{
		// model:		PU-GS, PU-G, PU-S, PU-0
		double alphax[4] = { 0.97, 0.90, 0.91, 0.93};
		double lambda[4] = { 0.04, 0.02, 0.07, 0.01 };
		double rho_on[4] = { 0.9, 0, 0.5, 0};
		double rho_off[4] = { 0.7, 0, 0.1, 0};

		// set up Laplace and sample. sample is contained in the lap_off processing.
		string f_on("Hou_ug.txt");
		string f_off("Hou_ue.txt");
		Soc::SetUpSocialLap(f_on, f_off, lunum[0], levnum[0], lgnum[0], luenum[0], lugnum[0]);

		string u_loc("Hou_ll_user.txt");
		string e_loc("Hou_ll_event.txt");
		Soc::ReadInLocation(u_loc, e_loc);

		// solver initial
		Mysol::BoolMatrix();

		// solver by model order
		bool d_flag = true;
		int disbest = 50;
		double AUC_model[4] = { 0 };
		double AUC_best = 0;
		double pat1 = 0, pat3 = 0, pat5 = 0, MAP = 0;
		for (int model = 0; model < 4; model++)
		{
			Soc::LapAdd(rho_on[model], rho_off[model]);
			if (model > 1)
				d_flag = false;

			//if (model == 0)
			//{
			//	int dis[6] = { 10, 20, 50, 100, 200, 500 };
			//	for (int i0 = 0; i0 < 6; i0++)
			//	{
			//		Mysol::Initial(dis[i0], alphax[model], lambda[model], d_flag);
			//		bool con_flag = Mysol::CoordianteDecent();
			//		if (!con_flag)
			//			cout << "model = " << model << ", dis = " << dis[i0] << endl;
			//		AUC_dis[i0] = Mysol::EvalutionAUC();
			//		Mysol::FreeWandH();

			//		if (AUC_best < AUC_dis[i0])
			//		{
			//			AUC_best = AUC_dis[i0];
			//			disbest = dis[i0];
			//		}
			//	}
			//	AUC_model[model] = AUC_best;
			//}
			// else
			{
				Mysol::Initial(disbest, alphax[model], lambda[model], d_flag);
				bool con_flag = Mysol::CCD_PlusPlus();
				if (!con_flag)
					cout << "model = " << model << ", dis = " << disbest << endl;
				AUC_model[model] = Mysol::EvalutionAUC();
				if (model == 0)
				{
					MAP = Mysol::EvalutionMAP(&pat1, &pat3, &pat5);
				}
				Mysol::FreeWandH(model);
			}
		}

		Mysol::FreeBoolMatrix();

		ofstream Hou("Hou.txt");
		for (int i0 = 0; i0 < 4; i0++)
			Hou << AUC_model[i0] << endl;
		Hou << "----------------------\n";
		Hou << "PU-SG\n";
		Hou << "p@1:\t" << pat1 << endl;
		Hou << "p@3:\t" << pat3 << endl;
		Hou << "p@5:\t" << pat5 << endl;
		Hou << "MAP:\t" << MAP << endl;
		Hou << "----------------------\n";
		MAP = Mybase::EvalutionRandom(&pat1, &pat3, &pat5);
		Hou << "Random\n";
		Hou << "p@1:\t" << pat1 << endl;
		Hou << "p@3:\t" << pat3 << endl;
		Hou << "p@5:\t" << pat5 << endl;
		Hou << "MAP:\t" << MAP << endl;
		MAP = Mybase::EvalutionMostPo(&pat1, &pat3, &pat5);
		Hou << "Most Popular\n";
		Hou << "p@1:\t" << pat1 << endl;
		Hou << "p@3:\t" << pat3 << endl;
		Hou << "p@5:\t" << pat5 << endl;
		Hou << "MAP:\t" << MAP << endl;
		MAP = Mybase::EvalutionLocAwa(&pat1, &pat3, &pat5);
		Hou << "Location-Aware\n";
		Hou << "p@1:\t" << pat1 << endl;
		Hou << "p@3:\t" << pat3 << endl;
		Hou << "p@5:\t" << pat5 << endl;
		Hou << "MAP:\t" << MAP << endl;
		MAP = Mybase::EvalutionUserKNN(&pat1, &pat3, &pat5);
		Hou << "user-knn\n";
		Hou << "p@1:\t" << pat1 << endl;
		Hou << "p@3:\t" << pat3 << endl;
		Hou << "p@5:\t" << pat5 << endl;
		Hou << "MAP:\t" << MAP << endl;
		Hou.close();

	}

	return 0;
}
