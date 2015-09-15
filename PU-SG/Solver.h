#ifndef SOLVER_H
#define SOLVER_H

class Mysol
{
private: 
	static int dis;
	static int lunum;
	static int levnum;

	static double alphax;
	static double lambda;
	//static double rho_on;
	//static double rho_off;
	static double accu;
	static int itmax;
	static double eta;
	static bool d_flag;

	static bool** uematrix;
	static double** Matrix_W;
	static double** Matrix_H;

public:
	static void BoolMatrix();
	static void FreeBoolMatrix();
	static void InitialWandH();
	static void FreeWandH(int model);
	static void Initial(int dis_in, double alphax_in, double lambda_in, bool flag_in);

	//---------------------------------
	static bool CoordianteDecent(); // whether is convergence
	static void UpdateWi(int wi);
	static void UpdateHi(int hi);
	static void ArryConst(double* arryin, double para, int arry_size);
	static double DotW_H(int wi, int hi, int a_size);
	static double FrobeniusNorm(double** matrix, int row, int column); 

	static void UpdateWi_v2(int wi);
	static void UpdateHi_v2(int hi);
	static double LapsAndW(double *laps_row, int w_j, int a_size);
	//--------------------------------
	static double EvalutionAUC();
	static double EvalutionAP(int u_node, double* pat1, double* pat3, double* pat5);
	static double EvalutionMAP(double* pat1, double* pat3, double* pat5);

	//-------------------------------------
	static bool CCD_PlusPlus();
	
};

#endif