#ifndef BASELINE_H
#define BASELINE_H

#include "Laplace.h"
#include <vector>

using std::vector;

class Mybase
{
private:
	static int lunum;
	static int levnum;

	static int dis;
	static int ratio;
	static double** Matrix_W;
	static double** Matrix_H;

	static vector<uxpair> vu_p;
	static vector<int> vindexp;

	static vector<uxpair> vu_n;
	static vector<int> vindexn;

public:
	static double EvalutionRandom(double* pat1, double* pat3, double* pat5);
	static double EvalutionMostPo(double* pat1, double* pat3, double* pat5);
	static double EvalutionLocAwa(double* pat1, double* pat3, double* pat5);
	static double EvalutionUserKNN(double* pat1, double* pat3, double* pat5);

	//------------- BPR-MF -------------
	static void MFInitial(int in_lunum, int in_levnum, int in_dis, int in_ratio = 10);
	static void MFFreeWandH();
	static void MFSetUpTrainPair();
	static bool MFTraining(int itermax = 1000, double drate = 0.01, double dlambda = 0.1, double daccu = 0.001);
	static double MFDot(const double *a, const double *b, int num);
	static double MFEvalution(double* pat1, double* pat3, double* pat5);

};

#endif