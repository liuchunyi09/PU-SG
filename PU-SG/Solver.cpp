#include "Laplace.h"
#include "Solver.h"
#include "Limits.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

int Mysol::dis    = 0;
int Mysol::lunum  = 0;
int Mysol::levnum = 0;

double Mysol::alphax    = 0.0;
double Mysol::lambda    = 0.0;
//double Mysol::rho_on  = 0.0;
//double Mysol::rho_off = 0.0;
double Mysol::accu      = 0.0;
int Mysol::itmax        = 0;
double Mysol::eta       = 0.0;
bool Mysol::d_flag      = true;

bool** Mysol::uematrix   = nullptr;

double** Mysol::Matrix_W = nullptr;
double** Mysol::Matrix_H = nullptr;

//==================================================

void Mysol::BoolMatrix()
{
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();
	int train_num = Soc::GetUETrainnum();

	uematrix = new bool*[lunum];
	for (int i0 = 0; i0 < lunum; i0++)
	{
		uematrix[i0] = new bool[levnum];
		memset(uematrix[i0], false, sizeof(bool)* levnum);
	}

	for (int i0 = 0; i0 < train_num; i0++)
	{
		uxpair temp = Soc::GetPair(i0);
		uematrix[temp.u][temp.x] = true;
	}
}

void Mysol::FreeBoolMatrix()
{
	for (int i0 = 0; i0 < lunum; i0++)
		delete[] uematrix[i0];
	delete[] uematrix;
}

void Mysol::InitialWandH()
{
	Matrix_W = new double*[lunum];
	for (int i0 = 0; i0 < lunum; i0++)
	{
		Matrix_W[i0] = new double[dis];
		memset(Matrix_W[i0], 0, sizeof(double)* dis);
	}

	Matrix_H = new double*[dis];
	for (int i0 = 0; i0 < dis; i0++)
	{
		Matrix_H[i0] = new double[levnum];
		memset(Matrix_H[i0], 0, sizeof(double)* levnum);
	}
}

void Mysol::FreeWandH(int model)
{
	if (model == 0)
	{
		ofstream aa("W.txt");
		for (int i0 = 0; i0 < lunum; i0++)
		{
			for (int i1 = 0; i1 < dis; i1++)
				aa << Matrix_W[i0][i1] << '\t';
			aa << endl;
		}
		aa.close();
	}

	for (int i0 = 0; i0 < lunum; i0++)
		delete[] Matrix_W[i0];
	delete[] Matrix_W;

	if (model == 0)
	{
		ofstream bb("H.txt");
		for (int i0 = 0; i0 < dis; i0++)
		{
			for (int i1 = 0; i1 < levnum; i1++)
				bb << Matrix_W[i0][i1] << '\t';
			bb << endl;
		}
		bb.close();
	}

	for (int i1 = 0; i1 < dis; i1++)
		delete[] Matrix_H[i1];
	delete[] Matrix_H;
}

void Mysol::Initial(int dis_in, double alphax_in, double lambda_in, bool flag_in)
{
	dis = dis_in;
	alphax = alphax_in;
	lambda = lambda_in;
	d_flag = flag_in;
	accu = ACCU;
	eta = ETA;
	itmax = ITMAX;
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();

	// set up ue matrix
	// BoolMatrix();

	// set up all-0 W and H
	InitialWandH();
}

bool Mysol::CoordianteDecent()
{
	double Wnorm_old = 0, Hnorm_old = 0;
	int itx = 0;
	for (itx = 0; itx < itmax; itx++)
	{
		// update W_{i~}
		for (int wi = 0; wi < lunum; wi++)
			//UpdateWi(wi);
			UpdateWi_v2(wi);

		// update H_{~i}
		for (int hi = 0; hi < levnum; hi++)
			//UpdateHi(hi);
			UpdateHi_v2(hi);

		// judge convergence
		double Wnorm = FrobeniusNorm(Matrix_W, lunum, dis);
		double Hnorm = FrobeniusNorm(Matrix_H, dis, levnum);

		if ((abs(Wnorm_old - Wnorm) / Wnorm < accu) &&
		        (abs(Hnorm_old - Hnorm) / Hnorm < accu))
			break;
		Wnorm_old = Wnorm;
		Hnorm_old = Hnorm;
	}
	if (itx == itmax)
		return false;
	else
		return true;
}

void Mysol::UpdateWi(int wi)
{
	double* w_temp = new double[dis];
	for (int i0 = 0; i0 < dis; i0++)
		w_temp[i0] = lambda * Matrix_W[wi][i0];

	for (int hi = 0; hi < levnum; hi++)
	{
		double coef = 0.0;
		int laps_be = -1;
		for (int i0 = 0; i0 < lunum; i0++)
		{
			double laps = Soc::LapsGet(wi, i0, &laps_be);
			coef += laps * DotW_H(i0, hi, dis);
		}

		if (uematrix[wi][hi])
			coef += alphax * (DotW_H(wi, hi, dis) - 1);
		else
		{
			double c = 0.0;
			if (d_flag)
				c = Soc::UEDistance(wi, hi);
			else
				c = 1;
			coef += (1 - alphax) * c * DotW_H(wi, hi, dis);
		}
		for (int ix = 0; ix < dis; ix++)
			w_temp[ix] += 2 * coef * Matrix_H[ix][hi];
	}

	for (int ix = 0; ix < dis; ix++)
		Matrix_W[wi][ix] -= eta * w_temp[ix];
}

void Mysol::UpdateHi(int hi)
{
	double* h_temp = new double[dis];
	for (int i0 = 0; i0 < dis; i0++)
		h_temp[i0] = lambda * Matrix_H[i0][hi];


	for (int wi = 0; wi < lunum; wi++)
	{
		double coef = 0.0;
		int laps_be = -1;
		for (int i0 = 0; i0 < lunum; i0++)
		{
			double laps = Soc::LapsGet(wi, i0, &laps_be);
			coef += laps * DotW_H(i0, hi, dis);
		}

		if (uematrix[wi][hi])
			coef += alphax * (DotW_H(wi, hi, dis) - 1);
		else
		{
			double c;
			if (d_flag)
				c = Soc::UEDistance(wi, hi);
			else
				c = 1;
			coef += (1 - alphax) * c * DotW_H(wi, hi, dis);
		}
		for (int ix = 0; ix < dis; ix++)
			h_temp[ix] += 2 * coef * Matrix_W[wi][ix];
	}

	for (int ix = 0; ix < dis; ix++)
		Matrix_H[ix][hi] -= eta * h_temp[ix];
}


void Mysol::ArryConst(double* arryin, double para, int arry_size)
{
	for (int i0 = 0; i0 < arry_size; i0++)
		arryin[i0] *= para;
}

double Mysol::DotW_H(int wi, int hi, int a_size)
{
	double res = 0;
	for (int i0 = 0; i0 < a_size; i0++)
		res += (Matrix_W[wi][i0] * Matrix_H[i0][hi]);
	return res;
}

double Mysol::FrobeniusNorm(double** matrix, int row, int column)
{
	double res = 0.0;
	for (int i0 = 0; i0 < row; i0++)
	{
		for (int i1 = 0; i1 < column; i1++)
			res += matrix[i0][i1] * matrix[i0][i1];
	}

	res = sqrt(res);
	return res;
}


void Mysol::UpdateWi_v2(int wi)
{
	double* w_temp = new double[dis];

	for (int i0 = 0; i0 < dis; i0++)
		w_temp[i0] = lambda * Matrix_W[wi][i0];

	double* laps = Soc::LapsGetRow(wi, lunum);

	for (int hi = 0; hi < levnum; hi++)
	{
		double coef = 0.0;
		if (uematrix[wi][hi])
			coef += alphax * (DotW_H(wi, hi, dis) - 1);
		else
		{
			double c = 0.0;
			if (d_flag)
				c = Soc::UEDistance(wi, hi);
			else
				c = 1;
			coef += (1 - alphax) * c * DotW_H(wi, hi, dis);
		}
		for (int ix = 0; ix < dis; ix++)
			w_temp[ix] += 2 * coef * Matrix_H[ix][hi] + 2 * LapsAndW(laps, ix, lunum);
	}

	for (int ix = 0; ix < dis; ix++)
		Matrix_W[wi][ix] -= eta * w_temp[ix];
	delete[] laps;
}

void Mysol::UpdateHi_v2(int hi)
{
	double* h_temp = new double[dis];
	for (int i0 = 0; i0 < dis; i0++)
		h_temp[i0] = lambda * Matrix_H[i0][hi];

	for (int wi = 0; wi < lunum; wi++)
	{
		double coef = 0.0;

		if (uematrix[wi][hi])
			coef += alphax * (DotW_H(wi, hi, dis) - 1);
		else
		{
			double c;
			if (d_flag)
				c = Soc::UEDistance(wi, hi);
			else
				c = 1;
			coef += (1 - alphax) * c * DotW_H(wi, hi, dis);
		}
		for (int ix = 0; ix < dis; ix++)
			h_temp[ix] += 2 * coef * Matrix_W[wi][ix];
	}

	for (int ix = 0; ix < dis; ix++)
		Matrix_H[ix][hi] -= eta * h_temp[ix];
}

double Mysol::LapsAndW(double *laps_row, int w_j, int a_size)
{
	double res = 0.0;
	for (int i0 = 0; i0 < a_size; i0++)
		res += laps_row[i0] * Matrix_W[i0][w_j];
	return res;
}

//---------------------------------------------------

double Mysol::EvalutionAUC()
{
	int test_num = Soc::GetUETestnum();

	double* positive = new double[test_num];

	// test set in Bool Matrix
	for (int i0 = 0; i0 < test_num; i0++)
	{
		uxpair temp = Soc::GetTestPair(i0);
		uematrix[temp.u][temp.x] = true;
		positive[i0] = DotW_H(temp.u, temp.x, dis);
	}

	double AUC_count = 0.0;
	for (int i0 = 0; i0 < lunum; i0++)
	{
		for (int i1 = 0; i1 < levnum; i1++)
		{
			if (!uematrix[i0][i1])
			{
				double scorex = DotW_H(i0, i1, dis);
				for (int i2 = 0; i2 < test_num; i2++)
				{
					if (positive[i2] > scorex)
						AUC_count += 1;
					else if (positive[i2] == scorex)
						AUC_count += 0.5;
				}
			}
		}
	}

	int neg_num = lunum * levnum - test_num;
	double AUC = AUC_count / (double(neg_num) * double(test_num));


	// test set remove
	for (int i0 = 0; i0 < test_num; i0++)
	{
		uxpair temp = Soc::GetTestPair(i0);
		uematrix[temp.u][temp.x] = false;
	}

	delete[] positive;
	return AUC;
}


bool Mysol::CCD_PlusPlus()
{
	double Wnorm_old = 0, Hnorm_old = 0;
	double Wnorm_max = -1, Hnorm_max = -1;
	int it = 0;
	for (it = 0; it < itmax; it++)
	{
		for (int k = 0; k < dis; k++)
		{
			double* u, *v;
			u = new double[lunum];
			for (int i0 = 0; i0 < lunum; i0++)
				u[i0] = 1;

			v = new double[levnum];
			for (int i0 = 0; i0 < levnum; i0++)
				v[i0] = 1;

			for (int wi = 0; wi < lunum; wi++)
			{
				double fenzi = 0;
				double fenmu = 0;
				for (int j = 0; j < levnum; j++)
				{
					if (uematrix[wi][j])
					{
						double x_wj = DotW_H(wi, j, dis) - Matrix_W[wi][k] * Matrix_H[k][j] - 1;
						fenzi += x_wj * v[j] * alphax;
						fenmu += v[j] * v[j] * alphax;
					}
					else
					{
						double y_wj = DotW_H(wi, j, dis) - Matrix_W[wi][k] * Matrix_H[k][j];
						double c_wj = 1;
						if (d_flag)
							c_wj = Soc::UEDistance(wi, j);
						fenzi += c_wj * y_wj * v[j] * (1 - alphax);
						fenmu += c_wj * v[j] * v[j] * (1 - alphax);
					}
				}
				double *lapxx = Soc::LapsGetRow(wi, lunum);
				for (int i0 = 0; i0 < lunum; i0++)
				{
					if (i0 == wi)
						continue;
					if (lapxx[i0] != 0)
						fenzi += lapxx[i0] * u[i0];
				}
				fenmu += lapxx[wi];

				delete[] lapxx;

				fenmu += lambda / 2;

				u[wi] = -fenzi / fenmu;
			}
			for (int i0 = 0; i0 < lunum; i0++)
				Matrix_W[i0][k] = u[i0];

			for (int hi = 0; hi < levnum; hi++)
			{
				double fenzi = 0;
				double fenmu = 0;
				for (int i = 0; i < lunum; i++)
				{
					if (uematrix[i][hi])
					{
						double x_wj = DotW_H(i, hi, dis) - Matrix_W[i][k] * Matrix_H[k][hi] - 1;
						fenzi += x_wj * u[i] * alphax;
						fenmu += u[i] * u[i] * alphax;
					}
					else
					{
						double y_wj = DotW_H(i, hi, dis) - Matrix_W[i][k] * Matrix_H[k][hi];
						double c_wj = 1;
						if (d_flag)
							c_wj = Soc::UEDistance(i, hi);
						fenzi += c_wj * y_wj * u[i] * (1 - alphax);
						fenmu += c_wj * u[i] * u[i] * (1 - alphax);
					}
				}
				fenmu += lambda / 2;
				v[hi] = -fenzi / fenmu;
			}
			for (int i0 = 0; i0 < lunum; i0++)
				Matrix_H[k][i0] = v[i0];

			delete[] u;
			delete[] v;
		}

		double Wnorm = FrobeniusNorm(Matrix_W, lunum, dis);
		double Hnorm = FrobeniusNorm(Matrix_H, dis, levnum);
		if (Wnorm > Wnorm_max)
			Wnorm_max = Wnorm;
		if (Hnorm > Hnorm_max)
			Hnorm_max = Hnorm;

		if ((abs(Wnorm - Wnorm_old) / Wnorm < accu) &&
		        (abs(Hnorm - Hnorm_old) / Wnorm < accu))
			break;
		Wnorm_old = Wnorm;
		Hnorm_old = Hnorm;
	}
	if (it == itmax)
		return false;
	else
		return true;
}

double Mysol::EvalutionAP(int u_node, double* pat1, double* pat3, double* pat5)
{
	double* rlist = new double[levnum];

	for (int i0 = 0; i0 < levnum; i0++)
	{
		if (uematrix[u_node][i0])
			rlist[i0] = 0;
		else
			rlist[i0] = DotW_H(u_node, i0, dis);
	}


	int* ulindex = new int[levnum];
	// memset(ulindex, 0, sizeof(int)* levnum);

	// Sort
	for (int i0 = 0; i0 < levnum; i0++)
	{
		double l_max = -1;
		int ixx = 0;
		for (int i1 = 0; i1 < levnum; i1++)
		{
			if (rlist[i1] > l_max)
			{
				l_max = rlist[i1];
				ixx = i1;
			}
		}
		ulindex[i0] = ixx;
		rlist[ixx] = -1;
	}
	delete[] rlist;

	int* reindex = new int[levnum];
	memset(reindex, 0, sizeof(int) * levnum);
	for (int i0 = 0; i0 < levnum; i0++)
		reindex[ulindex[i0]] = i0;

	delete[] ulindex;

	int nstart = Soc::GetTestIndex(u_node);
	int nend = Soc::GetTestIndex(u_node + 1);

	if (nend - nstart == 0)
	{
		*pat1 = 1;
		*pat3 = 1;
		*pat5 = 1;
		delete[] reindex;
		return 1;
	}

	int* prank = new int[nend - nstart];

	for (int i0 = nstart; i0 < nend; i0++)
	{
		uxpair temp = Soc::GetTestPair(i0);
		prank[i0 - nstart] = reindex[temp.x];
	}

	delete[] reindex;

	sort(prank, prank + nend - nstart);
	*pat1 = 1;
	*pat3 = 1;
	*pat5 = 1;
	if (prank[0] != 0)
		*pat1 = 0;

	double AP = 0.0;
	bool flag3 = true, flag5 = true;

	for (int i0 = 0; i0 < nend - nstart; i0++)
	{
		if (flag3 && (prank[i0] > 4))
		{
			(*pat3) = (i0 + 1) / 3;
			flag3 = false;
		}
		if (flag5 && (prank[i0] > 4))
		{
			(*pat5) = (i0 + 1) / 5;
			flag5 = false;
		}
		AP += double(i0 + 1) / double(prank[i0] + 1) / double(nend - nstart); 
	}


	delete[] prank;

	return AP;
}

double Mysol::EvalutionMAP(double* pat1, double* pat3, double* pat5)
{
	*pat1 = 0;
	*pat3 = 0;
	*pat5 = 0;

	double MAP = 0.0;
	double patx1 = 0, patx3 = 0, patx5 = 0, apx = 0;

	for (int i0 = 0; i0 < lunum; i0++)
	{
		apx = EvalutionAP(i0, &patx1, &patx3, &patx5);
		(*pat1) += patx1;
		(*pat3) += patx3;
		(*pat5) += patx5;
		MAP += apx;
	}

	*pat1 /= lunum;
	*pat3 /= lunum;
	*pat5 /= lunum;
	return MAP / lunum;
}