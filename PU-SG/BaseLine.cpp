#include "BaseLine.h"
#include "Laplace.h"
#include "Limits.h"
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

using namespace std;

int Mybase::lunum = 0;
int Mybase::levnum = 0;
int Mybase::dis = 0;
double** Mybase::Matrix_W = nullptr;
double** Mybase::Matrix_H = nullptr;

vector<uxpair> Mybase::vu_p;
vector<int> Mybase::vindexp;

vector<uxpair> Mybase::vu_n;
vector<int> Mybase::vindexn;


double Mybase::EvalutionRandom(double* pat1, double* pat3, double* pat5)
{
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();
	
	*pat1 = 0;
	*pat3 = 0;
	*pat5 = 0;
	double MAP = 0;
	
	for (int u = 0; u < lunum; u++)
	{
		int *rlist = new int[levnum];
		int countx = 0;
		vector<int> numpool(levnum, 0);
		for (int i0 = 0; i0 < levnum; i0++)
			numpool[i0] = i0;

		while (countx < levnum)
		{
			int xx = rand() % (levnum - countx);
			rlist[countx++] = numpool[xx];
			numpool.erase(numpool.begin() + xx);
		}

		int* reindex = new int[levnum];
		for (int i0 = 0; i0 < levnum; i0++)
			reindex[rlist[i0]] = i0;

		int trstart = Soc::GetPairIndex(u);
		int trend = Soc::GetPairIndex(u + 1);
		
		int* trnode = new int[trend - trstart];
		for (int i0 = 0; i0 < trend - trstart; i0++)
			trnode[i0] = reindex[Soc::GetPair(i0 + trstart).x]; // node rank

		sort(trnode, trnode + trend - trstart);

		for (int i0 = 0; i0 < trend - trstart; i0++)
		{
			for (int i1 = trnode[i0] - i0; i1 < levnum - 1 - i0; i1++)
				rlist[i1] = rlist[i1 + 1];
		}

		for (int i0 = 0; i0 < levnum - (trend - trstart); i0++)
			reindex[rlist[i0]] = i0;

		delete[] trnode;
		delete[] rlist;

		int nstart = Soc::GetTestIndex(u);
		int nend = Soc::GetTestIndex(u + 1);
		if (nend - nstart == 0)
		{
			continue; 
		}

		int* prank = new int[nend - nstart];

		for (int i0 = nstart; i0 < nend; i0++)
		{
			uxpair temp = Soc::GetTestPair(i0);
			prank[i0 - nstart] = reindex[temp.x];
		}

		delete[] reindex;
		
		sort(prank, prank + nend - nstart);

		double patx1 = 1, patx3 = 1, patx5 = 1, ap = 0;
		bool flag3 = true, flag5 = true;
		if (prank[0] != 0)
			patx1 = 0;

		for (int i0 = 0; i0 < nend - nstart; i0++)
		{
			if (flag3 && (prank[i0] > 2))
			{
				(*pat3) = (i0 + 1) / 3;
				flag3 = false;
			}
			if (flag5 && (prank[i0] > 4))
			{
				(*pat5) = (i0 + 1) / 5;
				flag5 = false;
			}
			ap += double(i0 + 1) / double(prank[i0] + 1) / double(nend - nstart);
		}

		delete[] prank;

		*pat1 += patx1;
		*pat3 += patx3;
		*pat5 += patx5;
		MAP += ap;
	}

	*pat1 /= lunum;
	*pat3 /= lunum;
	*pat5 /= lunum;
	MAP /= lunum;

	return MAP;
}

double Mybase::EvalutionMostPo(double* pat1, double* pat3, double* pat5)
{
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();

	*pat1 = 0;
	*pat3 = 0;
	*pat5 = 0;
	double MAP = 0;

	for (int u = 0; u < lunum; u++)
	{
		// generate list
		int *rlist = new int[levnum];
		
		vector<int> enumc(levnum, 0);
		int uenum = Soc::GetUETrainnum();
		for (int i0 = 0; i0 < uenum; i0++)
			enumc[Soc::GetPair(i0).x]++;

		for (int i0 = 0; i0 < levnum; i0++)
		{
			int imax = -1;
			int numx = 1;
			for (int i1 = 0; i1 < levnum; i1++)
			{
				if (enumc[i1] > imax)
				{
					imax = enumc[i1];
					numx = i1;
				}
			}
			rlist[i0] = numx;
			enumc[numx] = -1;
		}

		// remove exist user-event

		int* reindex = new int[levnum];
		for (int i0 = 0; i0 < levnum; i0++)
			reindex[rlist[i0]] = i0;

		int trstart = Soc::GetPairIndex(u);
		int trend = Soc::GetPairIndex(u + 1);

		int* trnode = new int[trend - trstart];
		for (int i0 = 0; i0 < trend - trstart; i0++)
			trnode[i0] = reindex[Soc::GetPair(i0 + trstart).x]; // node rank

		sort(trnode, trnode + trend - trstart);

		for (int i0 = 0; i0 < trend - trstart; i0++)
		{
			for (int i1 = trnode[i0] - i0; i1 < levnum - 1 - i0; i1++)
				rlist[i1] = rlist[i1 + 1];
		}

		for (int i0 = 0; i0 < levnum - (trend - trstart); i0++)
			reindex[rlist[i0]] = i0;

		delete[] trnode;
		delete[] rlist;

		int nstart = Soc::GetTestIndex(u);
		int nend = Soc::GetTestIndex(u + 1);
		if (nend - nstart == 0)
		{
			continue; 
		}

		int* prank = new int[nend - nstart];

		for (int i0 = nstart; i0 < nend; i0++)
		{
			uxpair temp = Soc::GetTestPair(i0);
			prank[i0 - nstart] = reindex[temp.x];
		}

		delete[] reindex;

		sort(prank, prank + nend - nstart);

		double patx1 = 1, patx3 = 1, patx5 = 1, ap = 0;
		bool flag3 = true, flag5 = true;
		if (prank[0] != 0)
			patx1 = 0;

		for (int i0 = 0; i0 < nend - nstart; i0++)
		{
			if (flag3 && (prank[i0] > 2))
			{
				(*pat3) = (i0 + 1) / 3;
				flag3 = false;
			}
			if (flag5 && (prank[i0] > 4))
			{
				(*pat5) = (i0 + 1) / 5;
				flag5 = false;
			}
			ap += double(i0 + 1) / double(prank[i0] + 1)  / double(nend - nstart);
		}
		
		delete[] prank;

		*pat1 += patx1;
		*pat3 += patx3;
		*pat5 += patx5;
		MAP += ap;
	}

	*pat1 /= lunum;
	*pat3 /= lunum;
	*pat5 /= lunum;
	MAP /= lunum;

	return MAP;
}

double Mybase::EvalutionLocAwa(double* pat1, double* pat3, double* pat5)
{
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();

	*pat1 = 0;
	*pat3 = 0;
	*pat5 = 0;
	double MAP = 0;

	for (int u = 0; u < lunum; u++)
	{
		// generate list
		int *rlist = new int[levnum];

		vector<double> ue_dis(levnum, 0);
		
		for (int i0 = 0; i0 < levnum; i0++)
			ue_dis[i0] = Soc::RealDistance(u, i0);

		for (int i0 = 0; i0 < levnum; i0++)
		{
			double dis_min = 10000000.0;
			int node = 0;
			for (int i1 = 0; i1 < levnum ; i1++)
			{
				if (dis_min > ue_dis[i1])
				{
					dis_min = ue_dis[i1];
					node = i1;
				}
			}
			rlist[i0] = node;
			ue_dis[node] = 10000000;
		}

		// remove exist user-event

		int* reindex = new int[levnum];
		for (int i0 = 0; i0 < levnum; i0++)
			reindex[rlist[i0]] = i0;

		int trstart = Soc::GetPairIndex(u);
		int trend = Soc::GetPairIndex(u + 1);

		int* trnode = new int[trend - trstart];
		for (int i0 = 0; i0 < trend - trstart; i0++)
			trnode[i0] = reindex[Soc::GetPair(i0 + trstart).x]; // node rank

		sort(trnode, trnode + trend - trstart);

		for (int i0 = 0; i0 < trend - trstart; i0++)
		{
			for (int i1 = trnode[i0] - i0; i1 < levnum - 1 - i0; i1++)
				rlist[i1] = rlist[i1 + 1];
		}

		for (int i0 = 0; i0 < levnum - (trend - trstart); i0++)
			reindex[rlist[i0]] = i0;

		delete[] trnode;
		delete[] rlist;

		int nstart = Soc::GetTestIndex(u);
		int nend = Soc::GetTestIndex(u + 1);
		if (nend - nstart == 0)
		{
			continue; 
		}

		int* prank = new int[nend - nstart];

		for (int i0 = nstart; i0 < nend; i0++)
		{
			uxpair temp = Soc::GetTestPair(i0);
			prank[i0 - nstart] = reindex[temp.x];
		}

		delete[] reindex;

		sort(prank, prank + nend - nstart);

		double patx1 = 1, patx3 = 1, patx5 = 1, ap = 0;
		bool flag3 = true, flag5 = true;
		if (prank[0] != 0)
			patx1 = 0;

		for (int i0 = 0; i0 < nend - nstart; i0++)
		{
			if (flag3 && (prank[i0] > 2))
			{
				(*pat3) = (i0 + 1) / 3;
				flag3 = false;
			}
			if (flag5 && (prank[i0] > 4))
			{
				(*pat5) = (i0 + 1) / 5;
				flag5 = false;
			}
			ap += double(i0 + 1) / double(prank[i0] + 1) / double(nend - nstart);
		}

		delete[] prank;

		*pat1 += patx1;
		*pat3 += patx3;
		*pat5 += patx5;
		MAP += ap;
	}

	*pat1 /= lunum;
	*pat3 /= lunum;
	*pat5 /= lunum;
	MAP /= lunum;

	return MAP;
}

double Mybase::EvalutionUserKNN(double* pat1, double* pat3, double* pat5)
{
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();

	*pat1 = 0;
	*pat3 = 0;
	*pat5 = 0;
	double MAP = 0;

	for (int u = 0; u < lunum; u++)
	{
		// generate list
		int *rlist = new int[levnum];

		double* simlar = Soc::LapsGetRow(u, lunum);
		
		vector<double> u_rating(levnum, 0);

		int ttrain_num = Soc::GetUETrainnum();
		for (int i0 = 0; i0 < ttrain_num; i0++)
		{
			uxpair temp = Soc::GetPair(i0);
			if (temp.u != u)
				u_rating[temp.x] -= simlar[temp.u];
		}

		delete[] simlar;

		for (int i0 = 0; i0 < levnum; i0++)
		{
			double rating_max = -1.0;
			int node = 0;
			for (int i1 = 0; i1 < levnum - i0; i1++)
			{
				if (rating_max < u_rating[i1])
				{
					rating_max = u_rating[i1];
					node = i1;
				}
			}
			rlist[i0] = node;
			u_rating[node] = -1;
		}

		// remove exist user-event

		int* reindex = new int[levnum];
		for (int i0 = 0; i0 < levnum; i0++)
			reindex[rlist[i0]] = i0;

		int trstart = Soc::GetPairIndex(u);
		int trend = Soc::GetPairIndex(u + 1);

		int* trnode = new int[trend - trstart];
		for (int i0 = 0; i0 < trend - trstart; i0++)
			trnode[i0] = reindex[Soc::GetPair(i0 + trstart).x]; // node rank

		sort(trnode, trnode + trend - trstart);

		for (int i0 = 0; i0 < trend - trstart; i0++)
		{
			for (int i1 = trnode[i0] - i0; i1 < levnum - 1 - i0; i1++)
				rlist[i1] = rlist[i1 + 1];
		}

		for (int i0 = 0; i0 < levnum - (trend - trstart); i0++)
			reindex[rlist[i0]] = i0;

		delete[] trnode;
		delete[] rlist;

		int nstart = Soc::GetTestIndex(u);
		int nend = Soc::GetTestIndex(u + 1);
		if (nend - nstart == 0)
		{
			continue; 
		}
		
		int* prank = new int[nend - nstart];

		for (int i0 = nstart; i0 < nend; i0++)
		{
			uxpair temp = Soc::GetTestPair(i0);
			prank[i0 - nstart] = reindex[temp.x];
		}

		delete[] reindex;

		sort(prank, prank + nend - nstart);

		double patx1 = 1, patx3 = 1, patx5 = 1, ap = 0;
		bool flag3 = true, flag5 = true;
		if (prank[0] != 0)
			patx1 = 0;

		for (int i0 = 0; i0 < nend - nstart; i0++)
		{
			if (flag3 && (prank[i0] > 2))
			{
				(*pat3) = (i0 + 1) / 3;
				flag3 = false;
			}
			if (flag5 && (prank[i0] > 4))
			{
				(*pat5) = (i0 + 1) / 5;
				flag5 = false;
			}
			ap += double(i0 + 1) / double(prank[i0] + 1) / double(nend - nstart);
		}

		delete[] prank;

		*pat1 += patx1;
		*pat3 += patx3;
		*pat5 += patx5;
		MAP += ap;
	}

	*pat1 /= lunum;
	*pat3 /= lunum;
	*pat5 /= lunum;
	MAP /= lunum;

	return MAP;
}

void Mybase::MFInitial(int in_lunum, int in_levnum, int in_dis, int in_ratio)
{
	lunum = in_lunum;
	levnum = in_levnum;
	dis = in_dis;
	ratio = in_ratio;

	Matrix_W = new double*[lunum];
	for (int i0 = 0; i0 < lunum; i0++)
	{
		Matrix_W[i0] = new double[dis];
		memset(Matrix_W[i0], 0, sizeof(double)* dis);
	}

	Matrix_H = new double*[levnum];
	for (int i0 = 0; i0 < levnum; i0++)
	{
		Matrix_H[i0] = new double[dis];
		memset(Matrix_H[i0], 0, sizeof(double)* dis);
	}
}

void Mybase::MFFreeWandH()
{
	for (int i0 = 0; i0 < lunum; i0++)
		delete[] Matrix_W[i0];
	delete[] Matrix_W;

	for (int i0 = 0; i0 < levnum; i0++)
		delete[] Matrix_H[i0];
	delete[] Matrix_H;
}

void Mybase::MFSetUpTrainPair()
{
	// set positive
	int ipn = Soc::GetUETrainnum();
	vu_p.resize(ipn);
	for (int i0 = 0; i0 < ipn; i0++)
		vu_p[i0] = Soc::GetPair(i0);

	vindexp.resize(lunum, 0);
	for (int i0 = 0; i0 < ipn; i0++)
		vindexp[vu_p[i0].u] = i0 + 1;
	for (int i0 = 1; i0 < lunum; i0++)
	{
		if (vindexp[i0] < vindexp[i0 - 1])
			vindexp[i0] = vindexp[i0 - 1];
	}

	// set negative
	bool* beex = new bool[levnum];
	vu_n.resize(ratio * ipn);
	int inn = 0;
	for (int i0 = 0; i0 < lunum; i0++)
	{
		memset(beex, false, sizeof(bool)* levnum);
		int nstart = 0;
		if (i0 != 0)
			nstart = vindexp[i0 - 1];
		int nend = vindexp[i0];

		for (int i1 = nstart; i1 < nend; i1++)
			beex[vu_p[i1].x] = true;

		double dposs = double(ratio * (nend - nstart)) / double(levnum - nend + nstart);
		int ineg_u = ratio * (nend - nstart);
		int countx = 0;
		while (countx < ineg_u)
		{
			for (int i1 = 0; i1 < levnum; i1++)
			{
				if (beex[i1])
					continue;
				if (double(rand()) / double(RAND_MAX) < dposs)
				{
					vu_n[inn + countx].u = i0;
					vu_n[inn + countx].x = i1;
					beex[i1] = true;
					dposs = double(ratio * (nend - nstart) - countx) / double(levnum - countx - nend + nstart);
					countx++;
				}
				if (countx == ineg_u)
					break;
			}
		}
		inn += countx;
	}

	delete[] beex;

	vindexn.resize(lunum, 0);
	for (int i0 = 0; i0 < ipn * ratio; i0++)
		vindexn[vu_n[i0].u] = i0 + 1;
	for (int i0 = 1; i0 < lunum; i0++)
	{
		if (vindexn[i0] < vindexn[i0 - 1])
			vindexn[i0] = vindexn[i0 - 1];
	}
}

bool Mybase::MFTraining(int itermax, double drate, double dlambda)
{
	for (int it = 0; it < itermax; it++)
	{
		for (int iu = 0; iu < lunum; iu++)
		{
			int nstart = 0;
			if (iu != 0)
				nstart = vindexp[iu - 1];
			int nend = vindexp[iu];

			int ipsample = rand() % (nend - nstart);
			int insample = rand() % (ratio * (nend - nstart));
			int nnstart = 0;
			if (iu != 0)
				nnstart = vindexp[iu - 1];

			double dx_ui = MFDot(Matrix_W[iu], Matrix_H[ipsample + nstart], dis);
			double dx_uj = MFDot(Matrix_W[iu], Matrix_H[insample + nnstart], dis);

			if (it == 0)
			{
				dx_ui = 1;
				dx_uj = 0;
			}

			double dx_uij = dx_ui - dx_uj;
			double dpara = exp(-dx_uij) / (1 + exp(-dx_uij));

			for (int i0 = 0; i0 < dis; i0++)
			{
				Matrix_W[iu][i0] *= (dlambda * drate + 1);
				Matrix_W[iu][i0] += drate * dpara * (Matrix_H[ipsample + nstart][i0] - Matrix_H[insample + nnstart][i0]);
			}

			for (int i0 = 0; i0 < dis; i0++)
			{
				Matrix_H[ipsample + nstart][i0] *= (dlambda * drate + 1);
				Matrix_H[insample + nnstart][i0] *= (dlambda * drate + 1);
				Matrix_H[ipsample + nstart][i0] += drate * dpara * Matrix_W[iu][i0];
				Matrix_H[insample + nnstart][i0] -= drate * dpara * Matrix_W[iu][i0];
			}
		}
	}
}

double Mybase::MFDot(const double *a, const double *b, int num)
{
	double dres = 0;
	for (int i0 = 0; i0 < num; i0++)
		dres += a[i0] * b[i0];

	return dres;
}

double Mybase::MFEvalution(double* pat1, double* pat3, double* pat5)
{
	lunum = Soc::GetUnum();
	levnum = Soc::GetEvnum();

	*pat1 = 0;
	*pat3 = 0;
	*pat5 = 0;
	double MAP = 0;

	for (int u = 0; u < lunum; u++)
	{
		// generate list
		int *rlist = new int[levnum];

		vector<double> score(levnum, 0);
		for (int i0 = 0; i0 < levnum; i0++)
			dscore[i0] = MFDot(Matrix_W[u], Matrix_H[i0], dis);
		for (int i0 = 0; i0 < levnum; i0++)
		{
			double dmaxs = 0.0;
			int imaxn = 0;
			for (int i1 = 0; i1 < levnum - i0; i1++)
			{
				if (dmaxs < dscore[i1])
				{
					dmaxs = dscore[i1];
					imaxn = i1;
				}
			}
			rlist[i0] = imaxn;
			score[imaxn] = -1;
		}

		delete[] dscore;

		// remove exist user-event

		int* reindex = new int[levnum];
		for (int i0 = 0; i0 < levnum; i0++)
			reindex[rlist[i0]] = i0;

		int trstart = Soc::GetPairIndex(u);
		int trend = Soc::GetPairIndex(u + 1);

		int* trnode = new int[trend - trstart];
		for (int i0 = 0; i0 < trend - trstart; i0++)
			trnode[i0] = reindex[Soc::GetPair(i0 + trstart).x]; // node rank

		sort(trnode, trnode + trend - trstart);

		for (int i0 = 0; i0 < trend - trstart; i0++)
		{
			for (int i1 = trnode[i0] - i0; i1 < levnum - 1 - i0; i1++)
				rlist[i1] = rlist[i1 + 1];
		}

		for (int i0 = 0; i0 < levnum - (trend - trstart); i0++)
			reindex[rlist[i0]] = i0;

		delete[] trnode;
		delete[] rlist;

		int nstart = Soc::GetTestIndex(u);
		int nend = Soc::GetTestIndex(u + 1);
		if (nend - nstart == 0)
		{
			continue; 
		}
		
		int* prank = new int[nend - nstart];

		for (int i0 = nstart; i0 < nend; i0++)
		{
			uxpair temp = Soc::GetTestPair(i0);
			prank[i0 - nstart] = reindex[temp.x];
		}

		delete[] reindex;

		sort(prank, prank + nend - nstart);

		double patx1 = 1, patx3 = 1, patx5 = 1, ap = 0;
		bool flag3 = true, flag5 = true;
		if (prank[0] != 0)
			patx1 = 0;

		for (int i0 = 0; i0 < nend - nstart; i0++)
		{
			if (flag3 && (prank[i0] > 2))
			{
				(*pat3) = (i0 + 1) / 3;
				flag3 = false;
			}
			if (flag5 && (prank[i0] > 4))
			{
				(*pat5) = (i0 + 1) / 5;
				flag5 = false;
			}
			ap += double(i0 + 1) / double(prank[i0] + 1) / double(nend - nstart);
		}

		delete[] prank;

		*pat1 += patx1;
		*pat3 += patx3;
		*pat5 += patx5;
		MAP += ap;
	}

	*pat1 /= lunum;
	*pat3 /= lunum;
	*pat5 /= lunum;
	MAP /= lunum;

	return MAP;
}