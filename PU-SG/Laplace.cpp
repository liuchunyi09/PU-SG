#include "Laplace.h"
#include "Limits.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;

int Soc::lunum = 0;
int Soc::levnum = 0;
//int Soc::lgnum = 0;
int Soc::lue_num = 0;
//int Soc::lug_num = 0;

int Soc::test_num = 0;
int Soc::train_num = 0;

vector<uxpair> Soc::pairs;
vector<int> Soc::pindex;

vector<Lap> Soc::lap_on;
vector<Lap> Soc::lap_off;
vector<int> Soc::on_index;
vector<int> Soc::off_index;

vector<Lap> Soc::laps;
vector<int> Soc::s_index;

vector<uxpair> Soc::test;
vector<int> Soc::test_index;

vector<Mydist> Soc::user_ll;
vector<Mydist> Soc::event_ll;

//======================================================================

bool Soc::SortPairsElement(uxpair a, uxpair b)
{
	if (a.u < b.u)
		return true;
	else if ((a.u == b.u) && (a.x < b.x))
		return true;
	else
		return false;
}

vector<Lap> Soc::Laplace(int unum, int xnum, string filename, int pairnum, bool off_flag)
{
	ifstream finc(filename);
	pairs.clear();
	pindex.clear();

	pairs.resize(pairnum);
	for (int i0 = 0; i0 < pairnum; i0++)
	{
		int u, x;
		finc >> u >> x;
		pairs[i0].u = u - 1;
		pairs[i0].x = x - 1;
	}

	finc.close();

	sort(pairs.begin(), pairs.end(), SortPairsElement); // sort the pairs

	if (off_flag)
	{
		pairs = Soc::SampleTest();
		pairnum = train_num;
	}



	// set up index
	pindex.resize(unum, 0);
	for (int i0 = 0; i0 < pairnum; i0++)
		pindex[pairs[i0].u] = i0 + 1;
	for (int i0 = 1; i0 < unum; i0++)
	{
		if (pindex[i0] < pindex[i0 - 1])
			pindex[i0] = pindex[i0 - 1];
	}

	vector<Lap> laplace;
	for (int i0 = 0; i0 < unum; i0++)
	{
		int number = 0;
		if (i0 == 0)
			number = pindex[i0];
		else 
			number = pindex[i0] - pindex[i0 - 1];

		if (number == 0)
			continue;

		for (int i1 = i0 + 1; i1 < unum; i1++)
		{
			if (pindex[i1] - pindex[i1 - 1] == 0)
				continue;

			int j0 = 0;
			if (i0 != 0)
				j0 = pindex[i0 - 1];
			int j1 = pindex[i1 - 1];
			int in = 0;
			int un = number + pindex[i1] - pindex[i1 - 1];
			Lap temp;
			while ((j0 < pindex[i0]) && (j1 < pindex[i1]))
			{
				if (pairs[j0].x < pairs[j1].x)
					j0++;
				else if (pairs[j0].x > pairs[j1].x)
					j1++;
				else
				{
					in++;
					j0++;
					j1++;
				}
			}
			if (in != 0)
			{
				temp.row = i0;
				temp.column = i1;
				temp.w = - double(in) / double(un - in);
				laplace.push_back(temp);
				temp.column = i0;
				temp.row = i1;
				laplace.push_back(temp);
			}
		}
	}

	// the diagonal
	vector<double> diag(unum, 0);
	for (int i0 = 0; i0 < laplace.size(); i0++)
		diag[laplace[i0].row] -= laplace[i0].w;

	for (int i0 = 0; i0 < unum; i0++)
	{
		Lap temp;
		temp.row = i0;
		temp.column = i0;
		temp.w = diag[i0];
		if (diag[i0] != 0)
			laplace.push_back(temp);
	}

	sort(laplace.begin(), laplace.end(), Soc::SortLap);
	return laplace;
}

bool Soc::SortLap(Lap a, Lap b)
{
	if (a.row < b.row)
		return true;
	else if ((a.row == b.row) && (a.column < b.column))
		return true;
	else
		return false;
}

void Soc::SetUpSocialLap(string name_on, string name_off, int unum, int evnum, int gnum, int ue_num, int ug_num)
{	
	lunum = unum;
	levnum = evnum;
	lue_num = ue_num;

	lap_on = Laplace(unum, gnum, name_on, ug_num, false);
	lap_off = Laplace(unum, evnum, name_off, ue_num, true);

	// pindex.clear();
	// pairs.clear();

	// set up index
	on_index.resize(unum, 0);
	for (int i0 = 0; i0 < lap_on.size(); i0++)
		on_index[lap_on[i0].row] = i0 + 1;

	for (int i0 = 1; i0 < unum; i0++)
	{
		if (on_index[i0] < on_index[i0 - 1])
			on_index[i0] = on_index[i0 - 1];
	}

	off_index.resize(unum, 0);
	for (int i0 = 0; i0 < lap_off.size(); i0++)
		off_index[lap_off[i0].row] = i0 + 1;

	for (int i0 = 1; i0 < unum; i0++)
	{
		if (off_index[i0] < off_index[i0 - 1])
			off_index[i0] = off_index[i0 - 1];
	}
}

void Soc::ReadInLocation(string name_user, string name_event)
{
	ifstream uin(name_user);
	ifstream ein(name_event);

	user_ll.resize(lunum);
	event_ll.resize(levnum);
	for (int i0 = 0; i0 < lunum; i0++)
		uin >> user_ll[i0].lon >> user_ll[i0].lat;

	for (int i0 = 0; i0 < levnum; i0++)
		ein >> event_ll[i0].lon >> event_ll[i0].lat;

	uin.close();
	ein.close();
}

double Soc::UEDistance(int usernode, int eventnode)
{
	double lon0 = user_ll[usernode].lon;
	double lat0 = user_ll[usernode].lat;
	double lon1 = event_ll[eventnode].lon;
	double lat1 = event_ll[eventnode].lat;

	double c = sin(lat0 * PI / 180) * sin(lat1 * PI / 180) * cos((lon0 - lon1) * PI / 180) + cos(lat0 * PI / 180) * cos(lat1 * PI / 180);
	if (c > 1)
		c = 1;
	if (c < -1)
		c = -1;

	double dis = acos(c) * EARTH;

	if (dis < 10)
		return 0;
	else if (dis > 100)
		return 1;
	else 
		return exp((dis - 100) / 50);
}

double Soc::RealDistance(int usernode, int eventnode)
{
	double lon0 = user_ll[usernode].lon;
	double lat0 = user_ll[usernode].lat;
	double lon1 = event_ll[eventnode].lon;
	double lat1 = event_ll[eventnode].lat;

	double c = sin(lat0 * PI / 180) * sin(lat1 * PI / 180) * cos((lon0 - lon1) * PI / 180) + cos(lat0 * PI / 180) * cos(lat1 * PI / 180);
	if (c > 1)
		c = 1;
	else if (c < -1)
		c = -1;
	double dis = acos(c) * EARTH;
	return dis;
}

void Soc::LapAdd(double rho_on, double rho_off)
{
	laps.clear();
	if ((rho_on == 0.0) && (rho_off != 0.0))
		laps = lap_off;
	else if ((rho_off == 0.0) && (rho_on != 0.0))
		laps = lap_on;
	else if ((rho_off == 0.0) && (rho_on == 0.0))
		laps.resize(0);
	else
	{
		laps.reserve(lap_on.size() + lap_off.size());

		auto it0 = lap_on.begin();
		auto it1 = lap_off.begin();

		Lap temp;
		while ((it0 != lap_on.end()) && (it1 != lap_off.end()))
		{
			if ((*it0).row < (*it1).row)
			{
				temp = *it0;
				temp.w = (*it0).w * rho_on;
				laps.push_back(temp);
				it0++;
			}
			else if ((*it0).row >(*it1).row)
			{
				temp = *it1;
				temp.w = (*it1).w * rho_off;
				laps.push_back(temp);
				it1++;
			}
			else
			{
				if ((*it0).column < (*it1).column)
				{
					temp = *it0;
					temp.w = (*it0).w * rho_on;
					laps.push_back(temp);
					it0++;
				}
				else if ((*it0).column >(*it1).column)
				{
					temp = *it1;
					temp.w = (*it1).w * rho_off;
					laps.push_back(temp);
					it1++;
				}
				else
				{
					temp = *it0;
					temp.w = (*it0).w * rho_on + (*it1).w * rho_off;
					laps.push_back(temp);
					it0++;
					it1++;
				}
			}
		}

		while (it0 != lap_on.end())
		{
			temp = *it0;
			temp.w = (*it0).w * rho_on;
			laps.push_back(temp);
			it0++;
		}

		while (it1 != lap_off.end())
		{
			temp = *it1;
			temp.w = (*it1).w * rho_off;
			laps.push_back(temp);
			it1++;
		}
	}

	// sort(laps.begin(), laps.end(), SortLap2);
	// index
	s_index.resize(lunum);
	for (int i0 = 0; i0 < laps.size(); i0++)
		s_index[laps[i0].row] = i0 + 1;

	for (int i0 = 1; i0 < lunum; i0++)
	{
		if (s_index[i0] < s_index[i0 - 1])
			s_index[i0] = s_index[i0 - 1];
	}
}

int Soc::GetUnum()
{
	return lunum;
}

int Soc::GetEvnum()
{
	return levnum;
}

int Soc::GetUETrainnum()
{
	return train_num;
}

int Soc::GetUETestnum()
{
	return test_num;
}

uxpair Soc::GetPair(int present)
{
	return pairs[present];
}

uxpair Soc::GetTestPair(int present)
{
	return test[present];
}

double Soc::LapsGet(int u_i, int u_j, int* laps_be)
{
	int nstart = 0;
	if (u_i != 0)
		nstart = s_index[u_i - 1];
	if (*laps_be != -1)
		nstart = *laps_be;
	for (int i0 = nstart; i0 < s_index[u_i]; i0++)
	{
		if (laps[i0].column == u_j)
		{
			*laps_be = i0;
			return laps[i0].w;
		}
		else if (laps[i0].column > u_j)
		{
			*laps_be = i0;
			return 0;
		}
	}
	return 0;
}

vector<uxpair> Soc::SampleTest()
{
	test_num = int(RATIO * lue_num);
	train_num = lue_num - test_num;

	test.resize(0);
	vector<uxpair> train(0);

	// srand(unsigned(time(nullptr)));

	int tr = 0, te = 0;
	for (int i0 = 0; i0 < lue_num; i0++)
	{
		if (rand() % lue_num < test_num)
		{
			test.push_back(pairs[i0]);
			te++;
		}
		else
		{
			tr++;
			train.push_back(pairs[i0]);
		}
	}

	test_num = te;
	train_num = tr;

	// set up test index
	test_index.resize(lunum, 0);
	for (int i0 = 0; i0 < te; i0++)
		test_index[test[i0].u] = i0 + 1;
	for (int i0 = 1; i0 < lunum; i0++)
	{
		if (test_index[i0] < test_index[i0 - 1])
			test_index[i0] = test_index[i0 - 1];
	}

	return train;
}

double* Soc::LapsGetRow(int wi, int lap_size)
{
	double* laps_row = new double[lap_size];
	memset(laps_row, 0, sizeof(double)* lap_size);
	int nstart = 0;
	if (wi != 0)
		nstart = s_index[wi - 1];
	for (int i0 = nstart; i0 < s_index[wi]; i0++)
		laps_row[laps[i0].column] = laps[i0].w;
	return laps_row;
}

bool Soc::SortLap2(Lap a, Lap b)
{
	if (a.row < b.row)
		return true;
	else if ((a.row == b.row) && (a.w > b.w))
		return true;
	else
		return false;
}

int Soc::GetTestIndex(int present)
{
	if (present == 0)
		return 0;
	else
		return test_index[present - 1];
}

int Soc::GetPairIndex(int present)
{
	if (present = 0)
		return 0;
	else
		return pindex[present - 1];
}