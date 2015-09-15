#ifndef LAPLACE_H
#define LAPLACE_H

#include <vector>

using std::vector;
using std::string;

struct uxpair
{
	int u, x;
};

struct Lap
{
	int row, column;
	double w;
};

struct Mydist
{
	double lon, lat;
};

class Soc
{
private:
	static int lunum;
	static int levnum;
	//static int lgnum;
	static int lue_num;
	//static int lug_num;
	static int test_num;
	static int train_num;

	// social
	static vector<uxpair> pairs;
	static vector<int> pindex; // train set will save here
	static vector<Lap> lap_on;
	static vector<Lap> lap_off;
	static vector<int> on_index;
	static vector<int> off_index;

	static vector<Lap> laps;
	static vector<int> s_index;

	static vector<uxpair> test;
	static vector<int> test_index;

	// distance
	static vector<Mydist> user_ll;
	static vector<Mydist> event_ll;

public:
	static bool SortPairsElement(uxpair a, uxpair b);
	static vector<Lap> Laplace(int unum, int xnum, string filename, int pairnum, bool off_flag);
	static bool SortLap(Lap a, Lap b);
	static void SetUpSocialLap(string name_on, string name_off, int unum, int evnum, int gnum, int ue_num, int ug_num);
	static void ReadInLocation(string name_user, string name_event);
	static double UEDistance(int usernode, int eventnode);
	static double RealDistance(int usernode, int eventnode);
	static void LapAdd(double rho_on, double rho_off);

	static int GetUnum();
	static int GetEvnum();
	static int GetUETrainnum();
	static int GetUETestnum();
	static uxpair GetPair(int present);
	static int GetPairIndex(int present);
	static uxpair GetTestPair(int present);
	static int GetTestIndex(int present);
	static double LapsGet(int u_i, int u_j, int* laps_be);
	static double* LapsGetRow(int wi, int lap_size);

	// sample test uepair
	static vector<uxpair> SampleTest();

	static bool SortLap2(Lap a, Lap b);

};

#endif