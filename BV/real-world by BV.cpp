#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include <tchar.h>
#include "fstream"
#include <string.h>
#include <time.h>
#include <windows.h>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

#define Separate_Run 1000000  
#define TryTime 10
#define T 3		//Player¡¯s temptation to defect 
#define H -1		//Hurt suffered by others as a result of an agent¡¯s defection 
#define P -9	//Cost of being punished 
#define E -2		// Enforcement cost,as known as cost of applying punishment 
#define LS 0.142857		//Learning Step  
#define ER 0.01	//Explore Rate 


struct Agent
{
	vector<vector<int>> neighbours;//constains players's neighbours
	vector<int> neighbours1;
	double boldness; //decide the posibility to defect
	double vengefulness;//decide the posibility to punish somebody
	double o;//opportunity of defect;
	double S;//opportunity of seeing an agent defects;
	double DS;//the defection score,equals to total temptation rewards minus total punishment;
	double PS;//the punishment score,equals to the enforcement cost due to punishment or metaounishment;
	double HS;
	double POS;//the punishment omisstion score,equals to the value  due to being metapunish;
	int LN;
	int degree_g;
	int degree1;
	int Punish;
	int learn;
	double TS;
	int LR; 

};
struct Agent Player[402511];

int sum = 0;
int flag = 0;
std::vector<vector<int>> ZH;

//randf() generate 0-1,randi(x)generate 0-(x-1)
/*************************** RNG (Random Number Generator)procedures ****************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti = NN + 1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{
	int i;
	for (i = 0; i < NN; i++)
	{
		mt[i] = seed & 0xffff0000;
		seed = 69069 * seed + 1;
		mt[i] |= (seed & 0xffff0000) >> 16;
		seed = 69069 * seed + 1;
	}
	mti = NN;
}
void lsgenrand(unsigned long seed_array[])
{
	int i;
	for (i = 0; i < NN; i++)
		mt[i] = seed_array[i];
	mti = NN;
}
double genrand()
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0, MATRIX_A };
	if (mti >= NN)
	{
		int kk;
		if (mti == NN + 1) sgenrand(4357);
		for (kk = 0; kk < NN - MM; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MM] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (; kk < NN - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MM - NN)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[NN - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[NN - 1] = mt[MM - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
	return y;
}

double randf() { return ((double)genrand() * 2.3283064370807974e-10); }
long randi(unsigned long LIM) { return((unsigned long)genrand() % LIM); }
void Initial(void)
{


	for (int i = 0; i < 402511; i++)
	{

		Player[i].boldness = randf();
		Player[i].vengefulness = randf();
		Player[i].S = randf();
		Player[i].DS = 0;
		Player[i].PS = 0;
		Player[i].HS = 0;
		Player[i].POS = 0;
		Player[i].TS = 0;
		Player[i].LN = 0;
		Player[i].Punish = 0;
		Player[i].learn = 0;
		Player[i].LR = 0;

	}
}
/*************************** paregraph ****************************************/
void generate_neigh(vector<int> v) {

	vector<int> edge_sta, edge_end, neigh;
	long agent = 0;

	for (int l = 0; l < v.size(); l++) {

		edge_sta.push_back(v[l]);
		edge_end.push_back(v[l]);

	}
	for (int l = 0; l < v.size() - 1; l++) {

		for (int l1 = l + 1; l1 < v.size(); l1++) {

			if (edge_sta[l] != edge_end[l1])
			{

				Player[edge_sta[l]].neighbours1.push_back(edge_end[l1]);

				Player[edge_end[l1]].neighbours1.push_back(edge_sta[l]);


			}
		}
	}
	for (int l = 0; l < v.size(); l++) {

		agent = v[l];

		for (int l1 = 0; l1 < v.size(); l1++) {

			if (agent != v[l1]) {
				neigh.push_back(v[l1]);
			}

		}

		Player[v[l]].neighbours.push_back(neigh);
		neigh.clear();

	}
}
vector<vector<int>> node_sum;
void Prodgraph()
{
	ifstream readinEdges;
	readinEdges.open("data_number.csv");

	if (!readinEdges.is_open())
		printf(" can not open real-word data");
	int i;
	string k1;
	long k2;
	int nid;
	vector<int> node_zh;
	vector<string>	pd;
	vector<vector<int>> node_sum;
	node_zh.clear();
	node_sum.clear();


	while (readinEdges)
	{

		readinEdges >> k1 >> k2;

		if (readinEdges) {

			pd.push_back(k1);

			if (pd.size() == 1) {

				node_zh.push_back(k2);

			}
			else {
				if (pd[pd.size() - 2] == pd[pd.size() - 1]) {

					node_zh.push_back(k2);

				}
				else {

					node_sum.push_back(node_zh);
					generate_neigh(node_zh);
					node_zh.clear();
					node_zh.push_back(k2);
					pd.clear();
					pd.push_back(k1);
				}

			}
		}


	}
	generate_neigh(node_zh);
	node_sum.push_back(node_zh);
	readinEdges.close();

	for (int l = 0; l < 402511; l++) {

		Player[l].degree1 = Player[l].neighbours1.size();

	}

}

/***************************End paregraph ****************************************/

void Interact(int x, int z)

{
	
	double sx, so;//possibility whether agent can see defector defects or see an agent violates the metanorm;

	double V;//opportunity of whether an agent punishes the defector.
	
	if (Player[x].boldness > Player[x].S)
	{
		Player[x].DS += T;

		for (int i = 0; i < Player[x].neighbours[z].size(); i++)
		{
			Player[Player[x].neighbours[z][i]].HS += H;

			sx = randf();

			if (sx < Player[x].S)
			{
				V = randf();

				if (V < Player[Player[x].neighbours[z][i]].vengefulness) 
				{

					Player[x].DS += P;
					Player[Player[x].neighbours[z][i]].PS += E;

				}
				else
				{
					for (int k = 0; k < Player[x].neighbours[z].size(); k++)
					{
						if (Player[x].neighbours[z][k] == x || Player[x].neighbours[z][k] == Player[x].neighbours[z][i])

							continue;

						else
						{
							so = randf();

							if (so < Player[Player[x].neighbours[z][i]].S)
							{
								V = randf();

								if (V < Player[Player[x].neighbours[z][k]].vengefulness)
								{

									Player[Player[x].neighbours[z][k]].PS += E;
									Player[Player[x].neighbours[z][i]].POS += P;

								}
							}
						}
					}
				}
			}
		}
	}

}


void Learn(int x)
{
	double Avg = 0;
	double LR;//Learn Rate;

	for (int i = 0; i < Player[x].degree1; i++)
	{

		Avg += Player[Player[x].neighbours1[i]].TS;

	}
	Avg /= Player[x].degree1;

	if (Player[x].TS < Avg) 
	{

		LR = randf();

		if (LR < ER)
		{
			Player[x].boldness = randf();
			Player[x].vengefulness = randf();
		}
		else
		{
			if (Player[x].DS < 0)
			{
				if (Player[x].boldness - LS <= 0)
					Player[x].boldness = 0;

				else

					Player[x].boldness = Player[x].boldness - LS;
			}
			else
			{

				if (Player[x].boldness + LS >= 1)
					Player[x].boldness = 1;

				else
					Player[x].boldness += LS;
			}
			if (Player[x].PS < Player[x].POS)
			{
				if (Player[x].vengefulness - LS <= 0)
					Player[x].vengefulness = 0;

				else
					Player[x].vengefulness = Player[x].vengefulness - LS;
			}
			else
			{
				if (Player[x].vengefulness + LS >= 1)
					Player[x].vengefulness = 1;

				else
					Player[x].vengefulness += LS;
			}
		}
	}
}


void main()
{
	sgenrand((unsigned)time(NULL));
	srand((unsigned)time(NULL));

	int z;
	Prodgraph();

	for (int i = 0; i < TryTime; i++)
	{

		Initial();//³õÊ¼»¯

		for (int k = 0; k < Separate_Run; k++) 
		{

			for (int l = 0; l < 402511; l++)
			{
				if (Player[l].neighbours1.size() != 0) {
					z = randi(Player[l].neighbours.size());
					Interact(l, z);
				}

			}
	
			for (int l = 0; l < 402511; l++)
			{

				Player[l].TS = Player[l].HS + Player[l].DS + Player[l].PS + Player[l].POS;

			}

			for (int l = 0; l < 402511; l++) {

				if (Player[l].neighbours1.size() != 0) {

					Learn(l);

				}

			}

		}

	}

}



