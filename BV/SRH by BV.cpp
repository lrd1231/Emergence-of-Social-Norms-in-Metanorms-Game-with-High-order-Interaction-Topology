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
#define G 5
#define L 4000
#define Separate_Run 1000000 //取1000000次的结果平均（每次实验进行1000000轮博弈）
#define TryTime 10//进行10次独立的实验
#define N 1000 // number of players
#define T 3		//Player’s temptation to defect 背叛获得的收益
#define H -1		//Hurt suffered by others as a result of an agent’s defection 背叛受到的伤害
#define P -9	//Cost of being punished 背叛者收到的惩罚
#define E -2		// Enforcement cost,as known as cost of applying punishment 惩罚成本
#define LS 0.142857		//Learning Step  个体学习率
#define ER 0.01	//Explore Rate 个体探索率
#define PP 0.5

vector<vector<int>> zuhe;
vector<vector<vector<int>>> zuhe_order;
vector<vector<int>>zuhe_order_yuansu1;
vector<int>zuhe_order_yuansu2;

struct Agent
{
	vector<vector<int>> neighbours;//constains players's neighbours
	int neighbours1[N];
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
struct Agent Player[N];

struct node
{
	unsigned long  value;
	int min = 0;
	int max = 0;
	int son_nodes[N];
	int son_sum = 0;
	int number;
	int n;
	int g;
	int yes[N];
	int yes_sum = 0;
	int no[N];
	int no_sum = 0;

};


int sum = 0;
int flag = 0;
std::vector<vector<int>> ZH;
std::vector<node> nodes;
int adj[N][N] = { 0 };

/*************************** paregraph ****************************************/

void make_graph(vector<vector<int>> v) {

	for (int i = 0; i < v.size(); i++) {

		for (int j = 0; j < v[0].size() - 1; j++) {

			for (int k = j + 1; k < v[0].size(); k++) {

				adj[v[i][j]][v[i][k]] = 1;
				adj[v[i][k]][v[i][j]] = 1;

			}
		}
	}
}
ofstream outfile6;

int generateNeighbor() {
	int ln = 0, l = 0;
	for (int i = 0; i < N; i++)
	{
		int m = 0;
		for (int j = 0; j < N; j++) {

			if (adj[i][j] == 1)
			{
				Player[i].neighbours1[m] = j;
				m++;
			}
		}
		Player[i].degree1 = m;
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {

			ln += adj[i][j];

		}

	}
	l = ln / 2;
	/*
	for (int i = 0; i < N; i++)
	{
		cout << "[";
		outfile << "[";
		for (int j = 0; j < N; j++) {
			if (j != N - 1) {
				cout << adj[i][j] << ", ";
				outfile << adj[i][j] << ", ";
			}
			else {
				cout << adj[i][j] << " ";
				outfile << adj[i][j] << " ";
			}

		}
		if (i != N - 1) {
			cout << "],";
			outfile << "],";
		}
		else {
			cout << "]";
			outfile << "]";
		}
		cout << endl;
		outfile << endl;
	}*/

	return l;

}

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
	for (int i = 0; i < N; i++)
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
		Player[i].degree_g = G - 1;
	}
}
void Interact(int x, int z)
{

	double sx, so;//possibility whether agent can see defector defects or see an agent violates the metanorm;
	double V;//opportunity of whether an agent punishes the defector. 

	sx = randf();
	if (Player[x].boldness > sx)
	{
		Player[x].DS += T;

		for (int i = 0; i < Player[x].degree_g; i++)
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
					for (int k = 0; k < Player[x].degree_g; k++)
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
		Player[x].LN++;
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


vector<int> yes;
vector<int> no;

void scalefree_hyper()
{
	double d0, d1;
	double ad0[N], ad1[N];//ad0存放每个节点的权重，ad1存放每个节点的概率分布
	d0 = 0;
	int s;
	for (int i = 0; i < N; i++)
	{
		d0 += pow(i + 1, -PP);
	}
	d1 = 0;
	for (int i = 0; i < N; i++)
	{
		ad0[i] = pow(i + 1, -PP) / d0;
		d1 += ad0[i];
		ad1[i] = d1;
	}

	yes.clear();
	for (int i1 = 0; i1 < G; i1++) {

		d0 = 1;

		while (d0 == 1)
		{
			d0 = randf();
		}
		for (int i2 = 0; i2 < N; i2++)
		{
			if (d0 < ad1[i2])
			{
				s = i2;
				break;
			}
		}
		int b0 = 0;
		for (int i2 = 0; i2 < yes.size(); i2++)
		{
			if (yes[i2] == s)
			{
				b0 = 1;
				break;
			}
		}
		if (b0)
		{
			i1--;
			continue;
		}
		yes.push_back(s);
	}


	/*
	for (int m = 0; m < yes.size(); m++) {
		cout << yes[m] << " ";
	}
	cout << endl;
	*/

}

int vis[N] = { 0 };
vector<vector<int>> nb;

void dfs(int v) {

	//cout << v << " ";
	vis[v] = 1;

	for (int i = 0; i < nb[v].size(); i++) {

		if (!vis[nb[v][i]]) {

			dfs(nb[v][i]);

		}
	}
}
void make_graph_h() {

	int i = 0;
	//k = randi(zuhe_order[l].size());
	for (int l = 0; l < N; l++) {

		for (int k = 0; k < zuhe_order[l].size(); k++) {

			for (int j = 0; j < zuhe_order[l][k].size() - 1; j++) {

				for (int m = j + 1; m < zuhe_order[l][k].size(); m++) {

					adj[zuhe_order[l][k][j]][zuhe_order[l][k][m]] = 1;

					adj[zuhe_order[l][k][m]][zuhe_order[l][k][j]] = 1;

				}
			}
		}
	}
}
void general_neighbors_h() {

	vector<int> neighborhood;
	int k;

	for (int l = 0; l < N; l++) {

		for (int k = 0; k < zuhe_order[l].size(); k++) {

			for (int j = 0; j < zuhe_order[l][k].size(); j++) {

				if (adj[l][zuhe_order[l][k][j]] == 1)
				{
					neighborhood.push_back(zuhe_order[l][k][j]);

				}
			}

			Player[l].neighbours.push_back(neighborhood);

			neighborhood.clear();
		}
	}

}

void main() {
	sgenrand((unsigned)time(NULL));
	srand((unsigned)time(NULL));
	int c = 0, sj = 0, p = 0, q = 0, l = 0, liantong = 100, lnn = 0, z, b0;
	vector<int>yuansu, nb1, zu;
	int i = 0;
	int flag = 0;

	while (liantong == 100) {
		i = 0;
		if (i == 0) {

			scalefree_hyper();

			for (int c = 0; c < G; c++) {

				yuansu.push_back(yes[c]);

			}

			zuhe.push_back(yuansu);
			i++;
			yes.clear();
			no.clear();
			yuansu.clear();
		}
		if (i != 0) {

			while (i < L)
			{

				scalefree_hyper();

				for (int n = 0; n < zuhe.size(); n++) {

					b0 = 0;

					for (int i2 = 0; i2 < yes.size(); i2++)
					{
						for (int c = 0; c < G; c++) {

							if (yes[i2] == zuhe[n][c])
							{
								b0++;

							}
						}

					}
					if (b0 == G) {

						break;

					}

				}
				if (b0 == G)
				{
					continue;
					yes.clear();

				}
				else {
					for (int c = 0; c < G; c++) {

						yuansu.push_back(yes[c]);

					}
					zuhe.push_back(yuansu);
					i++;
					yes.clear();
					yuansu.clear();

				}

			}

		}


		for (int l1 = 0; l1 < N; l1++) {

			memset(adj[l1], 0, sizeof(adj[l1]));

		}
		make_graph(zuhe);

		lnn = generateNeighbor();

		for (int l = 0; l < N; l++) {

			for (int l1 = 0; l1 < Player[l].degree1; l1++) {

				nb1.push_back(Player[l].neighbours1[l1]);

			}

			nb.push_back(nb1);
			nb1.clear();

		}
		dfs(0);

		for (int l = 0; l < N; l++) {
			if (vis[l] == 0) {

				liantong = 100;
				break;

			}
			else {

				liantong = 0;
			}

		}

		for (int l = 0; l < N; l++) {

			vis[l] = 0;

		}

	}


	for (int l = 0; l < L; l++) {
		for (int l1 = 0; l1 < L; l1++) {
			for (int k = 0; k < G; k++) {

				zuhe_order_yuansu2.clear();

				if (zuhe[l1][k] == l) {

					for (int j = 0; j < G; j++) {

						zuhe_order_yuansu2.push_back(zuhe[l1][j]);

					}

					zuhe_order_yuansu1.push_back(zuhe_order_yuansu2);
				}

			}
		}

		zuhe_order.push_back(zuhe_order_yuansu1);
		zuhe_order_yuansu1.clear();
	}

	for (int l1 = 0; l1 < N; l1++) {

		memset(adj[l1], 0, sizeof(adj[l1]));

	}
	make_graph_h();

	general_neighbors_h();

	for (i = 0; i < TryTime; i++)
	{

		Initial();
		double SUMB = 0.0, SUMV = 0.0;
		for (int k = 0; k < Separate_Run; k++) 
		{

			for (l = 0; l < N; l++)
			{
				if (zuhe_order[l].size() != 0) {
					for (int o = 0; o < 4; o++) {

						z = randi(zuhe_order[l].size());
						Interact(l, z);

					}
				}
				else {
					continue;
				}
			}

			for (l = 0; l < N; l++)
			{
				Player[l].TS = Player[l].HS + Player[l].DS + Player[l].PS + Player[l].POS;
			}

			for (l = 0; l < N; l++) {

				Learn(l);
		
			}
			zu.clear();

		}
		
	}

}



