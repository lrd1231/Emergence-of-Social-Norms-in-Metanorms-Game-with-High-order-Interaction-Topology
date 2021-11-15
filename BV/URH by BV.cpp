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
#define G 2
#define L 10000
#define Separate_Run 1000000 
#define TryTime 10
#define N 1000// number of players
#define T 3		//Player’s temptation to defect 
#define H -1		//Hurt suffered by others as a result of an agent’s defection 
#define P -9	//Cost of being punished 
#define E -2		// Enforcement cost,as known as cost of applying punishment 
#define LS 0.142857		//Learning Step  
#define ER 0.01	//Explore Rate 


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
//node nodes[1000000];
int flag = 0;
std::vector<vector<int>> ZH;
std::vector<node> nodes;
int adj[N][N] = { 0 };




/*************************** paregraph ****************************************/

void make_graph(vector<vector<int>> v) {

	for (int i = 0; i < v.size(); i++) {
		for (int j = 0; j < v[0].size() - 1; j++) {
			for (int k = j + 1; k < v[0].size(); k++) {

				adj[v[i][j] - 1][v[i][k] - 1] = 1;
				adj[v[i][k] - 1][v[i][j] - 1] = 1;

			}
		}
	}
}


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
void Initial(void)//初始化函数，形成网格结构中每个主体的各个参数初始状态值
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
		Player[i].Punish = 0;//
		Player[i].learn = 0;
		Player[i].LR = 0;
		Player[i].degree_g = G - 1;

	}
}

void Interact(int x, int z)

{

	double sx, so;//possibility whether agent can see defector defects or see an agent violates the metanorm;
	double V;//opportunity of whether an agent punishes the defector. 

	if (Player[x].boldness > Player[x].S)//
	{
		
		Player[x].DS += T;//
		for (int i = 0; i < Player[x].degree_g; i++)
		{
			Player[Player[x].neighbours[z][i]].HS += H;//
			
			sx = randf();
			if (sx < Player[x].S)
			{
				V = randf();//
				if (V < Player[Player[x].neighbours[z][i]].vengefulness) //
				{
					Player[x].DS += P; //
					Player[Player[x].neighbours[z][i]].PS += E;//
				}
				else//
				{
		
					for (int k = 0; k < Player[x].degree_g; k++)
					{
						if (Player[x].neighbours[z][k] == x || Player[x].neighbours[z][k] == Player[x].neighbours[z][i]) // 
							continue;
						else
						{
							so = randf();

							if (so < Player[Player[x].neighbours[z][i]].S)//
							{
								V = randf();

								if (V < Player[Player[x].neighbours[z][k]].vengefulness)//
								{

									Player[Player[x].neighbours[z][k]].PS += E;//
									Player[Player[x].neighbours[z][i]].POS += P;//

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
	for (int i = 0; i < Player[x].degree1; i++)//
	{
		Avg += Player[Player[x].neighbours1[i]].TS;
	}
	Avg /= Player[x].degree1;
	if (Player[x].TS < Avg) //
	{
		Player[x].LN++;
		LR = randf();
		if (LR < ER)//以探索率ER在学习
		{
			Player[x].boldness = randf();
			Player[x].vengefulness = randf();
		}
		else
		{
			if (Player[x].DS < 0)//
			{
				if (Player[x].boldness - LS <= 0)//  
					Player[x].boldness = 0;
				else
					Player[x].boldness = Player[x].boldness - LS;
			}
			else//
			{
				if (Player[x].boldness + LS >= 1)
					Player[x].boldness = 1;
				else
					Player[x].boldness += LS;
			}
			if (Player[x].PS < Player[x].POS)//
			{
				if (Player[x].vengefulness - LS <= 0)
					Player[x].vengefulness = 0;
				else
					Player[x].vengefulness = Player[x].vengefulness - LS;
			}
			else//
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
void new_one() {
	int da = N;
	int xiao = G;
	double x, board;

	int last_number = 0;
	for (int i = 1; i <= G; i++) {
		if ((xiao >= da) || (xiao == 1))
		{
			break;
		}

		double random_number = randf();

		for (int j = 1; j <= da - xiao + 1; j++) {

			x = j / (double(da - xiao) + 1);
			board = 1 - pow(1 - x, xiao);

			if (random_number <= board) {

				int choose_son = j;

				for (int k = 1; k < j; k++) {

					no.push_back(last_number + k);

				}

				last_number = last_number + j;
				yes.push_back(last_number);

				da = da - j;
				xiao = xiao - 1;

				break;
			}
		}
	}

	for (int i = 1; i <= G - yes.size(); i++) {

		int random_int = (rand() % (N - last_number)) + 1;
		yes.push_back(random_int + last_number);

	}

	//for (int m = 0; m < yes.size(); m++) {
	//	cout << yes[m]<<" ";
	//}
	//cout << endl;


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

					adj[zuhe_order[l][k][j] - 1][zuhe_order[l][k][m] - 1] = 1;
					adj[zuhe_order[l][k][m] - 1][zuhe_order[l][k][j] - 1] = 1;

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

				if (adj[l][zuhe_order[l][k][j] - 1] == 1)
				{
					neighborhood.push_back(zuhe_order[l][k][j] - 1);

				}
			}

			Player[l].neighbours.push_back(neighborhood);

			neighborhood.clear();
		}
	}
}

/********************** END of RNG ************************************/

void main()
{
	sgenrand((unsigned)time(NULL));
	srand((unsigned)time(NULL));
	int c = 0,  p = 0, q = 0,  liantong = 100, lnn = 0, z, h = 22, flag = 0;
	vector<int>yuansu, nb1, zu;
	double i = 0;

	while (liantong == 100) {

		nb.clear();
		nb1.clear();
		zuhe.clear();
		i = 0;

		if (i == 0) {

			new_one();

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
				new_one();

				for (int n = 0; n < zuhe.size(); n++) {

					p = 0; q = 0; flag = 0;

					while (zuhe[n][p] == yes[q] && q < yes.size()) {

						flag++;
						q++;
						p++;
						if (p >= G || q >= G) {
							break;
						}

					}
					if (flag == G) {

						yes.clear();
						no.clear();
						yuansu.clear();
						break;

					}

				}
				if (flag != G) {

					for (int c = 0; c < G; c++) {

						yuansu.push_back(yes[c]);

					}
					zuhe.push_back(yuansu);
					i++;
					yes.clear();
					no.clear();
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

	for (int l = 0; l < N; l++) {
		for (int l1 = 0; l1 < L; l1++) {
			for (int k = 0; k < G; k++) {

				zuhe_order_yuansu2.clear();

				if (zuhe[l1][k] == l + 1) {

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
	
		for (int k = 0; k < Separate_Run; k++) 
		{
	
			for (int l = 0; l < N; l++)
			{
				for (int o = 0; o < 4; o++) {

					z = randi(zuhe_order[l].size());

					Interact(l, z);

				}

			}
			/*计算每个个体的总得分TS*/


			for (int l = 0; l < N; l++)
			{

				Player[l].TS = Player[l].HS + Player[l].DS + Player[l].PS + Player[l].POS;

			}
			
			for (int l = 0; l < N; l++) {

				Learn(l);

			}

			


			

		}


	}
}



