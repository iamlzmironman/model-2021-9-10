#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
using namespace std;

#define L 100     //lattice  size
#define SIZE 10000  //number of players
#define MC_STEPS  10000 //run-time in MCS
#define First_STEPS  1000 //run-time in MCS
#define K  0.1  //noise effects
#define TryTime 2
#define S_steps 0.01
#define RANDOMIZE   220
double a;
double b;

struct Agent
{
	double T;
	short Strat;//contains players' strategies,1 cooperator,0 defector
	int neighbours[4];//constains players's neighbours
	double S;//标记博弈类型：<0囚徒，>0雪堆
	int age;
	double sum_payoff;
	double aver_payoff;
};
struct Agent Player[SIZE];
ofstream outfile, outResults;//outfile 用于输出最后的态，outResults用于输出每次的结果
long int Numsg[2][2];//count the number of four kinds of people
double Result[4][TryTime];
int index[SIZE], shufl[SIZE];
double Cost;
//randf() generate 0-1,randi(x)generate 0-(x-1)
/*************************** RNG procedures ****************************************/
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

/********************** END of RNG ************************************/
void Prodgraph(void)//defines nerghbors on a square lattice
{
	int iu, ju;
	long int player1, player2;
	int i, j;
	for (i = 0; i < L; i++)
	{
		for (j = 0; j < L; j++)
		{
			player1 = L * i + j;       /* consider a site >> a player ,i行j列元素*/
			iu = i + 1;  ju = j;           /* south */
			if (iu == L) iu = 0;
			player2 = L * iu + ju;   /* the location of player2，标记player1的邻居 */
			Player[player1].neighbours[0] = player2;      /* the east link of player1 ends at player 2 */
			iu = i;      ju = j + 1;       /* east */
			if (ju == L) ju = 0;
			player2 = L * iu + ju;
			Player[player1].neighbours[1] = player2;
			iu = i - 1;  ju = j;           /* north*/
			if (i == 0) iu = L - 1;
			player2 = L * iu + ju;
			Player[player1].neighbours[2] = player2;
			iu = i;     ju = j - 1;       /* west */
			if (j == 0) ju = L - 1;
			player2 = L * iu + ju;
			Player[player1].neighbours[3] = player2;
		}//j
	}//i
}//prodgraph

void Initial(void)
{
	int i, j, playerX;
	Numsg[0][0] = 0;
	Numsg[0][1] = 0;
	Numsg[1][0] = 0;
	Numsg[1][1] = 0;	//NumStrat[0] = 0;NumStrat[1] = 0;
	for (i = 0; i < L; i++)
	{
		for (j = 0; j < L; j++)
		{
			playerX = L * i + j;// consider a site >> a player
			Player[playerX].age = 0;
			Player[playerX].sum_payoff = 0;
			Player[playerX].aver_payoff = 0;
			if (randf() <= 0.5)
			{
				Player[playerX].T = b;
				Player[playerX].S = a;
			}
			else
			{
				Player[playerX].T = b;
				Player[playerX].S = -a;
			}
			if (randf() <= 0.5)
			{
				Player[playerX].Strat = 1;
			}
			else
			{
				Player[playerX].Strat = 0;
			}
			if (Player[playerX].S <= 0 && Player[playerX].Strat == 0)//背叛、囚徒
				Numsg[0][0] += 1;
			if (Player[playerX].S <= 0 && Player[playerX].Strat == 1)//合作、囚徒
				Numsg[1][0] += 1;
			if (Player[playerX].S > 0 && Player[playerX].Strat == 0)//背叛、雪堆
				Numsg[0][1] += 1;
			if (Player[playerX].S > 0 && Player[playerX].Strat == 1)//合作、雪堆
				Numsg[1][1] += 1;

		}//j
	}//i

}
/************************新增乱序*************************/
void shuffle(int index[], int aftershuf[]) {
	for (int i = 0; i < SIZE; i++)
		aftershuf[i] = index[i];
	for (int i = 0; i < SIZE; i++)
	{
		int j = rand() % SIZE;
		int temp = aftershuf[i];
		aftershuf[i] = aftershuf[j];
		aftershuf[j] = temp;
	}
}
/****************************end****************************/
double Game(int x)
{
	int y;
	double payoff = 0.0;
	int strat = Player[x].Strat;
	for (int i = 0; i < 4; i++)
	{
		y = Player[x].neighbours[i];
		payoff += (double)strat * Player[y].Strat + Player[x].S * strat * (1 - (double)Player[y].Strat) + Player[x].T * (1 - (double)strat) * Player[y].Strat;
	}
	return payoff;
}

void Learning(int playerX, double payoffX, double payoffY)
{
	double probability = 1 / (1 + exp((payoffX - payoffY) / K));
	if (randf() < probability)
	{
		//NumStrat[Player[playerX].Strat]--;
		if (Player[playerX].S <= 0)
		{
			Numsg[Player[playerX].Strat][0]--;
			Numsg[1 - Player[playerX].Strat][0]++;
		}
		else
		{
			Numsg[Player[playerX].Strat][1]--;
			Numsg[1 - Player[playerX].Strat][1]++;
		}
		Player[playerX].Strat = 1 - Player[playerX].Strat;
		//NumStrat[Player[playerX].Strat]++;
	}

}
void strategyUpdating(int playerX)
{
	int playerY;
	int stratX, stratY;
	double payoffX, payoffY;
	playerY = Player[playerX].neighbours[(int)randi(4)];
	stratX = Player[playerX].Strat;
	stratY = Player[playerY].Strat;
	payoffX = Game(playerX);
	Player[playerX].age += 1;
	Player[playerX].sum_payoff += payoffX;
	Player[playerX].aver_payoff = Player[playerX].sum_payoff / Player[playerX].age;

	if ((stratX + stratY) == 1)
	{

		payoffY = Game(playerY);
		Learning(playerX, payoffX, payoffY);
	}//end of a source selection	
}

/********************************修改博弈竞争阶段****************************************/

void  gameCompetition(int playerSource, int playerTarget)
{
	double payoffSource, payoffTarget, cost;
	//如果选中的是刚刚已经接受过策略输出的节点，直接返回，不进行任何操作
	if (Player[playerSource].age == 0 || Player[playerTarget].age == 0)
		return;
	//获取历史平均收益
	payoffSource = Player[playerSource].aver_payoff;
	payoffTarget = Player[playerTarget].aver_payoff;
	if (Player[playerTarget].S != Player[playerSource].S)//双方认知不同
	{
		cost = Cost;
		double probability = 1 / (1 + exp(((payoffTarget + cost) - payoffSource) / K));
		if (randf() < probability)
		{
			if (Player[playerSource].S <= 0)
			{
				Numsg[Player[playerSource].Strat][0]++;
				Numsg[Player[playerTarget].Strat][1]--;
			}
			else
			{
				Numsg[Player[playerSource].Strat][1]++;
				Numsg[Player[playerTarget].Strat][0]--;
			}
			//博弈输出
			Player[playerTarget].T = Player[playerSource].T;
			Player[playerTarget].S = Player[playerSource].S;
			Player[playerTarget].Strat = Player[playerSource].Strat;//策略推广
			Player[playerTarget].age = 0;
			Player[playerTarget].sum_payoff = 0;
			Player[playerTarget].aver_payoff = 0;
		}
	}
}

/****************************************博弈竞争修改结束**********************************************/

double AveragePI(int g)
{
	double x = 0.0;
	for (int i = 0; i < TryTime; i++)
	{
		x += Result[g][i];

	}
	x = x / TryTime;
	return x;
}
void ResultOut(double pi)
{
	outfile << b << "\t" << a << "\t" << pi << endl;
	cout << b << "\t" << a << "\t" << pi << endl;
}
int main()
{
	//sgenrand((unsigned)time(NULL));//initialize RNG
	sgenrand(RANDOMIZE); // initialize the RNG
	int steps, i, playerX, playerSource, playerTarget;
	int test;
	double P_dp, P_ds, P_cp, P_cs;
	double AverP_dp, AverP_ds, AverP_cp, AverP_cs, AverP_c;
	outResults.open("Results__研究大小cost=0.xls", ios::app);
	outfile.open("data__研究大小cost=0.xls");
	/****************new*******************/
	//初始化index数组
	for (int i = 0; i < SIZE; i++)
		index[i] = i;
	/*****************end********************/
	//初始化个体的邻居
	Prodgraph();
	outResults << "花费" << "\t" << "|S|" << '\t' << "T" << "\t" << "背叛-囚徒" << '\t' << "背叛-残雪" << '\t' << "合作-背叛" << '\t' << "合作-残雪" << '\t' << "合作者" << '\t' << "共存" << endl;
	cout << "花费" << "\t" << "|S|" << '\t' << "T" << "\t" << "背叛-囚徒" << '\t' << "背叛-残雪" << '\t' << "合作-背叛" << '\t' << "合作-残雪" << '\t' << "合作者" << '\t' << "共存" << endl;
	for (Cost = 0; Cost <= 0; Cost += 1)//cost取0,0.5,1,2,3,4,10
	{
		for (b = 1; b <= 2.01; b += 0.01)
		{
			for (a = 0; a <= 1.01; a += 0.01)//θ
			{
				for (test = 0; test < TryTime; test++)
				{
					Initial();//inital strategy distribution
					for (steps = 0; steps < MC_STEPS; steps++)
					{
						/*****************************模型修改****************************/

						//所有player按照乱序顺序（不重复随机选择节点）进行策略更新
						shuffle(index, shufl);
						for (i = 0; i < SIZE; i++)
						{
							playerX = shufl[i];//choose a source site
							strategyUpdating(playerX);
						}
						//不重复随机选择节点，对任意一个邻居进行价值观输出
						shuffle(index, shufl);
						for (i = 0; i < SIZE; i++)
						{
							playerSource = shufl[i];//源节点，输出博弈价值观的节点
							playerTarget = Player[playerSource].neighbours[(int)randi(4)];//随机选择一个邻居节点作为目标节点，接受价值观输出
							gameCompetition(playerSource, playerTarget);
						}//end of elementary MC step

						/***********************模型修改结束***********************/
						if (Numsg[0][0] + Numsg[1][0] == 0 || Numsg[0][1] + Numsg[1][1] == 0)
						{
							break;
						}
					}//end of steps
					P_dp = (double)Numsg[0][0] / SIZE;
					P_ds = (double)Numsg[0][1] / SIZE;
					P_cp = (double)Numsg[1][0] / SIZE;
					P_cs = (double)Numsg[1][1] / SIZE;
					Result[0][test] = P_dp;
					Result[1][test] = P_ds;
					Result[2][test] = P_cp;
					Result[3][test] = P_cs;
				}//end of test
				AverP_dp = AveragePI(0);
				AverP_ds = AveragePI(1);
				AverP_cp = AveragePI(2);
				AverP_cs = AveragePI(3);
				AverP_c = AverP_cp + AverP_cs;
				if (AverP_dp > 0 && AverP_ds > 0 && AverP_cp > 0 && AverP_cs > 0)
				{
					outResults << Cost << '\t' << a << '\t' << b << "\t" << AverP_dp << '\t' << AverP_ds << '\t' << AverP_cp << '\t' << AverP_cs << '\t' << AverP_c << '\t' << 1 << endl;
					cout << Cost << '\t' << a << '\t' << b << "\t" << AverP_dp << '\t' << AverP_ds << '\t' << AverP_cp << '\t' << AverP_cs << '\t' << AverP_c << '\t' << 1 << endl;
				}
				else
				{
					outResults << Cost << '\t' << a << '\t' << b << "\t" << AverP_dp << '\t' << AverP_ds << '\t' << AverP_cp << '\t' << AverP_cs << '\t' << AverP_c << '\t' << 0 << endl;
					cout << Cost << '\t' << a << '\t' << b << "\t" << AverP_dp << '\t' << AverP_ds << '\t' << AverP_cp << '\t' << AverP_cs << '\t' << AverP_c << '\t' << 0 << endl;
				}
			}//end of a
		}
	}
	outResults.close();
	outfile.close();
	return 0;
}