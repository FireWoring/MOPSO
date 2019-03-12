#ifndef _MOPSO_H_
#define _MOPSO_H_

#include "canshu.h"
#include "Coordinate.h"
#include "Particle.h"
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <math.h>
#include <iterator>

#define MAXITERATION 10  //迭代次数
#define PARTICLENUMM  40  //粒子群粒子数
#define REPNUM       40  //rep的粒子数
#define GRIDNUM      10   //栅格个数
#define DIMENDION    2    //栅格维度

using namespace std;

class MOPSO
{
public:
	vector<Particle> rep;
	vector<Particle> pop;
	Particle leader;

	MOPSO();
	void getCost(vector<pair<double, double>> &cost, vector<Particle>& p);
	void domination(Particle &x, Particle &y, int &flag);
	void determineDomination(vector<Particle>& p);
	void addTORep();
	void createGrid(vector<Particle>& p);
	void findGridIndex(Particle &particle);
	int  rouletteWheelSelection(vector<double> &p);
	void selectLeader();
	void updateVolicityAndPosition(Particle &p);
	Coordinate mutate(Coordinate x, double pm);
	void deleteOneRepMember();
	void MOPSOFunc();
private:
	static int Cmax, Cmin;
	static double Umax, Umin;
	static double AlphaHmax, AlphaHmin;
	static double AlphaCmax, AlphaCmin;
	static double BetaHmax, BetaHmin;
	static double BetaCmax, BetaCmin;
	static double Vcmax, Vcmin;
	static double Vumax, Vumin;
	static double VAlphaHmax, VAlphaHmin;
	static double VAlphaCmax, VAlphaCmin;
	static double VBetaHmax, VBetaHmin;
	static double VBetaCmax, VBetaCmin;
	static double c1, c2;//学习因子
	static double w;
	static double wdamp; //惯性因子下降速率
	static double alpha; //Inflation rate

	//栅格的横纵坐标刻度数组
	vector<double> x;
	vector<double> y;

	//leader选择压力
	static int beta;
	static int gamma;

	//突变率
	static double mu;

	Particle *p[PARTICLENUMM];
};

#endif