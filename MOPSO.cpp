#define _CRT_SECURE_NO_WARNINGS
#include "MOPSO.h"
#include <Windows.h>
#include <fstream>

int MOPSO::Cmin = 1;
int MOPSO::Cmax = 15;
double MOPSO::Umin = 2.0;
double MOPSO::Umax = 3.0;
double MOPSO::AlphaHmin = 3.0;
double MOPSO::AlphaHmax = 4.5;
double MOPSO::AlphaCmin = 1.0;
double MOPSO::AlphaCmax = 3.0;
double MOPSO::BetaHmin = 5.0;
double MOPSO::BetaHmax = 8.0;
double MOPSO::BetaCmin = 3.0;
double MOPSO::BetaCmax = 4.5;
double MOPSO::Vcmax = MOPSO::Cmax - MOPSO::Cmin;
double MOPSO::Vcmin = 0 - MOPSO::Vcmax;
double MOPSO::Vumax = MOPSO::Umax - MOPSO::Umin;
double MOPSO::Vumin = 0 - MOPSO::Vumax;
double MOPSO::VAlphaHmax = MOPSO::AlphaHmax - MOPSO::AlphaHmin;
double MOPSO::VAlphaHmin = 0 - MOPSO::VAlphaHmax;
double MOPSO::VAlphaCmax = MOPSO::AlphaCmax - MOPSO::AlphaCmin;
double MOPSO::VAlphaCmin = 0 - MOPSO::VAlphaCmax;
double MOPSO::VBetaHmax = MOPSO::BetaHmax - MOPSO::BetaHmin;
double MOPSO::VBetaHmin = 0 - MOPSO::VBetaHmax;
double MOPSO::VBetaCmax = MOPSO::BetaCmax - MOPSO::BetaCmin;
double MOPSO::VBetaCmin = 0 - MOPSO::VBetaCmax;

double MOPSO::c1 = 1.0;
double MOPSO::c2 = 1.0;
double MOPSO::w = 0.4;

double MOPSO::wdamp = 0.99;
double MOPSO::alpha = 0.1;

int MOPSO::beta = 2;
int MOPSO::gamma = 2;
double MOPSO::mu = 0.1;

MOPSO::MOPSO()
{
	int randC;
	double randU, randAlphaH, randAlphaC, randBetaH, randBetaC;
	srand((unsigned)time(NULL));
	for (int i = 0; i < PARTICLENUMM; i++)
	{
		randC = (int)((rand() / (double)RAND_MAX)*(Cmax - Cmin)) + Cmin;
		randU = ((rand() / (double)RAND_MAX)*(Umax - Umin)) + Umin;
		randAlphaH = ((rand() / (double)RAND_MAX)*(AlphaHmax - AlphaHmin)) + AlphaHmin;
		randAlphaC = ((rand() / (double)RAND_MAX)*(AlphaCmax - AlphaCmin)) + AlphaCmin;
		randBetaH = ((rand() / (double)RAND_MAX)*(BetaHmax - BetaHmin)) + BetaHmin;
		randBetaC = ((rand() / (double)RAND_MAX)*(BetaCmax - BetaCmin)) + BetaCmin;

		canshu can(NNUM, randC, randU, LAMDA, randAlphaH, randAlphaC, randBetaH, randBetaC);

		p[i] = new Particle(randC, randU, randAlphaH, randAlphaC, randBetaH, randBetaC);
		pop.push_back(*(new Particle(randC, randU, randAlphaH, randAlphaC, randBetaH, randBetaC)));
		/*cout << "c:" << pop[i].getCorrdinate().c << " u:" << pop[i].getCorrdinate().u << 
			" AlphaH:" << pop[i].getCorrdinate().alphaH << " AlphaC:" << pop[i].getCorrdinate().alphaC <<
			" BetaH:" << pop[i].getCorrdinate().betaH << " BetaC:" << pop[i].getCorrdinate().betaC << endl;*/
	}
}

void MOPSO::getCost(vector<pair<double, double>> &results, vector<Particle>& p)
{
	pair<double, double> tmp;
	//计算每一个粒子的函数值，存在pair数组中
	for (int i = 0; i < p.size(); ++i)
	{
		tmp.first = p.at(i).f->F(p.at(i).answer);
		tmp.second = p.at(i).f->TNom(p.at(i).answer);
		results.push_back(tmp);
		//cout << results[i].first << " " << results[i].second << endl;
	}
}

void MOPSO::domination(Particle &x, Particle &y, int &flag)
{
	//flag=0说明两者互不支配
	//flag=1说明x支配y
	//flag=2说明y支配x
	if (!(x.f->F(x.answer)) <= y.f->F(y.answer) && x.f->TNom(x.answer) <= y.f->TNom(y.answer))
		flag = 2;
	else if ((!(x.f->F(x.answer) >= y.f->F(y.answer) && x.f->TNom(x.answer) >= y.f->TNom(y.answer))))
		flag = 1;
	else
		flag = 0;
}

void MOPSO::determineDomination(vector<Particle>& p)
{
	vector<pair<double, double>> results;
	vector<int> hasDominate;
	getCost(results, p);
	//将支配解的函数值置为0
	for (int i = 0; i < p.size() - 1; ++i)
	{
		for (int j = i; j < p.size(); ++j)
		{
			if (results[i].first >= results[j].first&&results[i].second > results[j].second || results[i].first > results[j].first&&results[i].second >= results[j].second)
				hasDominate.push_back(i);

			if (results[j].first >= results[i].first&&results[j].second > results[i].second || results[j].first > results[i].first&&results[j].second >= results[i].second)
				hasDominate.push_back(j);
		}
	}
	//将非支配解放入到REF中
	//cout << p.size() << endl;
	unique(hasDominate.begin(), hasDominate.end());
	for (size_t i = 0; i < p.size(); ++i)
	{
		if (find(hasDominate.begin(), hasDominate.end(), i) == hasDominate.end())
		{
			vector<Particle>::const_iterator ite = rep.begin();
			bool findFlag = false;
			for (; ite != rep.end(); ++ite)
			{
				if (ite->getCorrdinate() == p.at(i).getCorrdinate())
				{
					findFlag = true;
					break;
				}
			}
			if (findFlag == false)
				rep.push_back(*(new Particle(p.at(i).f->class_can.c, p.at(i).f->class_can.u, p.at(i).f->class_can.alphaH, p.at(i).f->class_can.alphaC, p.at(i).f->class_can.betaH, p.at(i).f->class_can.betaC)));
		}
	}

	/////////////////////测试代码
	/*vector<Particle>::iterator i;
	int k = 0;
	for (i = rep.begin(); i != rep.end(); ++i)
	{
	cout << k++ << endl;
	cout << "D:" << i->getCorrdinate().d << " N:" << i->getCorrdinate().N << " K:" << i->getCorrdinate().K << " Ub:" << i->getCorrdinate().ub << " Uv:" << i->getCorrdinate().uv << endl;
	}*/
	//////////////////////////
}

void MOPSO::addTORep()
{
	//计算更新之后pop的非支配解
	determineDomination(pop);

	//将rep拷贝到一个缓存空间，求新的rep
	vector<Particle> tmp(rep.begin(), rep.end());
	rep.clear();
	//cout << "tmp size():" << tmp.size() << endl;
	determineDomination(tmp);
	//cout << "rep size    ():" << rep.size() << endl;
}

void MOPSO::createGrid(vector<Particle>& p)
{
	vector<pair<double, double>> results;
	getCost(results, p);

	//求最大最小的函数值，来栅格坐标的最大最小值
	pair<double, double> minResult, maxResult;
	minResult.first = results[0].first;
	minResult.second = results[0].second;
	maxResult.first = results[0].first;
	maxResult.second = results[0].second;
	for (int i = 1; i < PARTICLENUMM; ++i)
	{
		//求最小的坐标
		if (minResult.first > results[i].first)
			minResult.first = results[i].first;

		if (minResult.second > results[i].second)
			minResult.second = results[i].second;

		//求最大的坐标
		if (maxResult.first < results[i].first)
			maxResult.first = results[i].first;

		if (maxResult.second < results[i].second)
			maxResult.second = results[i].second;
	}

	minResult.first -= alpha*(maxResult.first - minResult.first);
	minResult.second -= alpha*(maxResult.second - minResult.second);
	maxResult.first += alpha*(maxResult.first - minResult.first);
	maxResult.second += alpha*(maxResult.second - minResult.second);

	//计算栅格的横纵坐标刻度
	double everX = (maxResult.first - minResult.first) / 10;
	double everY = (maxResult.second - minResult.second) / 10;
	//cout << "x坐标每一个刻度的长度：" << everX << " y坐标每一个刻度的长度：" << everY << endl;
	x.resize(0);
	y.resize(0);
	for (int i = 0; i < GRIDNUM; ++i)
	{
		x.push_back(minResult.first + i*everX);
		y.push_back(maxResult.first + i*everY);
	}
}

void MOPSO::findGridIndex(Particle & particle)
{
	//初始化粒子坐标
	int flag = 0;

	for (int i = 0; i < x.size(); ++i)
	{
		if (particle.f->F(particle.answer) < x.at(i))
		{
			particle.gridSubIndex.x = i;
			break;
		}
	}

	for (int i = 0; i < y.size(); ++i)
	{
		if (particle.f->TNom(particle.answer) < y.at(i))
		{
			particle.gridSubIndex.y = i;
			break;
		}
	}

	//计算粒子索引
	particle.gridIndex = (particle.gridSubIndex.x + 1)* (GRIDNUM + 2) + particle.gridSubIndex.y;
}

int MOPSO::rouletteWheelSelection(vector<double>& p)
{
	/*for (size_t i = 0; i < p.size(); ++i)
	{
	cout << p.at(i) << " ";
	}
	cout << endl;*/

	int i;
	double r = rand() / (double)RAND_MAX;
	vector<double> C;
	if (p.size() > 0)
	{
		C.push_back(*p.begin());

		for (i = 1; i < p.size(); ++i)
			C.push_back(*(p.begin() + i) + C.at(i - 1));

		for (i = 0; i < p.size(); ++i)
			if (r <= *(C.begin() + i))
				break;

		return i;
	}
	return -1;
}

void MOPSO::selectLeader()
{
	volatile size_t i;
	//所有存储库成员的网格索引
	vector<int> GI;
	for (i = 0; i < rep.size(); ++i)
		GI.push_back(rep[i].gridIndex);
	//将网格索引归一
	unique(GI.begin(), GI.end());

	//构建索引值对应的索引数组
	int* indexNum = new int[GI.size()];
	for (i = 0; i < GI.size(); ++i)
		indexNum[i] = 0;
	for (i = 0; i < GI.size(); ++i)
	{
		for (size_t j = 0; j < rep.size(); ++j)
		{
			if (GI.at(i) == rep[j].gridIndex)
			{
				++indexNum[i];
			}
		}
	}

	// 选择概率
	vector<double> p;
	double sum = 0.0;
	for (i = 0; i < GI.size(); ++i)
	{
		double tmp = exp((double)(-beta*indexNum[i]));
		p.push_back(tmp);
		sum += tmp;
	}
	for (i = 0; i < GI.size(); ++i)
	{
		p.at(i) /= sum;
	}

	delete[] indexNum;
	//indexNum = NULL;

	//选择相应的粒子在GI中的索引
	int index = rouletteWheelSelection(p);

	//从索引数组中找到选择粒子所在的栅格索引
	if (index == -1)
		return;
	int GIindex = GI.at(index);

	//取出栅格索引对应的粒子群
	vector<int> selectIndex;
	for (i = 0; i < rep.size(); ++i)
	{
		if (GIindex == rep.at(i).gridIndex)
			selectIndex.push_back(i);
	}

	//从上述粒子群中任意选择一个
	int particleIndex;
	if (selectIndex.size() == 0)
	{
		particleIndex = 0;
	}
	else
	{
		particleIndex = rand() % selectIndex.size();
	}

	//取出对应的粒子返回
	leader = rep.at(selectIndex[particleIndex]);
}

void MOPSO::updateVolicityAndPosition(Particle &p)
{
	//更新速度
	p.Vc = (int)(w*p.Vc + (double)c1*rand() / RAND_MAX*(p.getPbest().c - p.getCorrdinate().c) + (double)c2*rand() / RAND_MAX*(leader.getPbest().c - p.getCorrdinate().c));
	p.Vu = w*p.Vu + c1*rand() / (double)RAND_MAX*(p.getPbest().u - p.getCorrdinate().u) + c2*rand() / (double)RAND_MAX*(leader.getPbest().u - p.getCorrdinate().u);
	p.ValphaH = w*p.ValphaH + c1*rand() / (double)RAND_MAX*(p.getPbest().alphaH - p.getCorrdinate().alphaH) + c2*rand() / (double)RAND_MAX*(leader.getPbest().alphaH - p.getCorrdinate().alphaH);
	p.ValphaC = w*p.ValphaC + c1*rand() / (double)RAND_MAX*(p.getPbest().alphaC - p.getCorrdinate().alphaC) + c2*rand() / (double)RAND_MAX*(leader.getPbest().alphaC - p.getCorrdinate().alphaC);
	p.VbetaH = w*p.VbetaH + c1*rand() / (double)RAND_MAX*(p.getPbest().betaH - p.getCorrdinate().betaH) + c2*rand() / (double)RAND_MAX*(leader.getPbest().betaH - p.getCorrdinate().betaH);
	p.VbetaC = w*p.VbetaC + c1*rand() / (double)RAND_MAX*(p.getPbest().betaC - p.getCorrdinate().betaC) + c2*rand() / (double)RAND_MAX*(leader.getPbest().betaC - p.getCorrdinate().betaC);

	//更新位置
	int newC = min(max((int)(p.getCorrdinate().c + p.Vc), Cmin), Cmax);
	double newU = p.getCorrdinate().u + p.Vu > Umin ? p.getCorrdinate().u + p.Vu : Umin;
	newU = newU < Umax ? newU : Umax;
	double newAlphaH = p.getCorrdinate().alphaH + p.ValphaH > AlphaHmin ? p.getCorrdinate().alphaH + p.ValphaH : AlphaHmin;
	newAlphaH = newAlphaH < AlphaHmax ? newAlphaH : AlphaHmax;
	double newAlphaC = p.getCorrdinate().alphaC + p.ValphaC > AlphaCmin ? p.getCorrdinate().alphaC + p.ValphaC : AlphaCmin;
	newAlphaC = newAlphaC < AlphaCmax ? newAlphaC : AlphaCmax;
	double newBetaH = p.getCorrdinate().betaH + p.VbetaH > BetaHmin ? p.getCorrdinate().betaH + p.VbetaH : BetaHmin;
	newBetaH = newBetaH < BetaHmax ? newBetaH : BetaHmax;
	double newBetaC = p.getCorrdinate().betaC + p.VbetaC > BetaCmin ? p.getCorrdinate().betaC + p.VbetaC : BetaCmin;
	newBetaC = newBetaC < BetaCmax ? newBetaC : BetaCmax;
	p.setCoordinate(newC, newU, newAlphaH, newAlphaC, newBetaH, newBetaC);
	//cout << newC << " " << newU << " " << newAlphaH << " " << newAlphaC << " " << newBetaH << " " << newBetaC << endl;
}

Coordinate MOPSO::mutate(Coordinate x, double pm)
{
	//6代表自变量是六维的
	int j = rand() % 6 + 1;

	//变异范围
	int newC = (int)(pm*(Cmax - Cmin));
	double newU = (double)pm*(Umax - Umin);
	double newAlphaH = (double)pm*(AlphaHmax - AlphaHmin);
	double newAlphaC = (double)pm*(AlphaCmax - AlphaCmin);
	double newBetaH = (double)pm*(BetaHmax - BetaHmin);
	double newBetaC = (double)pm*(BetaCmax - BetaCmin);

	int lbInt, ubInt;
	double lbDouble = 0.0, ubDouble = 0.0;
	//cout << pm << endl;
	//cout << x.d << " " << x.N << " " << x.K << " " << x.ub << " " << x.uv << endl;
	switch (j)
	{
	case 1:
		lbInt = x.c - newC;
		ubInt = x.c + newC;
		if (lbInt < Cmin) lbInt = Cmin;
		if (ubInt > Cmax) ubInt = Cmax;
		x.c = (int)(rand() % (ubInt - lbInt)) + lbInt;
		break;
	case 2:
		lbDouble = x.u - newU;
		ubDouble = x.u + newU;
		if (lbDouble < Umin) lbDouble = Umin;
		if (ubDouble > Umax) ubDouble = Umax;
		x.u = (rand() / (double)RAND_MAX) * (ubDouble - lbDouble) + lbDouble;
		break;
	case 3:
		lbDouble = x.alphaH - newAlphaH;
		ubDouble = x.alphaH + newAlphaH;
		if (lbDouble < AlphaHmin) lbDouble = AlphaHmin;
		if (ubDouble > AlphaHmax) ubDouble = AlphaHmax;
		x.alphaH = (rand() / (double)RAND_MAX) * (ubDouble - lbDouble) + lbDouble;
		break;
	case 4:
		lbDouble = x.alphaC - newAlphaC;
		ubDouble = x.alphaC + newAlphaC;
		if (lbDouble < AlphaCmin) lbDouble = AlphaCmin;
		if (ubDouble > AlphaCmax) ubDouble = AlphaCmax;
		x.alphaC = (rand() / (double)RAND_MAX) * (ubDouble - lbDouble) + lbDouble;
		break;
	case 5:
		lbDouble = x.betaH - newBetaH;
		ubDouble = x.betaH + newBetaH;
		if (lbDouble < BetaHmin) lbDouble = BetaHmin;
		if (ubDouble > BetaHmax) ubDouble = BetaHmax;
		x.betaH = (rand() / (double)RAND_MAX) * (ubDouble - lbDouble) + lbDouble;
		break;
	case 6:
		lbDouble = x.betaC - newBetaC;
		ubDouble = x.betaC + newBetaC;
		if (lbDouble < BetaCmin) lbDouble = BetaCmin;
		if (ubDouble > BetaCmax) ubDouble = BetaCmax;
		x.betaC = (rand() / (double)RAND_MAX) * (ubDouble - lbDouble) + lbDouble;
		break;
	}
	//cout << x.c <<" " <<x.u <<" "<<x.alphaH<<" "<<x.alphaC<<" "<<x.betaH <<" "<<x.betaC << endl;
	return x;
}

void MOPSO::deleteOneRepMember()
{
	volatile size_t i;
	//所有存储库成员的网格索引
	vector<int> GI;
	for (i = 0; i < rep.size(); ++i)
		GI.push_back(rep[i].gridIndex);
	//将网格索引归一
	unique(GI.begin(), GI.end());

	//构建索引值对应的索引数组
	int* indexNum = new int[GI.size()];
	for (i = 0; i < GI.size(); ++i)
		indexNum[i] = 0;
	for (i = 0; i < GI.size(); ++i)
	{
		for (size_t j = 0; j < rep.size(); ++j)
		{
			if (GI.at(i) == rep[j].gridIndex)
			{
				++indexNum[i];
			}
		}
	}

	// 选择概率
	vector<double> p;
	double sum = 0.0;
	for (i = 0; i < GI.size(); ++i)
	{
		double tmp = exp((double)(gamma*indexNum[i]));
		p.push_back(tmp);
		sum += tmp;
	}
	for (i = 0; i < GI.size(); ++i)
	{
		p.at(i) /= sum;
	}

	delete[] indexNum;
	//indexNum = NULL;

	//选择相应的粒子在GI中的索引
	int index = rouletteWheelSelection(p);

	//从索引数组中找到选择粒子所在的栅格索引
	if (index == -1)
		return;
	int GIindex = GI.at(index);

	//取出栅格索引对应的粒子群
	vector<int> selectIndex;
	for (size_t i = 0; i < rep.size(); ++i)
	{
		if (GIindex == rep.at(i).gridIndex)
			selectIndex.push_back(i);
	}

	//从上述粒子群中任意选择一个
	srand((unsigned int)time(NULL));
	int particleIndex = selectIndex.at(rand() % selectIndex.size());

	//删除对应的粒子
	rep.erase(rep.begin() + particleIndex);
}

void MOPSO::MOPSOFunc()
{
	//计算是否支配，并将非支配解放置到Rep当中
	determineDomination(pop);
	//cout << "rep size():" << rep.size() << endl;

	//创建栅格，然后初始化粒子的坐标和索引
	createGrid(pop);
	for (size_t i = 0; i < rep.size(); ++i)
	{
		findGridIndex(rep[i]);
		/*cout << "第i个粒子：" << i << endl;
		cout << "x,y坐标：" << rep[i].gridSubIndex.x << " " << rep[i].gridSubIndex.y << endl;
		cout << "粒子索引：" << rep[i].gridIndex << endl << endl;*/
	}

	//main loop
	for (int i = 0; i < MAXITERATION; ++i)
	{
		for (int j = 0; j < PARTICLENUMM; ++j)
		{
			cout << "mopso main loop : " << i << " " << j << endl;
			//寻找leader
			selectLeader();
			//更新速度位置
			updateVolicityAndPosition(pop.at(j));
			//使用突变
			double pm = pow((1 - i / MAXITERATION), (1 / mu));
			Particle newSol;
			int flag;
			srand((unsigned int)time(NULL));
			if ((rand() / (double)RAND_MAX) < pm)
			{
				Coordinate tmp = mutate(pop.at(j).getCorrdinate(), pm);
				newSol.setCoordinate(tmp);
				//cout << "xxx" << newSol.getCorrdinate().d << " " << newSol.getCorrdinate().N << " " << newSol.getCorrdinate().K << " " << newSol.getCorrdinate().ub << " " << newSol.getCorrdinate().uv << endl;
				domination(newSol, pop.at(j), flag);
				//cout << "xxxxxxxxxx" << endl;
				if (flag == 1)
					pop.at(j).setCoordinate(newSol.getCorrdinate());
				else if (flag == 0)
				{
					if (rand() / RAND_MAX < 0.5)
						pop.at(j).setCoordinate(newSol.getCorrdinate());
				}
				//cout << "1: " << pop.at(j).getCorrdinate().d << " " << pop.at(j).getCorrdinate().N << " " << pop.at(j).getCorrdinate().K << " " << pop.at(j).getCorrdinate().ub << " " << pop.at(j).getCorrdinate().uv << endl;
			}
			//cout << "2: " << pop.at(j).getPbest().d << " " << pop.at(j).getPbest().N << " " << pop.at(j).getPbest().K << " " << pop.at(j).getPbest().ub << " " << pop.at(j).getPbest().uv << endl;
			//更新pBest
			Particle tmp = pop.at(j).convertPbestToParticle();
			domination(tmp, pop.at(j), flag);
			if (flag == 1)
				pop.at(j).setPbest(pop.at(j).getCorrdinate());
			else if (flag == 0)
			{
				srand((unsigned int)time(NULL));
				if ((rand() / (double)RAND_MAX) < 0.5)
					pop.at(j).setPbest(pop.at(j).getCorrdinate());
			}
			//cout << "3: " << pop.at(j).getPbest().d << " " << pop.at(j).getPbest().N << " " << pop.at(j).getPbest().K << " " << pop.at(j).getPbest().ub << " " << pop.at(j).getPbest().uv << endl;
		}
		//cout << rep.size() << endl;

		//将新的非支配向量加入到Rep当中
		addTORep();

		//重新分配栅格
		createGrid(pop);

		//更新粒子的栅格的索引和坐标
		for (size_t k = 0; k < rep.size(); ++k)
		{
			findGridIndex(rep[k]);
			/*cout << "第i个粒子：" << k << endl;
			cout << "x,y坐标：" << rep[k].gridSubIndex.x << " " << rep[k].gridSubIndex.y << endl;
			cout << "粒子索引：" << rep[k].gridIndex << endl << endl;*/
		}

		//判断Rep是否已满,如果粒子数超过Rep的容量，需要删除相应密集栅格中的粒子
		if (rep.size() > REPNUM)
		{
			size_t num = rep.size() - REPNUM;
			for (size_t j = 0; j < num; ++j)
			{
				deleteOneRepMember();
			}
		}

		//更新w
		w = w*wdamp;

		//cout << "rep size:" << rep.size() << endl;
		//将得到的粒子写入文件
		char str[10] = { 0 };
		_itoa(i, str, 10);
		char filename[20] = "./repx_";
		strcat(filename, str);
		strcat(filename, ".txt");
		ofstream out(filename);
		for (int k = 0; k < rep.size() ;k++)
		{
			out << rep.at(k).f->F(rep.at(k).answer) << " " << rep.at(k).f->TNom(rep.at(k).answer) << endl;
		}
		out.close();
	}

	//for (size_t i = 0; i < rep.size(); ++i)
	//cout << i << endl << rep[i].getCorrdinate().d << " " << rep[i].getCorrdinate().N << " " << rep[i].getCorrdinate().K << " " << rep[i].getCorrdinate().ub << " " << rep[i].getCorrdinate().uv << " " << endl << endl;

}

