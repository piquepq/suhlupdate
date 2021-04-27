#pragma once
#include <vector>
#include <iostream>
#include "CZLPVector.h"
using std::cin;
using std::cout;
using std::endl;
using std::vector;
class suhlupdate
{
public:
	suhlupdate(int RowNum, int ColNum);

	void copy_perm(vector<int> perm, vector<int> permback);
	void copy_L(vector<double> Lvalues, vector<int> Lindex, vector<int> Lstart, vector<double> LRvalues, vector<int> LRindex, vector<int> LRstart, vector<int> Lpivotindex, vector<int> Lpivotlookup);
	void copy_U(vector<double> Upivotvalues, vector<double> Uvalues, vector<int> Uindex, vector<int> Ustart,
				vector<int> Uend, vector<double> URvalues, vector<int> URindex, vector<int> URstart, vector<int> URend, vector<int> URspace, vector<int> Upivotlookup,
				vector<int> Upivotindex);
	void btran(CZLPVector &aq);
	void ftran(CZLPVector &aq);
	void btranl(CZLPVector &aq);
	void btranu(CZLPVector &aq);
	void ftranl(CZLPVector &aq);
	void ftranu(CZLPVector &aq);
	void ftranr(CZLPVector &aq);
	void btranr(CZLPVector &aq);

	int search_last_nonzero_position(CZLPVector &aq);

	void create_two_pointer();
	void update_two_pointer();
	//Ѱ����Ҫ���뵽��һ��
	int search_column(int position);

	void update(CZLPVector &aq, int ColOut);

	void expandrow(int row);
	void show();

private:
	//����Ĵ�С
	const int RowNum;
	const int ColNum;
	//Ŀǰ�����˶��ٲ�
	int Steps;
	//�Ƿ���Ҫ���¸���
	bool flag;
	//���о���P
	vector<int> perm;
	vector<int> permback;
	//L����

	vector<int> Lpivotindex;
	vector<int> Lpivotlookup;
	vector<double> Lvalues;
	vector<int> Lindex;
	vector<int> Lstart;

	vector<double> LRvalues;
	vector<int> LRindex;
	vector<int> LRstart;

	//U����
	vector<double> Upivotvalues;
	vector<int> Upivotindex;
	vector<int> Upivotlookup;
	vector<double> Uvalues;
	vector<int> Uindex;
	vector<int> Ustart;
	vector<int> Uend;

	vector<double> URvalues;
	vector<int> URindex;
	vector<int> URstart;
	vector<int> URend;
	vector<int> URspace;

	vector<int> cptr;
	vector<int> rptr;

	//R����,ÿһ��Ϊһ��eta����
	vector<double> Rpivotvalues;
	vector<int> RpivotIndex;
	vector<double> Rvalues;
	vector<int> Rindex;
	vector<int> Rstart;

	vector<double> workarray;
};
