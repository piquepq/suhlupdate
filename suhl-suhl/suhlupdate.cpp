#include "suhlupdate.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
using std::copy;
using std::vector;

#define TINY 1e-5

suhlupdate::suhlupdate(int RowNum, int ColNum) : RowNum(RowNum), ColNum(ColNum)
{
	Steps = 0;
	workarray.assign(RowNum, 0);
}

void suhlupdate::copy_perm(vector<int> perm, vector<int> permback)
{
	this->perm = perm;
	this->permback = permback;
}
void suhlupdate::copy_L(vector<double> Lvalues, vector<int> Lindex, vector<int> Lstart, vector<double> LRvalues, vector<int> LRindex, vector<int> LRstart, vector<int> Lpivotindex, vector<int> Lpivotlookup)
{
	this->Lvalues = Lvalues;
	this->Lindex = Lindex;
	this->Lstart = Lstart;
	this->LRvalues = LRvalues;
	this->LRindex = LRindex;
	this->LRstart = LRstart;
	this->Lpivotindex = Lpivotindex;
	this->Lpivotlookup = Lpivotlookup;
}
void suhlupdate::copy_U(vector<double> Upivotvalues, vector<double> Uvalues, vector<int> Uindex, vector<int> Ustart, vector<int> Uend, vector<double> URvalues, vector<int> URindex, vector<int> URstart, vector<int> URend, vector<int> URspace, vector<int> Upivotlookup,
						vector<int> Upivotindex)
{
	this->Upivotvalues = Upivotvalues;
	this->Uvalues = Uvalues;
	this->Uindex = Uindex;
	this->Ustart = Ustart;
	this->Uend = Uend;
	this->URvalues = URvalues;
	this->URindex = URindex;
	this->URstart = URstart;
	this->URend = URend;
	this->URspace = URspace;
	this->Upivotlookup = Upivotlookup;
	this->Upivotindex = Upivotindex;
}
void suhlupdate::btran(CZLPVector &aq)
{
	btranu(aq);
	btranl(aq);
}
void suhlupdate::ftran(CZLPVector &aq)
{
	ftranl(aq);
	ftranu(aq);
}
void suhlupdate::btranl(CZLPVector &rhs)
{

	int RHScount = 0;
	int *RHSindex = &rhs.Index[0];
	double *RHSarray = &rhs.Array[0];

	const int *LRstart = &this->LRstart[0];
	const int *LRindex = this->LRindex.size() > 0 ? &this->LRindex[0] : NULL;
	const double *LRvalue = this->LRvalues.size() > 0 ? &this->LRvalues[0] : NULL;

	for (int i = RowNum - 1; i >= 0; i--)
	{
		int pivotRow = Lpivotindex[i];
		const double pivotX = RHSarray[pivotRow];
		if (fabs(pivotX) > 0)
		{
			RHSindex[RHScount++] = pivotRow;
			RHSarray[pivotRow] = pivotX;
			const int start = LRstart[i];
			const int end = LRstart[i + 1];
			for (int k = start; k < end; k++)
				RHSarray[LRindex[k]] -= pivotX * LRvalue[k];
		}
		else
			RHSarray[pivotRow] = 0;
	}

	rhs.set_nonzeronum(RHScount);
	rhs.set_ispacked(true);
	rhs.pack();
}
void suhlupdate::btranu(CZLPVector &rhs)
{

	int RHScount = 0;
	int *RHSindex = &rhs.Index[0];
	double *RHSarray = &rhs.Array[0];

	const int *URstart = &this->URstart[0];
	const int *URend = &this->URend[0];
	const int *URindex = &this->URindex[0];
	const double *URvalue = &this->URvalues[0];

	int UpivotCount = Upivotindex.size();
	for (int iLogic = 0; iLogic < UpivotCount; iLogic++)
	{
		if (Upivotindex[iLogic] == -1)
			continue;

		const int pivotRow = Upivotindex[iLogic];
		double pivotX = RHSarray[pivotRow];
		if (fabs(pivotX) > 0)
		{
			pivotX /= Upivotvalues[iLogic];
			RHSindex[RHScount++] = pivotRow;
			RHSarray[pivotRow] = pivotX;
			const int start = URstart[iLogic];
			const int end = URend[iLogic];
			for (int k = start; k < end; k++)
				RHSarray[URindex[k]] -= pivotX * URvalue[k];
		}
		else
			RHSarray[pivotRow] = 0;
	}

	rhs.set_nonzeronum(RHScount);
	rhs.set_ispacked(true);
	rhs.pack();

	btranr(rhs);
}
void suhlupdate::ftranl(CZLPVector &rhs)
{

	int RHScount = 0;
	int *RHSindex = &rhs.Index[0];
	double *RHSarray = &rhs.Array[0];

	const int *Lstart = &this->Lstart[0];
	const int *Lindex = this->Lindex.size() > 0 ? &this->Lindex[0] : NULL;
	const double *Lvalue = this->Lvalues.size() > 0 ? &this->Lvalues[0] : NULL;
	// Transform
	for (int i = 0; i < RowNum; i++)
	{
		int pivotRow = Lpivotindex[i];
		const double pivotX = RHSarray[pivotRow];
		if (fabs(pivotX) > 0)
		{
			RHSindex[RHScount++] = pivotRow; //？？？？index都是变换前的
			const int start = Lstart[i];
			const int end = Lstart[i + 1];
			for (int k = start; k < end; k++)
				RHSarray[Lindex[k]] -= pivotX * Lvalue[k]; //
		}
		else
			RHSarray[pivotRow] = 0;
	}

	rhs.set_nonzeronum(RHScount);
	rhs.set_ispacked(true);
	rhs.pack();
}
void suhlupdate::ftranu(CZLPVector &rhs)
{
	ftranr(rhs);

	int RHScount = 0;
	int *RHSindex = &rhs.Index[0];
	double *RHSarray = &rhs.Array[0];

	const int *Ustart = &this->Ustart[0];
	const int *Uend = &this->Uend[0];
	const int *Uindex = this->Uindex.size() > 0 ? &this->Uindex[0] : NULL;
	const double *Uvalue = this->Uvalues.size() > 0 ? &this->Uvalues[0] : NULL;

	int UpivotCount = Upivotindex.size();
	for (int iLogic = UpivotCount - 1; iLogic >= 0; iLogic--)
	{

		if (Upivotindex[iLogic] == -1)
			continue;

		const int pivotRow = Upivotindex[iLogic];
		double pivotX = RHSarray[pivotRow];
		if (fabs(pivotX) > 0)
		{
			pivotX /= Upivotvalues[iLogic];
			RHSindex[RHScount++] = pivotRow;
			RHSarray[pivotRow] = pivotX;
			const int start = Ustart[iLogic];
			const int end = Uend[iLogic];
			for (int k = start; k < end; k++)
				RHSarray[Uindex[k]] -= pivotX * Uvalue[k];
		}
		else
			RHSarray[pivotRow] = 0;
	}
	rhs.set_nonzeronum(RHScount);
	rhs.set_ispacked(true);
	rhs.pack();
}
void suhlupdate::ftranr(CZLPVector &vector)
{
	const int PFpivotCount = RpivotIndex.size();
	int *RpivotIndex = NULL;
	if (this->RpivotIndex.size() > 0) //看是不是第一次更新
		RpivotIndex = (int *)&this->RpivotIndex[0];

	const int *PFstart = this->Rstart.size() > 0 ? &this->Rstart[0] : NULL;
	const int *PFindex = this->Rindex.size() > 0 ? &this->Rindex[0] : NULL;
	const double *PFvalue = this->Rvalues.size() > 0 ? &this->Rvalues[0] : NULL;

	int RHScount = vector.get_nonzeronum();
	int *RHSindex = &vector.Index[0];
	double *RHSarray = &vector.Array[0];

	for (int i = 0; i < PFpivotCount; i++)
	{
		int iRow = RpivotIndex[i];
		double value0 = RHSarray[iRow];
		double value1 = value0;
		const int start = PFstart[i];
		const int end = PFstart[i + 1];
		for (int k = start; k < end; k++)
			value1 -= RHSarray[PFindex[k]] * PFvalue[k];

		if (value0 || value1)
		{
			if (value0 == 0)
				RHSindex[RHScount++] = iRow;
			RHSarray[iRow] =
				(fabs(value1) < TINY) ? 0 : value1;
		}
	}

	vector.set_nonzeronum(RHScount);
	vector.set_ispacked(true);
	vector.pack();
}
void suhlupdate::btranr(CZLPVector &vector)
{

	const int PFpivotCount = RpivotIndex.size();
	const int *PFpivotIndex =
		this->RpivotIndex.size() > 0 ? &this->RpivotIndex[0] : NULL;
	const int *PFstart = this->Rstart.size() > 0 ? &this->Rstart[0] : NULL;
	const int *PFindex = this->Rindex.size() > 0 ? &this->Rindex[0] : NULL;
	const double *PFvalue = this->Rvalues.size() > 0 ? &this->Rvalues[0] : NULL;

	int RHScount = vector.get_nonzeronum();
	int *RHSindex = &vector.Index[0];
	double *RHSarray = &vector.Array[0];

	for (int i = PFpivotCount - 1; i >= 0; i--)
	{
		int pivotRow = PFpivotIndex[i];
		double pivotX = RHSarray[pivotRow];
		if (pivotX)
		{
			const int start = PFstart[i];
			const int end = PFstart[i + 1];
			for (int k = start; k < end; k++)
			{
				int iRow = PFindex[k];
				double value0 = RHSarray[iRow];
				double value1 = value0 - pivotX * PFvalue[k];
				if (value0 == 0)
					RHSindex[RHScount++] = iRow;
				RHSarray[iRow] =
					(fabs(value1) < TINY) ? 0 : value1;
			}
		}
	}

	vector.set_nonzeronum(RHScount);
	vector.set_ispacked(true);
	vector.pack();
}

int suhlupdate::search_last_nonzero_position(CZLPVector &aq)
{

	int i = aq.Array.size() - 1;
	while (aq.Array[i] == 0)
	{
		i--;
	}
	return i;
}


//还要加排除出基出基列
int suhlupdate::search_column(int aqLen)
{	
	if (aqLen == RowNum - 1)
		return Ustart.size() ;
	int res = 0;
	for (int i = 0; i < Upivotindex.size(); i++)
	{
		
		if (Upivotindex[i] != -1) 
		{
			res++;
		}
		if (res == aqLen) 
		{
			break;
		}
	}
	//因为要插在该列后面一列
	res += 1;
	return res;
}

void suhlupdate::expandrow(int row)
{

	//将其初始化为0向量
	for (int i = 0; i < RowNum; i++)
	{
		workarray[i] = 0;
	}
	//将第p行填充进workarray
	assert(row < RowNum);
	for (int i = URstart[row]; i < URend[row]; i++)
	{
		int position = Upivotlookup[Uindex[i]];
		int res = 0;
		for (int i = 0; i < position; i++) {
			if (Upivotindex[i] = -1) {
				res++;
			}
		}
		position -= res;
		workarray[position] = URvalues[i];
	}
}
/*
Uindex[i]=j，指该元素变换前在矩阵的第j个位置

Upivotindex[i]=j,指该元素在变换后的上三角矩阵U第i列在原矩阵第j列

Upivotlookup[j]=i




*/
void suhlupdate::update(CZLPVector &aq, int ColOut)
{
	//ColOut为矩阵在B中第几列，查找变换后在第几列
	//标记删除出基列,查找变换后在Start第几个
	int out = Upivotlookup[ColOut];
	double pivot ;
	
	//求得新入基向量的最后一个元素在哪里
	int aqLen = aq.Index[aq.get_nonzeronum() - 1];

	// 找到排列矩阵，看具体需要该行插入到第几列
//???????
	int insertposition = search_column(aqLen);

	//标记删除，U中几列对应原矩阵第几列
	Upivotindex[out] = -1;
	double alpha = Upivotvalues[out];

	//将要删除的行展开，方便后面进行消元
	expandrow(out);

	// 搜索删除列存储中的对应元素
	for (int k = Ustart[out]; k < Uend[out]; k++)
	{

		//Uindex对应的是变换前的坐标
		int iLogic = Upivotlookup[Uindex[k]];
		int iFind = URstart[iLogic];
		int iLast = --URend[iLogic];
		for (; iFind <= iLast; iFind++)
			if (URindex[iFind] == ColOut)
				break;
		URindex[iFind] = URindex[iLast];
		URvalues[iFind] = URvalues[iLast]; 
		URspace[iLogic]++;
	}


	//删除out行的元素，先在列存储里删除，后面行元素会改变的元素再最后添加,可能会形成fill-in
	for (int k = URstart[out]; k < URend[out]; k++)
	{
		int iLogic = Upivotlookup[URindex[k]];

		int iFind = Ustart[iLogic];
		int iLast = --Uend[iLogic];
		for (; iFind <= iLast; iFind++)
			if (Uindex[iFind] == ColOut)

				break;
		assert(iFind < Uindex.size());
		assert(iLast < Uindex.size());
		Uindex[iFind] = Uindex[iLast];
		Uvalues[iFind] = Uvalues[iLast];
	}

	// 将aq插入到U中,直接insert到insertposition列之前，pivot值在求出消元部分后最后加入 
	Ustart.insert(Ustart.begin() + insertposition, Uindex.size());
	for (int i = 0; i < aq.Index.size(); i++)
	{
		if (aq.Index[i] != ColOut)
		{
			Uindex.push_back(aq.Index[i]);
			Uvalues.push_back(aq.PackedArray[i]);
			continue;
		}
		if (aq.Index[i] == ColOut) {
			pivot = aq.PackedArray[i];
			continue; 
		}
	}
	Uend.insert(Uend.begin() + insertposition, Uindex.size());
	int UstartX = Ustart[insertposition];
	int UendX = Uend[insertposition];

	// 将aq的元素插入到UR中
	for (int k = UstartX; k < UendX; k++)
	{
		int iLogic = Upivotlookup[Uindex[k]];

		// 如果空间不够将该行复制到最后
		if (URspace[iLogic] == 0)
		{

			int row_start = URstart[iLogic];
			int row_count = URend[iLogic] - row_start;
			int new_start = URindex.size();
			int new_space = row_count * 1.1 + 5;

			URindex.resize(new_start + new_space);
			URvalues.resize(new_start + new_space);

			int iFrom = row_start;
			int iEnd = row_start + row_count;
			int iTo = new_start;
			copy(&URindex[iFrom], &URindex[iEnd], &URindex[iTo]);
			copy(&URvalues[iFrom], &URvalues[iEnd], &URvalues[iTo]);

			URstart[iLogic] = new_start;
			URend[iLogic] = new_start + row_count;
			URspace[iLogic] = new_space - row_count;
		}
		URspace[iLogic]--;
		int iPut = URend[iLogic]++;
		URindex[iPut] = ColOut;
		URvalues[iPut] = Uvalues[k];
	}
	//建立新向量消元出基部分
	//cout << insertposition - out << endl;
	double *eliminationPart = new double[insertposition - out];
	memset(eliminationPart, 0, (insertposition - out) * sizeof(double));
	//用回代方法求消元部分
	for (int i = out + 1; i < insertposition; i++)
	{
		double pivotvalue = Upivotvalues[i];
		double res = workarray[i] / pivotvalue;
		eliminationPart[i - out] = res;
		for (int k = URstart[i]; k < URend[k]; k++)
		{
			workarray[Upivotlookup[URindex[k]]] -= URvalues[k];
		}
		res /= pivotvalue;
		pivot += res * aq.Array[i];
		eliminationPart[i - out] = res;
	}


	CZLPVector ep(insertposition - out);
	ep.set(eliminationPart, insertposition - out);
	delete[] eliminationPart;

	// 将第p行与新的储存UR的部分，有问题？？？
	URstart.insert(URstart.begin() + insertposition, URindex.size());
	for (int i = insertposition + 1; i < workarray.size(); i++)
	{
		if (workarray[i] != 0)
		{
			URindex.push_back(Upivotindex[i]);
			URvalues.push_back(workarray[i]);
		}
	}
	URend.insert(URend.begin() + insertposition, URindex.size());
	URspace.insert(URspace.begin() + insertposition, URspace[out] + URend[out] - URstart[out]);

	// 将新pivot更新
	Upivotlookup[ColOut] = insertposition;
	Upivotindex.insert(Upivotindex.begin() + insertposition, ColOut);
	Upivotvalues.insert(Upivotvalues.begin() + insertposition, pivot);

	// 将新的消元部分储存在eta矩阵中
	//将原来第p行的对应于新的pivot


	//对付R矩阵
	for (int i = 0; i < ep.Array.size(); i++)
	{
		//不用储存对角元，因为是1
		if (i != out)
		{
			Rvalues.push_back(ep.Array[i]);
			Rindex.push_back(i
			
			); //要变换前的value
		}
	}

	//看是否需要重新分解
	Steps++;
	if (Steps > 100)
	{
		flag = true;
	}
}
void suhlupdate::show()
{
	std::ofstream matrix;
	matrix.open("matrix.txt");
	int m = Ustart.size();
	int *colstart = &Ustart[0];
	int *rowindex = &Uindex[0];
	int *end = &Uend[0];
	double *value = &Uvalues[0];
	std::vector<double> dense(m * m, 0);
	matrix << m << endl;
	for (int j = 0; j < m; ++j)
	{
		for (int i = colstart[j]; i < end[j]; ++i)
		{
			int row, column;
			row = rowindex[i];
			column = j;
			dense[row * m + column] = value[i];
		}
	}
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			matrix << std::setprecision(10) << dense[j + i * m] << " ";
		}
		matrix << endl;
	}
	for (auto i : Upivotindex)
	{
		matrix << i << " ";
	}
	matrix << endl;
	for (auto i : Rindex)
	{
		matrix << i << " ";
	}
	matrix << endl;
	for (auto i : Rvalues)
	{
		matrix << i << " ";
	}
	matrix << endl;
	matrix.close();
}


void suhlupdate::create_two_pointer() {
	cptr.assign(Uindex.size()*1.2, -1);

	rptr.assign(URindex.size()*1.2, -1);

	//cptr与uindex以及uvalues的元素一一对应，每次改变时3者同时改变,值为-1的元素表示为空

	for(int i=0;i<RowNum;i++)
	{
	//一列一列的处理
		for (int j = Ustart[i]; j < Uend[i]; j++) {
			int row = Uindex[j];
			for (int k = URstart[row]; k < URend[row]; k++) {
				if (URindex[k] == i) {
					cptr[j] = k;
					break;
				}
			}
		}
	}
	for (int i = 0; i < RowNum; i++)
	{
		//一行一行的处理
		for (int j = URstart[i]; j < URend[i]; j++) {
			int col = URindex[j];
			for (int k = Ustart[col]; k < Uend[col]; k++) {
				if (Uindex[k] == i) {
					rptr[j] = k;
					break;
				}
			}
		}
	}
}
void update_two_pointer();