#include <iostream>
#include "CZLPVector.h"
#include "suhlupdate.h"

int main()
{
	int ColNum = 7, RowNum = 5;
	suhlupdate suhl(RowNum, ColNum);
	vector<int> LStart({ 0, 4, 7, 9, 10, 10 });
	vector<int> LIndex({ 1, 2, 3, 4, 2, 3, 4, 3, 4, 4 });
	vector<double> LValue({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
	vector<int> LRStart({ 0, 0, 1, 3, 6, 10 });
	vector<int> LRIndex({ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 });
	vector<double> LRValue({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
	vector<int> LPivot({ 0, 1, 2, 3, 4 });
	vector<int> LPivotLookup({ 0, 1, 2, 3, 4 });
	suhl.copy_L(LValue, LIndex, LStart, LRValue, LRIndex, LRStart, LPivot, LPivotLookup);
	vector<int> URStart({ 0, 4, 7, 9, 10 });
	vector<int> UREnd({ 4, 7, 9, 10, 10 });
	vector<int> URIndex({ 1, 2, 3, 4, 2, 3, 4, 3, 4, 4 });
	vector<double> URValue({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
	vector<int> UStart({ 0, 0, 1, 3, 6 });
	vector<int> UEnd({ 0, 1, 3, 6, 10 });
	vector<int> UIndex({ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 });
	vector<double> UValue({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
	vector<int> UPivot({ 0, 1, 2, 3, 4 });
	vector<int> UPivotLookup({ 0, 1, 2, 3, 4 });
	vector<double> UPivotValue({ 1, 1, 1, 1, 1 });
	vector<int> URspace(5, 0);
	suhl.copy_U(UPivotValue, UValue, UIndex, UStart, UEnd, URValue, URIndex, URStart, UREnd, URspace, UPivotLookup, UPivot);

	vector<double> res = { 0, 1, 0, 1, 0 };
	double *aq1 = new double[RowNum];
	for (int i = 0; i < res.size(); i++)
	{
		aq1[i] = res[i];
	}
	CZLPVector aq(5);
	aq.set(aq1, 5);
	//suhl.ftranl(aq);
	int position, insertPosition;
	// position = suhl.search_last_nonzero_position(*aq);
	// insertPosition = suhl.search_column(position);
	suhl.update(aq, 3);
	// cout << insertPosition;
	suhl.show();
	delete[] aq1;
	return 0;
}
