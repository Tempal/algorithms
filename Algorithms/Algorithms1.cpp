#include "pch.h"
#include "Algorithms1.h"

#include <stdio.h>
#include <cmath>
#include<iostream>
#include <fstream> 
#include <string>
#include<algorithm>
#include<map>
#include<stack>
#include<vector>
using namespace std;


Algorithms1::Algorithms1()
{
	for (int i = 0; i < iSize; ++i){
		pCount[i] = 0;
		pValue[i] = 0;
	}
}


Algorithms1::~Algorithms1()
{
}
void Algorithms1::writeToFile(string s,int* p,int iSize) {
	ofstream oFile(s);
	for (int i = 0; i < iSize; i++) {
		oFile << *p << endl;
	}
	oFile.close();
}
//int getValueFromP(int* pValue, int n) {
//	int iSum = 0;
//	for(int i=0;i<n;++i)
//		iSum += pValue[i]*pow(10,i);
//	return iSum;
//}
//n为位数
void Algorithms1::pageCountCore(int* pValue, int n) {
	if (n < 1 || n>10) return;
	//原先的做法，考虑过渡，主要是没有必要分解，并且函数getValueFromP不需要独立
	//if (n == 1) {
	//	for (int i = 0; i <= pValue[n - 1]; ++i)
	//		pCount[i] += 1;
	//	return;
	//}
	//if (pValue[n - 1] == 0) {
	//	int iStep = 0;
	//	while (n > 1 && pValue[n - 1] == 0) {
	//		--n;
	//		++iStep;
	//	}			
	//	if (n == 1 && pValue[n - 1] == 0) {
	//		pCount[0] += iStep + 1;
	//		return;
	//	}			
	//	else {
	//		pCount[0] = iStep * getValueFromP(pValue, n - iStep);
	//		pageCountCore(pValue, n - iStep);
	//	}
	//}
	//else {
	//	for (int i = 0; i < iSize; ++i) {
	//		//从0算起，所以最后面的乘法不能为pValue[n - 1] -1
	//		pCount[i] += pBase[n - 1] * (pValue[n - 1] );			
	//	}
	//	for (int i = 0; i < pValue[n - 1]; ++i)
	//		pCount[i] += pow(10, n - 1);
	//	//从0算起，所以要加1
	//	pCount[pValue[n - 1]] += getValueFromP(pValue, n - 1)+1;
	//	pageCountCore(pValue, n - 1);
	//}

	if (n < 1 || n>10)
		return;
	//最开始的一块，算出所有的NXXXX出现的次数。
	for (int i = 0; i < pValue[n - 1]; ++i) {
		pCount[i] += pow(10, (n - 1));
	}
	//稳定块，算出所有稳定块出现的次数
	for (int i = 0; i < iSize; ++i) {
		pCount[i] += pBase[n - 1] * pValue[n - 1];
	}
	for (int i = 0; i < n - 1; ++i) {
		pCount[pValue[n - 1]] += pValue[i] * pow(10, i);
	}
	++pCount[pValue[n - 1]];//算上0
	pageCountCore(pValue, n - 1);
}
void Algorithms1::GetPageCount(int m, std::string s) {
	for (int i = 0; i < iSize; ++i) {
		pCount[i] = 0;
		pValue[i] = 0;
	}
	int n = log10(m)+1;
	int iCount = 0,iBase=10;
	while (m / iBase) {
		pValue[iCount]=m%iBase;
		m=m/iBase;
		++iCount;
	}
	pValue[iCount] = m;
	pageCountCore(pValue,n);
	int iReduce = 0;
	for (int i = 0; i < n; i++)
		iReduce += pow(10,i);
	pCount[0] -= iReduce;
	writeToFile(s,pCount,iSize);
}
//计算Cm
int Algorithms1::Permutations(const int m, const int n) {
	int iValue = 1;
	for (int i = m; i <= n; ++i)
		iValue *= i;
	for (int i = 1; i <= m - n + 1; ++i)
		iValue /= i;
	return iValue;
}

int Algorithms1::getIndex(string s) {
	const int pBase[] = {0,26,351,2951,17901,83681};
	int iIndex = pBase[s.length()-1],iSize=s.length(),iBase=26;
	for (int i = 0; i < iSize; ++i) {
		if (i + 1 < iSize&&s[i] >= s[i + 1])
			return 0;
		int iPos = s[i] - 'a' +1;
		for (int j = 1; j < iPos; ++j) {
			iIndex += Permutations(26 - j, iSize - i - 1);
		}		
	}
	++iIndex;
	return iIndex;
}