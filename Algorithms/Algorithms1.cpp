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
	primes();
	cnt = vecPrime.size();
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
//计算全排列
int Algorithms1::Permutations(const int m, const int n) {
	int iValue = 1;
	for (int i = m; i >= m-n+1; --i)
		iValue *= i;
	for (int i = 1; i <= n; ++i)
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
		int iFrom = 0;
		if (i == 0)
			iFrom = 1;
		else
			iFrom = s[i - 1] - 'a' + 2;//注意从前面的一个字符加1开始算起，到当前字符减一为止
		for (int j = iFrom; j < iPos; ++j) {
			iIndex += Permutations(26 - j, iSize - i - 1);//选择基数有从当前后面的的位数
		}		
	}
	++iIndex;//当前的字符
	return iIndex;
}

//很耗时间的质数求解
int Algorithms1::div(int x) {
	int iCount = 0, iQrt = (int)sqrt(x);
	if (iQrt*iQrt == x)
		++iCount;
	for (int i = 1; i < iQrt; ++i) {
		if (x%i == 0)
			iCount += 2;
	}
	return iCount;
}

void Algorithms1::primes() {
	//bool bPrime[65536];
	//for (int i = 0; i < 65536; ++i)
	//	bPrime[i] = true;
	//for (int j = 2; j < 65536; ++j) {
	//	int iStep = j;
	//	if (bPrime[iStep] == true) {
	//		while (iStep + j < 65536) {
	//			bPrime[iStep + j] = false;
	//			iStep += iStep + j;
	//		}
	//	}
	//}
	//for (int ii = 2; ii < 65536; ++ii) {
	//	if (bPrime[ii] == true)
	//		vecPrime.push_back(ii);
	//}
	bool get[MAXP + 1];
	long i;
	for (int i = 2; i <= MAXP; i++)
		get[i] = true;
	for(i=2;i<MAXP;i++)
		if (get[i]) {
			long j = i + i;
			while (j <= MAXP) {
				get[j] = false;
				j += i;
			}
		}
	long ii, j;
	for (ii = 2, j = 0; ii <= MAXP; ii++)
		if (get[ii])prim[++j] = ii;
	PCOUNT = j;
}
int Algorithms1::div2(int x) {
	int iSum = 1;
	int iSize = vecPrime.size();
	for (int i = 0; i < 65536 &&i< iSize&& vecPrime[i] <= sqrt(x); ++i) {
		int iCount = 0;
		while (x%vecPrime[i] == 0) {
			x = x / vecPrime[i];
			++iCount;
		}
		if (iCount!=0)
			iSum = iSum * (iCount + 1);
	}
	//额外超过质数65536，剪枝，不用求后面的质数情况
	if (x != 1)
		iSum *= 2;
	return iSum;
}
//int Max = 0;
//long long qpow(long long x, long long n) {
//	long long result = 1;
//	while (n) {
//		if (n & 1) {
//			result *= x;
//		}
//		x *= x;
//		n /= 2;
//	}
//	return result;
//}
//void Algorithms1::dfs(int point, int cnt1, long long now, long long num,int a,int b) {
//	if (now > b) {
//		return;
//	}
//	long long k = num * (cnt1 + 1);
//	if (now >= a) {
//		if (k > Max) Max = k;
//	}
//	for (int i = point; i < cnt && vecPrime[i] <= b; ++i) {
//		long long tmp = now * vecPrime[i];
//		if (tmp > b) return;
//		if (i == point) {
//			dfs(i, cnt1 + 1, tmp, num,a,b);
//		}
//		else {
//			long long k = num * (cnt1 + 1);
//			long long rest = b / now;
//			long long q = log(rest) / log(vecPrime[point]);
//			long long r2q = qpow(2, q);
//			if (r2q*k <= Max) return;
//			if ((a - 1) / now == b / now) return;
//			else if (now<a&&b / now>vecPrime[cnt - 1]) {
//				if (2 * k > Max) Max = 2 * k;
//			}
//			dfs(i, 1, tmp, k,a,b);
//		}
//	}
//}


//from是当前循环的质因数的位置，tot是num包含的质因数的个数，num是当前循环的值，low是最低值，up是最高值，并且原区间为[num*low,num*up]
void Algorithms1::search(long from, long tot, long num, long low, long up) {
	if(num>1)
		if ((tot > max) || (tot == max) && (num < numb)) {
			max = tot;
			numb = num;
		}
	//这句话要不要效果都一样
	if ((low == up) && (low > num))search(from, tot * 2, num*low, 1, 1);

	for (long i = from; i <= PCOUNT; i++) {
		//当前乘积prim[i]超过up，剪掉
		if (prim[i] > up)return;
		else {
			long j = prim[i], x = low - 1, y = up, n = num,t = tot, m = 1;
			while (true) {
				//如果当前乘积为prim[i] < a，而(low - 1) / prim[i] = up / prim[i], 则[a, b]不存在k的整数倍，直接剪掉。
				m++;
				t += tot;
				x /= j;
				y /= j;
				if (x == y)break;
				n *= j;
				search(i+1,t,n,x+1,y);
			}
			//当前的质因数每加1次，次数要乘以2，当前最大的质数的个数*当前质因数的个数tot都小于当前的max，就不用再循环到下一步。
			//因为后续的prime[i]出现的个数绝对比这个小。
			//up/numb=prim[i]^m,剩下因子数最多为m,个数为2^m
			m = 1 << m;
			if (tot < max / m) return;
		}
	}
}

int Algorithms1::getMaxDiv(const int iFrom, const int iTo) {

	//int iMax = div2(iFrom);
	//for (int i = iFrom + 1; i <= iTo; ++i) {
	//	cout << i << endl;
	//	int iTemp = div2(i);
	//	if (iMax < iTemp)
	//		iMax = iTemp;
	//}
	//return iMax;
	//Max = 0; 
	//dfs(0, 0, 1, 1,iFrom,iTo);
	if ((iFrom == 1) && (iTo == 1)) {
		max = 1;
		numb = 1;
	}
	else
	{
		max = 2;
		numb = 1;
		search(1,1,1,iFrom,iTo);
	}

	return max;
}

//求出最大和最小的下标
template<class T>
int mini(int n, T*x) {
	T temp = x[1];
	int iMin = 1;
	for (int i = 1; i <= n; ++i) {
		if (temp > x[i]) {
			temp = x[i];
			iMin = i;
		}
	}
	return iMin;
}

template<class T>
int maxi(int n, T* x) {
	T temp = x[1];
	int iMax = 1;
	for (int i = 1; i <= n; ++i) {
		if (temp < x[i]) {
			temp = x[i];
			iMax = i;
		}
	}
	return iMax;
}

double Algorithms1::maxgap(int n, double *x) {
	double minx = x[mini(n, x)], maxx = x[maxi(n, x)];
	//用n一2个等间距点分割区间[minx, maxx].
	//产生n一1个桶，每个桶i中用high〔i]和low「i}
	//分别存储分配给桶1的数中的最大数和最小数
	int *count = new int[iMaxSize];
	double *low = new double[iMaxSize], *high = new double[iMaxSize];
	//桶初始化
	for (int i = 1; i <= n - 1; i++) {
		count[i] = 0;
		low[i] = maxx;
		high[i] = minx;
	}
	//将n个数置于n-1个桶里面
	for (int i = 1; i <= n; i++) {
		int bucket = int((n - 1)*(x[i] - minx) / (maxx - minx)) + 1;
		count[bucket]++;
		if (x[i] < low[bucket])low[bucket] = x[i];
		if (x[i] > high[bucket])high[bucket] = x[i];
	}
	//此时，除了maxx和minx外的n一2个数被置于n一1个桶中
	//由鸽舍原理即知，至少有一个桶是空的.
	//这意味着最大间隙不会出现在同一桶中的两个数之间。
	//对每一个桶做一次线性扫描即可找出最大间隙，
	double tmp = 0, left = high[1];
	for (int i = 2; i <= n - 1; i++) {
		if (count[i]) {
			double thisgap = low[i] - left;
			if (thisgap > tmp)tmp = thisgap;
			left = high[i];
		}
	}
	delete count;
	delete low;
	delete high;
	return tmp;
}


//换列
void Algorithms1::tran_col(int p[iMaxSize1][iMaxSize1], int i,int j,int iRow) {
	for (int k = 0; k < iRow; k++) {
		/*如果是自己和自己交换时，就会变成0！！！
		p[k][j] = p[k][j] + p[k][i];
		p[k][i] = p[k][j] - p[k][i];
		p[k][j] = p[k][j] - p[k][i];*/
		int temp = p[k][j];
		p[k][j] = p[k][i];
		p[k][i] = temp;
	}
	if(i!=j)
		iCount++;
}
//反转一行
void Algorithms1::tran_row(int p[iMaxSize1][iMaxSize1], const int line, int iCol) {
	for (int i = 0; i < iCol; i++)
		p[line][i] = p[line][i] ^ 1;
	iCount++;
}
//赋值,p2的值赋予p1
void Algorithms1::copy(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1],int iRow,int iCol) {
	for (int i = 0; i < iRow; ++i) {
		for (int j = 0; j < iCol; ++j)
			p1[i][j] = p2[i][j];
	}
}
//判断两个数组是否相同
bool Algorithms1::isSame(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol) {
	for (int i = 0; i < iRow; ++i) {
		for (int j = 0; j < iCol; ++j)
			if(p1[i][j] != p2[i][j])
				return false;
	}
	return true;
}
//判断两个数组的某一列是否相同
bool Algorithms1::isSameCol(int p1[iMaxSize1][iMaxSize1], int i, int p2[iMaxSize1][iMaxSize1], int j, int iRow) {
	for (int k = 0; k < iRow; ++k) {
		if (p1[k][i] != p2[k][j])
			return false;
	}
	return true;
}
int Algorithms1::transf(int p1[iMaxSize1][iMaxSize1],int p2[iMaxSize1][iMaxSize1],int iRow,int iCol){
	int iMin = 0;
	int iBest = iCol + iRow + 1;
	int p3[iMaxSize1][iMaxSize1] = {};
	//保持现场
	copy(p3,p1,iRow,iCol);
	bool bSameCol = false;
	for (int i = 0; i < iCol; ++i) {
		copy(p1, p3, iRow, iCol);
		iCount = 0;
		tran_col(p1,0,i,iRow);
		//始终以第一列为准
		for (int j = 0; j < iRow; ++j) {
			if (p1[j][0] != p2[j][0])
				tran_row(p1, j, iCol);
		}
		
		//必须要从0开始，防止只有一列
		for (int j = 0; j < iCol; ++j) {
			bSameCol = false;
			for (int k = j; k < iCol; ++k) {
				//如果同一列相同就不用动，不同才互换。因为类似001111与101110，如果0与最后的1想换才一次，但是和前面几个要好多次反复
				if (isSameCol(p1, k, p2, j, iRow)) {
					if (k != j && isSameCol(p1, k, p2, k, iRow))
						continue;
					bSameCol = true;
					tran_col(p1, j, k, iRow);
					break;
				}
			}
			if (!bSameCol)
				break;
		}
		if (bSameCol&&iMin > iCount)
			iMin = iCount;
	}
	if (iMin == iRow + iCol + 1)
		return -1;
	else
		return iMin;
}