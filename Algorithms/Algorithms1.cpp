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
//nΪλ��
void Algorithms1::pageCountCore(int* pValue, int n) {
	if (n < 1 || n>10) return;
	//ԭ�ȵ����������ǹ��ɣ���Ҫ��û�б�Ҫ�ֽ⣬���Һ���getValueFromP����Ҫ����
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
	//		//��0�������������ĳ˷�����ΪpValue[n - 1] -1
	//		pCount[i] += pBase[n - 1] * (pValue[n - 1] );			
	//	}
	//	for (int i = 0; i < pValue[n - 1]; ++i)
	//		pCount[i] += pow(10, n - 1);
	//	//��0��������Ҫ��1
	//	pCount[pValue[n - 1]] += getValueFromP(pValue, n - 1)+1;
	//	pageCountCore(pValue, n - 1);
	//}

	if (n < 1 || n>10)
		return;
	//�ʼ��һ�飬������е�NXXXX���ֵĴ�����
	for (int i = 0; i < pValue[n - 1]; ++i) {
		pCount[i] += pow(10, (n - 1));
	}
	//�ȶ��飬��������ȶ�����ֵĴ���
	for (int i = 0; i < iSize; ++i) {
		pCount[i] += pBase[n - 1] * pValue[n - 1];
	}
	for (int i = 0; i < n - 1; ++i) {
		pCount[pValue[n - 1]] += pValue[i] * pow(10, i);
	}
	++pCount[pValue[n - 1]];//����0
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
//����ȫ����
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
			iFrom = s[i - 1] - 'a' + 2;//ע���ǰ���һ���ַ���1��ʼ���𣬵���ǰ�ַ���һΪֹ
		for (int j = iFrom; j < iPos; ++j) {
			iIndex += Permutations(26 - j, iSize - i - 1);//ѡ������дӵ�ǰ����ĵ�λ��
		}		
	}
	++iIndex;//��ǰ���ַ�
	return iIndex;
}

//�ܺ�ʱ����������
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
	//���ⳬ������65536����֦�������������������
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


//from�ǵ�ǰѭ������������λ�ã�tot��num�������������ĸ�����num�ǵ�ǰѭ����ֵ��low�����ֵ��up�����ֵ������ԭ����Ϊ[num*low,num*up]
void Algorithms1::search(long from, long tot, long num, long low, long up) {
	if(num>1)
		if ((tot > max) || (tot == max) && (num < numb)) {
			max = tot;
			numb = num;
		}
	//��仰Ҫ��ҪЧ����һ��
	if ((low == up) && (low > num))search(from, tot * 2, num*low, 1, 1);

	for (long i = from; i <= PCOUNT; i++) {
		//��ǰ�˻�prim[i]����up������
		if (prim[i] > up)return;
		else {
			long j = prim[i], x = low - 1, y = up, n = num,t = tot, m = 1;
			while (true) {
				//�����ǰ�˻�Ϊprim[i] < a����(low - 1) / prim[i] = up / prim[i], ��[a, b]������k����������ֱ�Ӽ�����
				m++;
				t += tot;
				x /= j;
				y /= j;
				if (x == y)break;
				n *= j;
				search(i+1,t,n,x+1,y);
			}
			//��ǰ��������ÿ��1�Σ�����Ҫ����2����ǰ���������ĸ���*��ǰ�������ĸ���tot��С�ڵ�ǰ��max���Ͳ�����ѭ������һ����
			//��Ϊ������prime[i]���ֵĸ������Ա����С��
			//up/numb=prim[i]^m,ʣ�����������Ϊm,����Ϊ2^m
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

//���������С���±�
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
	//��nһ2���ȼ���ָ�����[minx, maxx].
	//����nһ1��Ͱ��ÿ��Ͱi����high��i]��low��i}
	//�ֱ�洢�����Ͱ1�����е����������С��
	int *count = new int[iMaxSize];
	double *low = new double[iMaxSize], *high = new double[iMaxSize];
	//Ͱ��ʼ��
	for (int i = 1; i <= n - 1; i++) {
		count[i] = 0;
		low[i] = maxx;
		high[i] = minx;
	}
	//��n��������n-1��Ͱ����
	for (int i = 1; i <= n; i++) {
		int bucket = int((n - 1)*(x[i] - minx) / (maxx - minx)) + 1;
		count[bucket]++;
		if (x[i] < low[bucket])low[bucket] = x[i];
		if (x[i] > high[bucket])high[bucket] = x[i];
	}
	//��ʱ������maxx��minx���nһ2����������nһ1��Ͱ��
	//�ɸ���ԭ��֪��������һ��Ͱ�ǿյ�.
	//����ζ������϶���������ͬһͰ�е�������֮�䡣
	//��ÿһ��Ͱ��һ������ɨ�輴���ҳ�����϶��
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


//����
void Algorithms1::tran_col(int p[iMaxSize1][iMaxSize1], int i,int j,int iRow) {
	for (int k = 0; k < iRow; k++) {
		/*������Լ����Լ�����ʱ���ͻ���0������
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
//��תһ��
void Algorithms1::tran_row(int p[iMaxSize1][iMaxSize1], const int line, int iCol) {
	for (int i = 0; i < iCol; i++)
		p[line][i] = p[line][i] ^ 1;
	iCount++;
}
//��ֵ,p2��ֵ����p1
void Algorithms1::copy(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1],int iRow,int iCol) {
	for (int i = 0; i < iRow; ++i) {
		for (int j = 0; j < iCol; ++j)
			p1[i][j] = p2[i][j];
	}
}
//�ж����������Ƿ���ͬ
bool Algorithms1::isSame(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol) {
	for (int i = 0; i < iRow; ++i) {
		for (int j = 0; j < iCol; ++j)
			if(p1[i][j] != p2[i][j])
				return false;
	}
	return true;
}
//�ж����������ĳһ���Ƿ���ͬ
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
	//�����ֳ�
	copy(p3,p1,iRow,iCol);
	bool bSameCol = false;
	for (int i = 0; i < iCol; ++i) {
		copy(p1, p3, iRow, iCol);
		iCount = 0;
		tran_col(p1,0,i,iRow);
		//ʼ���Ե�һ��Ϊ׼
		for (int j = 0; j < iRow; ++j) {
			if (p1[j][0] != p2[j][0])
				tran_row(p1, j, iCol);
		}
		
		//����Ҫ��0��ʼ����ֹֻ��һ��
		for (int j = 0; j < iCol; ++j) {
			bSameCol = false;
			for (int k = j; k < iCol; ++k) {
				//���ͬһ����ͬ�Ͳ��ö�����ͬ�Ż�������Ϊ����001111��101110�����0������1�뻻��һ�Σ����Ǻ�ǰ�漸��Ҫ�ö�η���
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