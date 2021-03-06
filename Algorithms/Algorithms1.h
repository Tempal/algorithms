#pragma once
#include <string>
#include<vector>
class Algorithms1
{
private:
	//1-1
	const int pBase[10] = { 0, 1, 20, 300, 4000, 50000, 600000, 7000000, 80000000, 900000000 };
	const int iSize = 10;
	int pCount[10];
	int pValue[10];
	int cnt;

	void pageCountCore(int* pValue,int n);
	void writeToFile(std::string s, int* p, int iSize);
	//1-2
	const int iLength = 6;
	char pChar[6];
	//记录了1-1000000000所有的质因数
	std::vector<int> vecPrime;
	//很耗时间的质数求解
	int div(int x);
	void primes();
	//正常的求解
	int div2(int x);
	static const long int MAXP = 100000;
	long prim[MAXP];
	long int max, numb, PCOUNT;//max存放最多约束个数，numb存放约数个数最多的数
	void search(long from, long tot, long num, long low, long up);
	const int iMaxSize = 200015;
	static const int iMaxSize1 = 150;
	//次数
	int iCount = 0;
public:
	Algorithms1();
	//一本书的页码从白然数1开始顺序编码直到白然数，i。书的页码按照通常的习惯编排，
	//每个页码都不含多余的前导数字0。例如，第6页用数字6表示，而不是06或006等。数
	//字计数问题要求对给定书的总页码，i，计算出书的全部页码中分别用到多少次数字0, 1,
	//2，…，9。
	//n为数字，s为输出
	void GetPageCount(int m, std::string s);
	int Permutations(const int m, const int n);
	//升序排序字符串1-2
	//在数据加密和数据压缩中常需要对特殊的字符串进行编码。给定的字母表 A 由 26 个小写英文字母组成 A={a,b,…,z}。
	//该字母表产生的升序字符串是指字符串中字母按照从左到右出现的次序与字母在字母表中出现的次序相同，
	//且每个字符最多出现 1 次。
	int getIndex(std::string s);


	//正整数 x 的约数是能整除 x 的正整数。正整数 x 的约数个数记为 div(x)。例如，1，2，
	//5，10 都是正整数 10 的约数，且 div(10) = 4。设 a 和 b 是 2 个正整数，a≤b，找出 a 和 b
	//	之间约数个数最多的数 x
	int getMaxDiv(const int iFrom, const int iTo);
	//void dfs(int point, int cnt1, long long now, long long num, int a, int b);

	//最大间隙问题：给定 n 个实数X_1,X_2,X_3,X_4,X_5,X_6,X_n 
	//求这 n 个数在实轴上相邻 2 个数之间的最大差值。假设对任何实数的下取整函数耗时O(1) ，设计解最大间隙问题的线性时间算法。
	double maxgap(int n, double *x);
	//有mXn(m<=100,n<=100)个金币在桌面上排成一个m行n列的金币阵列。每一枚金 币或正面朝上或背面朝上。用数字表示金币状态，
	//0 表示金币正面朝上，1 表示背面朝上。 金币阵列游戏的规则是： （1）每次可将任一行金币翻过来放在原来的位置上； .
	//（2）每次可任选 2 列，交换这 2 列金币的位置。 
	int transf(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol);
	//换列
	void tran_col(int p[iMaxSize1][iMaxSize1], int i, int j, int iRow);
	//反转一行
	void tran_row(int p[iMaxSize1][iMaxSize1], const int line, int iCol);
	//赋值,p2的值赋予p1
	void copy(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol);
	//判断两个数组是否相同
	bool isSame(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol);
	//判断两个数组的某一列是否相同
	bool isSameCol(int p1[iMaxSize1][iMaxSize1], int i, int p2[iMaxSize1][iMaxSize1], int j, int iRow);
	~Algorithms1();
};

