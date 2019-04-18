#pragma once
#include <string>
class Algorithms1
{
private:
	const int pBase[10] = { 0, 1, 20, 300, 4000, 50000, 600000, 7000000, 80000000, 900000000 };
	const int iSize = 10;
	long long int pCount[10];
	int pValue[10];
	void pageCountCore(int* pValue,int n);
	void writeToFile(std::string s);
public:
	Algorithms1();
	//一本书的页码从白然数1开始顺序编码直到白然数，i。书的页码按照通常的习惯编排，
	//每个页码都不含多余的前导数字0。例如，第6页用数字6表示，而不是06或006等。数
	//字计数问题要求对给定书的总页码，i，计算出书的全部页码中分别用到多少次数字0, 1,
	//2，…，9。
	//n为数字，s为输出
	void GetPageCount(int m, std::string s);
	~Algorithms1();
};

