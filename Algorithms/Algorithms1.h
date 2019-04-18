#pragma once
#include <string>
class Algorithms1
{
private:
	//1-1
	const int pBase[10] = { 0, 1, 20, 300, 4000, 50000, 600000, 7000000, 80000000, 900000000 };
	const int iSize = 10;
	int pCount[10];
	int pValue[10];
	void pageCountCore(int* pValue,int n);
	void writeToFile(std::string s, int* p, int iSize);
	//1-2
	const int iLength = 6;
	char pChar[6];
public:
	Algorithms1();
	//һ�����ҳ��Ӱ�Ȼ��1��ʼ˳�����ֱ����Ȼ����i�����ҳ�밴��ͨ����ϰ�߱��ţ�
	//ÿ��ҳ�붼���������ǰ������0�����磬��6ҳ������6��ʾ��������06��006�ȡ���
	//�ּ�������Ҫ��Ը��������ҳ�룬i����������ȫ��ҳ���зֱ��õ����ٴ�����0, 1,
	//2������9��
	//nΪ���֣�sΪ���
	void GetPageCount(int m, std::string s);
	int Permutations(const int m, const int n);
	//���������ַ���
	int getIndex(std::string s);
	~Algorithms1();
};

