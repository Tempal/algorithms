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
	//��¼��1-1000000000���е�������
	std::vector<int> vecPrime;
	//�ܺ�ʱ����������
	int div(int x);
	void primes();
	//���������
	int div2(int x);
	static const long int MAXP = 100000;
	long prim[MAXP];
	long int max, numb, PCOUNT;//max������Լ��������numb���Լ������������
	void search(long from, long tot, long num, long low, long up);

public:
	Algorithms1();
	//һ�����ҳ��Ӱ�Ȼ��1��ʼ˳�����ֱ����Ȼ����i�����ҳ�밴��ͨ����ϰ�߱��ţ�
	//ÿ��ҳ�붼���������ǰ������0�����磬��6ҳ������6��ʾ��������06��006�ȡ���
	//�ּ�������Ҫ��Ը��������ҳ�룬i����������ȫ��ҳ���зֱ��õ����ٴ�����0, 1,
	//2������9��
	//nΪ���֣�sΪ���
	void GetPageCount(int m, std::string s);
	int Permutations(const int m, const int n);
	//���������ַ���1-2
	//�����ݼ��ܺ�����ѹ���г���Ҫ��������ַ������б��롣��������ĸ�� A �� 26 ��СдӢ����ĸ��� A={a,b,��,z}��
	//����ĸ������������ַ�����ָ�ַ�������ĸ���մ����ҳ��ֵĴ�������ĸ����ĸ���г��ֵĴ�����ͬ��
	//��ÿ���ַ������� 1 �Ρ�
	int getIndex(std::string s);


	//������ x ��Լ���������� x ���������������� x ��Լ��������Ϊ div(x)�����磬1��2��
	//5��10 ���������� 10 ��Լ������ div(10) = 4���� a �� b �� 2 ����������a��b���ҳ� a �� b
	//	֮��Լ������������ x
	int getMaxDiv(const int iFrom, const int iTo);
	//void dfs(int point, int cnt1, long long now, long long num, int a, int b);
	~Algorithms1();
};

