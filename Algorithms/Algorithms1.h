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
	const int iMaxSize = 200015;
	static const int iMaxSize1 = 150;
	//����
	int iCount = 0;
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

	//����϶���⣺���� n ��ʵ��X_1,X_2,X_3,X_4,X_5,X_6,X_n 
	//���� n ������ʵ�������� 2 ����֮�������ֵ��������κ�ʵ������ȡ��������ʱO(1) ����ƽ�����϶���������ʱ���㷨��
	double maxgap(int n, double *x);
	//��mXn(m<=100,n<=100)��������������ų�һ��m��n�еĽ�����С�ÿһö�� �һ����泯�ϻ��泯�ϡ������ֱ�ʾ���״̬��
	//0 ��ʾ������泯�ϣ�1 ��ʾ���泯�ϡ� ���������Ϸ�Ĺ����ǣ� ��1��ÿ�οɽ���һ�н�ҷ���������ԭ����λ���ϣ� .
	//��2��ÿ�ο���ѡ 2 �У������� 2 �н�ҵ�λ�á� 
	int transf(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol);
	//����
	void tran_col(int p[iMaxSize1][iMaxSize1], int i, int j, int iRow);
	//��תһ��
	void tran_row(int p[iMaxSize1][iMaxSize1], const int line, int iCol);
	//��ֵ,p2��ֵ����p1
	void copy(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol);
	//�ж����������Ƿ���ͬ
	bool isSame(int p1[iMaxSize1][iMaxSize1], int p2[iMaxSize1][iMaxSize1], int iRow, int iCol);
	//�ж����������ĳһ���Ƿ���ͬ
	bool isSameCol(int p1[iMaxSize1][iMaxSize1], int i, int p2[iMaxSize1][iMaxSize1], int j, int iRow);
	~Algorithms1();
};

