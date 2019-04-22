// Algorithms.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//代码运行文件

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Algorithms1.h"
using namespace std;



int main()
{
	string s = "D:\\workspace\\VS2017\\solutions\\ch1\\prog13\\test\\";
	string s1 = "D:\\workspace\\VS2017\\solutions\\ch1\\prog13\\temp\\";
	Algorithms1 ag;
	for (char i = 0; i < 11; ++i) {
		string sOut = s1 +"countOut" +to_string(i);
		string sIn = s + "divi"+to_string(i) +"\.in" ;
		ifstream fInput(sIn);
		int iValue = 0;
		//char pS[100];
		//fInput.getline(pS, 100);
		ofstream oFile(sOut);
		//while (!fInput.eof()) {
			//string s1;
			//fInput.getline(pS,100);
			//iValue=atoi(pS);
			//s1 = pS;
			//ag.GetPageCount(iValue,sOut);
			//pageCount(iValue);
			/*ofstream oFile(sOut);
			for (int i = 0; i < iSize; i++) {
				oFile << pCount[i] << endl;
			}
			oFile.close();*/
			//int iIndex=ag.getIndex(s1);
			int iFrom, iTo;
			fInput >> iFrom >> iTo;
			//1000000 2000000
			int iResult=ag.getMaxDiv(iFrom, iTo);
			//int iResult = ag.getMaxDiv(121, 121);
			oFile << iResult << endl;
		//}
		//ag.getIndex("ahou");
		oFile.close();
		fInput.close();
	}

    std::cout << "Hello World!\n"; 
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
