// Algorithms.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//代码运行文件

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Algorithms1.h"
using namespace std;



const int iMaxSize = 200015;
int main()
{
	string s = "D:\\workspace\\VS2017\\solutions\\ch1\\prog15\\test\\";
	string s1 = "D:\\workspace\\VS2017\\solutions\\ch1\\prog15\\temp\\";
	Algorithms1 ag;
	for (char i = 0; i < 7; ++i) {
		string sOut = s1 +"gapOut" +to_string(i);
		string sIn = s + "gap"+to_string(i) +"\.in" ;
		ifstream fInput(sIn);
		ofstream oFile(sOut);
		int iValue = 0;
		char pS[100];
		fInput.getline(pS, 100);
		int iTime = atoi(pS);
		double *dIn=new double[iMaxSize];
		memset(dIn, 0, iMaxSize*sizeof(double));
		for (int i = 1; i <= iTime; ++i)
			fInput >> dIn[i];
		double dMax= ag.maxgap(iTime, dIn);
		delete dIn;
		oFile << dMax;



		//while (iTime > 0) {
		//	int p1[iMaxSize][iMaxSize] = {}, p2[iMaxSize][iMaxSize] = {};
		//	int iRow, iCol;
		//	fInput >> iRow >> iCol;
		//	int iMin = iRow + iCol+1;
		//	for (int i = 0; i < iRow; ++i) {
		//		for (int j = 0; j < iCol; ++j) {
		//			fInput>>p1[i][j];
		//		}
		//	}
		//	for (int i = 0; i < iRow; ++i) {
		//		for (int j = 0; j < iCol; ++j) {
		//			fInput >> p2[i][j];
		//		}
		//	}
		//	int iBest = iCol + iRow + 1;
		//	int p3[iMaxSize][iMaxSize] = {};
		//	//保持现场
		//	copy(p3,p1,iRow,iCol);
		//	bool bSameCol = false;
		//	for (int i = 0; i < iCol; ++i) {
		//		copy(p1, p3, iRow, iCol);
		//		iCount = 0;
		//		tran_col(p1,0,i,iRow);
		//		//始终以第一列为准
		//		for (int j = 0; j < iRow; ++j) {
		//			if (p1[j][0] != p2[j][0])
		//				tran_row(p1, j, iCol);
		//		}
		//		
		//		//必须要从0开始，防止只有一列
		//		for (int j = 0; j < iCol; ++j) {
		//			bSameCol = false;
		//			for (int k = j; k < iCol; ++k) {
		//				//如果同一列相同就不用动，不同才互换。因为类似001111与101110，如果0与最后的1想换才一次，但是和前面几个要好多次反复
		//				if (isSameCol(p1, k, p2, j, iRow)) {
		//					if (k != j && isSameCol(p1, k, p2, k, iRow))
		//						continue;
		//					bSameCol = true;
		//					tran_col(p1, j, k, iRow);
		//					break;
		//				}
		//			}
		//			if (!bSameCol)
		//				break;
		//		}
		//		if (bSameCol&&iMin > iCount)
		//			iMin = iCount;
		//	}
		//	if(iMin== iRow + iCol + 1)
		//		oFile << -1 << endl;
		//	else
		//		oFile << iMin << endl;
		//	--iTime;
		//}
		


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
			//int iFrom, iTo;
			//fInput >> iFrom >> iTo;
			//1000000 2000000
			//int iResult=ag.getMaxDiv(iFrom, iTo);
			//int iResult = ag.getMaxDiv(121, 121);
			//oFile << iResult << endl;
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
