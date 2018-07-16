#include "NavSolJps2.h"

CONSOLE_APP_MAIN
{
	Vector <InputData> inDatVect;
	Vector <Ephemeris> eph;
	
	//Зададим массив входящих данных
	double inData[2][7];
	inData[0][0] = 20;
	inData[0][1] = 4;
	inData[0][2] = 21;
	inData[0][3] = 5;
	inData[0][4] = 22;
/*	
	inData[1][0] = 20317.7206562;
	inData[1][1] = 22309.8467601915;
	inData[1][2] = 19305.7438067749;
	inData[1][3] = 19344.8207480251;
	inData[1][4] = 23093.3796233308;
*/
	
	//Попробуем использовать модельные измерения
	double xx = -3782627.48060055 / 1000;
	double yy = 706759.897227085 / 1000;
	double zz = 6330779.72431116 / 1000;

	//Эфемериды аппаратов
	double ephem[5][3];
	ephem[0][0] = 8852.877527; ephem[0][1] = 19812.499427; ephem[0][2] = 13455.272608;	//20-й аппарат
	ephem[1][0] = -12198.743176; ephem[1][1] = 6387.390091; ephem[1][2] = 21491.421479;	//4-й аппарат
	ephem[2][0] = 10976.098464; ephem[2][1] = 2307.666364; ephem[2][2] = 22904.680281;	//21
	ephem[3][0] = 5644.592226; ephem[3][1] = 12957.107325; ephem[3][2] = 21231.360978;	//5
	ephem[4][0] = 5752.717420; ephem[4][1] = -16816.611187; ephem[4][2] = 18235.134882;	//22
	
//	inData[1][0] = pow((pow((xx - ephem[0][0]), 2) + pow((yy - ephem[0][1]), 2) + pow((zz - ephem[0][2]), 2)), 0.5);
//	inData[1][1] = pow((pow((xx - ephem[1][0]), 2) + pow((yy - ephem[1][1]), 2) + pow((zz - ephem[1][2]), 2)), 0.5);
//	inData[1][2] = pow((pow((xx - ephem[2][0]), 2) + pow((yy - ephem[2][1]), 2) + pow((zz - ephem[2][2]), 2)), 0.5);
//	inData[1][3] = pow((pow((xx - ephem[3][0]), 2) + pow((yy - ephem[3][1]), 2) + pow((zz - ephem[3][2]), 2)), 0.5);
//	inData[1][4] = pow((pow((xx - ephem[4][0]), 2) + pow((yy - ephem[4][1]), 2) + pow((zz - ephem[4][2]), 2)), 0.5);

	inData[1][0] = 20317.7206562;
	inData[1][1] = 22309.8467601915;
	inData[1][2] = 19305.7438067749;
	inData[1][3] = 19344.8207480251;
	inData[1][4] = 23093.3796233308;
	
	//Заполним массивы входящих данных
	for (int i = 0; i < 5; i++){
		inDatVect.Add(InputData(inData[0][i], inData[1][i], 0.0));
		eph.Add(Ephemeris(0, ephem[i][0], ephem[i][1], ephem[i][2]));
	}
	
	LOG("Data test output");
	for (int i = 0; i < inDatVect.GetCount(); i++){
		RDUMP(inDatVect[i].psRangLen);
	}
	
	LOG("In the procedure");
	//Найдём координаты станции
	//Coord statPos = NavProblemPsRang(inDatVect, eph);
	Coord statPos = NewNavProb(inDatVect, eph);
}

//Новый алгоритм
Coord NewNavProb(Vector <InputData> &inData, Vector <Ephemeris> &eph){
	Coord outDat;
	
	//Зададим начальные данные
	double x[4] = {0, 0, 0, 0};
	Vector <double> roApr;
	Vector <hLine> H;
	
	//Вычислим вектор приближенных дальностей
	roApr.Clear();
	for (int i = 0; i < eph.GetCount(); i++){
		roApr.Add(pow((pow(eph[i].x - x[0], 2) + pow(eph[i].y - x[1], 2) + pow(eph[i].z - x[2], 2)), 0.5));
	}
	
	//Вычислим матрицу топоцентрических направлений
	H.Clear();
	for (int i = 0; i < roApr.GetCount(); i++){
		double hX = (eph[i].x - x[0]) / roApr[i];
		double hY = (eph[i].y - x[1]) / roApr[i];
		double hZ = (eph[i].z - x[2]) / roApr[i];
		H.Add(hLine(hX, hY, hZ));
	}
	
	//Вычислим произведение транспонированной матрицы и исходной
	double HtH[4][4] = {{0,0,0,0},
						{0,0,0,0},
						{0,0,0,0},
						{0,0,0,0}};
	HtHPowCalc(H, HtH);
	//Матрица всегда бует квадратной и 4х4
	
	//Приведём матрицу в нужный вид
	double HtHMatr[100][100];
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			HtHMatr[i][j]  = HtH[i][j];
		}
	}
	
	//Находим обратную к HtH матрицу
	double revHtH[100][100];
	inverse(HtHMatr, revHtH, 4, det(HtHMatr,4));

	LOG("Original matr test output");
	for (int i = 0; i < 4; i++){
		LOG(AsString(HtHMatr[i][0]) << " " << AsString(HtHMatr[i][1]) << " " << AsString(HtHMatr[i][2]) << " " << AsString(HtHMatr[i][3]));
	}
	
	LOG("\nReverce matr test output");
	for (int i = 0; i < 4; i++){
		LOG(AsString(revHtH[i][0]) << " " << AsString(revHtH[i][1]) << " " << AsString(revHtH[i][2]) << " " << AsString(revHtH[i][3]));
	}
	
	//Умножим обратную матрицу на транспонированную матрицу входящих данных
	double mH[100][100];
	for (int i = 0; i < H.GetCount(); i++){
		mH[0][i] = H[i].xx;
		mH[1][i] = H[i].yy;
		mH[2][i] = H[i].zz;
		mH[3][i] = H[i].one;
	}
	MatrMultiply(revHtH, mH);
	
	//Вычислим результирующий вектор координат
	//Для этого преобразуем Вектор измерений в матрицу
	
	double obsVect[100][100];
	for (int i = 0; i < 100; i++){
		for (int j = 0; j < 100; j++){
			if (j == 0 && i < inData.GetCount()){
				obsVect[i][j] = inData[i].psRangLen;
			}else{
				obsVect[i][j] = 0;
			}
		}
	}
	
	MatrMultiply(mH, obsVect);
	
	//Теперь умножим матрицу на mH вектор измерений obsVect
	//MatrMultiply(mH, obsVect);

	//Выведем результат
	LOG("coord test output");
	String outStr = "";
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 2; j++){
			outStr += AsString(obsVect[i][j]);
			outStr += " ";
		}
		outStr += "\n";
	}
	LOG(outStr);
	
	RDUMP(pow(pow(obsVect[0][0], 2) + pow(obsVect[1][0], 2) + pow(obsVect[2][0], 2) , 0.5));
	
	return outDat;
}

void HtHPowCalc(Vector <hLine> &H, double HtH[4][4]){
	for (int i = 0; i < H.GetCount(); i++){
		HtH[0][0] += H[i].xx * H[i].xx;
		HtH[0][1] += H[i].xx * H[i].yy;
		HtH[0][2] += H[i].xx * H[i].zz;
		HtH[0][3] += H[i].xx * H[i].one;
		HtH[1][1] += H[i].yy * H[i].yy;
		HtH[2][2] += H[i].zz * H[i].zz;
		HtH[3][3] += H[i].one * H[i].one;
		HtH[1][2] += H[i].yy * H[i].zz;
		HtH[1][3] += H[i].yy * H[i].one;
		HtH[2][3] += H[i].zz * H[i].one;
	}
	HtH[1][0] = HtH[0][1];
	HtH[2][0] = HtH[0][2];
	HtH[2][1] = HtH[1][2];
	HtH[3][0] = HtH[0][3];
	HtH[3][1] = HtH[1][3];
	HtH[3][2] = HtH[2][3];
}

void SqrMatrVectMultiply(float m[100][100], float v[100], int order){
	float resVect[100];
	

	for (int i = 0; i < order; i++){
		resVect[i] = 0;
		for (int j = 0; j < order; j++){
			resVect[i] += m[j][i] * v[i];
		}
	}
	
	LOG("result vect is");
	for (int i = 0; i < order; i++){
		v[i] = resVect[i];
		RDUMP(v[i]);
	}
}

void MatrMultiply(double a[100][100], double b[100][100]){
	double res[100][100];
	
	for (int i = 0; i < 100; i++){
		for (int j = 0; j < 100; j++){
			res[i][j] = 0;
		}
	}
	
	for (int i = 0; i < 100; i++){
		for (int j = 0; j < 100; j++){
			for (int k = 0; k < 100; k++){
				res[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	
	for (int i = 0; i < 100; i++){
		for (int j = 0; j < 100; j++){
			b[i][j] = res[i][j];
		}
	}
}