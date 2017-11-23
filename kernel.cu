#include <cuda.h>
#include <cuda_runtime.h>

#include <cstdio>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <conio.h>
#include <stddef.h>
#include <fstream>
#include <sstream>
#include <tchar.h>

#include <gmp_util.h>

typedef long big_int;
using mpfr::mpreal;

__shared__ int dimension, ciclo, x, cram;
mpfr::mpreal temp, aux, temp1, det00;

std::vector < std::vector <mpreal> > matrixChi;
std::vector < std::vector <mpreal> > matrixChi2;

void input(mpreal array[][3], mpreal array1[][1]);
mpreal determinent(mpreal array[][3]);
mpreal calculate(mpreal array[][3], int a, int b, int c);
mpreal copy(mpreal array[][3], mpreal array1[][1], int a);
void comp_copy(mpreal array[][3], mpreal array1[][3]);
mpreal determinent4(mpreal array[][4]);
std::vector<std::vector<mpreal> > newdet(mpreal array[][4], int col);
void input4(mpreal array[][4], mpreal array1[][1]);
mpreal copy4(mpreal array[][4], mpreal array1[][1], int a);
void comp_copy4(mpreal array[][4], mpreal array1[][4]);

__host__ __device__ void f() {
#ifdef __CUDA_ARCH__
	printf("Hebra en CUDA: %d\n", threadIdx.x);
#else
	printf("CUDA Funcionando!\n");
#endif
}

__global__ void kernel() {
	f();
}


int main()
{

	const int digits = 200;
	mpreal::set_default_prec(mpfr::digits2bits(digits));
	mpreal overflow = std::numeric_limits<mpreal>::max();



	ciclo = 1;

	kernel << <1, 1 >> >();
	if (cudaDeviceSynchronize() != cudaSuccess) {
		fprintf(stderr, "CUDA Fallo\n");
	}
	f();

	std::cout << "Corriendo variables con " << digits << " bits de largo\n\n";

	std::string s;
	x = 0;
	std::ifstream myReadFile;
	myReadFile.open("matrix.txt");
	if (myReadFile.is_open()) {
		myReadFile >> s;
		std::stringstream geek(s);
		dimension = 0;
		geek >> dimension;
		std::cout << "Dimensiones del problema: " << dimension << "\n";
	}
	myReadFile >> s;
	std::stringstream geek(s);
	cram = 0;
	geek >> cram;
	//generando CHI
	matrixChi.resize(dimension);
	for (int i = 0; i < dimension; i++)
	{
		matrixChi[i].resize(dimension + 1);
	}
	matrixChi2.resize(dimension - 1);
	for (int i = 0; i < dimension - 1; i++)
	{
		matrixChi2[i].resize(dimension);
	}
	//numeros desde file=matrix.txt

	if (myReadFile.is_open()) {
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension + 1; j++) {
				myReadFile >> s;
				std::stringstream geek(s);
				x = 0;
				geek >> x;
				matrixChi[i][j] = x;
			}
		}
	}
	myReadFile.close();
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension + 1; j++)
			if (j == dimension) std::cout << " | " << (mpreal)(matrixChi[i][j]) << " ";
			else std::cout << matrixChi[i][j] << " ";
			std::cout << "\n";
	}


	while (dimension > cram) {
		std::cout << "Generando Matrix de nivel " << dimension << "\n";


		//sacando pivote
		std::cout << "Pivote es " << matrixChi[0][0] << "\n";
		for (int i = 0; i < dimension - 1; i++)
			for (int j = 0; j < dimension; j++)
				matrixChi2[i][j] = (mpreal)(matrixChi[0][0] * matrixChi[i + 1][j + 1] - matrixChi[i + 1][0] * matrixChi[0][j + 1]);

		std::cout << "Condensacion a ciclo " << ciclo << "\n";
		ciclo++;
		for (int i = 0; i < dimension - 1; i++) {
			for (int j = 0; j < dimension; j++)
				if (j == dimension - 1) std::cout << " | " << (mpreal)(matrixChi2[i][j]) << " ";
				else std::cout << (mpreal)(matrixChi2[i][j]) << " ";
				std::cout << "\n";
		}

		//generando CHI
		dimension--;
		matrixChi.resize(dimension);
		for (int i = 0; i < dimension; i++)
		{
			matrixChi[i].resize(dimension + 1);
		}

		//poner datos de chi2 en chi
		for (int i = 0; i < dimension; i++)
			for (int j = 0; j < dimension + 1; j++)
				matrixChi[i][j] = (mpreal)matrixChi2[i][j];

		//generando la que se determina
		matrixChi2.resize(dimension - 1);
		for (int i = 0; i < dimension - 1; i++)
		{
			matrixChi2[i].resize(dimension);
		}
	}
	//comienza cramer

	std::cout << "\nDeterminantes de Cramer\n\n";

	bool sahi = true;
	while (sahi)
	{
		if (cram == 3 && cram != 4) {
			//long double permite valores de 1.8 × 10^308
			mpreal matrix[3][3];
			mpreal matrix1[3][1];
			mpreal reserve[3][3];
			mpreal detr[3];
			int sp1 = 0, teen = 1;
			int cont = 0;
			char in;
			//pasa los valores de matrixchi a la matriz de valores y resultados
			input(matrix, matrix1);
			//se respalda la matriz
			comp_copy(reserve, matrix);
			//se calcula el determinante general o coeficiente de la matriz
			det00 = determinent(matrix);
			while (sp1<3)
			{
				detr[cont] = copy(matrix, matrix1, sp1);
				comp_copy(matrix, reserve);
				cont++;
				sp1++;
			}
			cont = 0;
			while (cont<3)
			{
				std::cout << "x" << teen << " = " << (mpreal)detr[cont] << " /" << (mpreal)det00 << " [" << (mpreal)detr[cont] / (mpreal)det00 << "]" << std::endl;
				cont++;
				teen++;
			}
			mpreal x1 = ((mpreal)detr[0] / (mpreal)det00);
			mpreal x2 = ((mpreal)detr[1] / (mpreal)det00);
			mpreal x3 = ((mpreal)detr[2] / (mpreal)det00);
			mpreal err_total = 0;
			cont = 0;
			while (cont<3)
			{
				//calculo de errores
				mpreal resultado = matrix[cont][0] * x1 + matrix[cont][1] * x2 + matrix[cont][2] * x3;
				mpreal error = 100 - ((matrix1[cont][0] / resultado) * 100);
				std::cout << "[ESPERADO / OBTENIDO] para Ecuacion " << cont + 1 << "\n";
				std::cout << "[" << (mpreal)matrix1[cont][0] << " / " << (mpreal)resultado << "] - ";
				std::cout << "Error Relativo del " << abs((mpreal)error) << "%\n\n";
				err_total += abs((mpreal)error);
				cont++;
			}
			std::cout << "Error Promedio del " << (mpreal)err_total / 3 << "%\n";

			std::cout << "Finalizado, presione X para terminar\n" << overflow;
			std::cin >> in;
			if (in == 'x' || in == 'X')
				return 1;
		}
		else if (cram == 4) {
			mpreal matrix[4][4];
			mpreal matrix1[4][1];
			mpreal reserve[4][4];
			mpreal detr[4];
			int sp1 = 0, teen = 1;
			int cont = 0;
			char in;
			//pasa los valores de matrixchi a la matriz de valores y resultados
			input4(matrix, matrix1);
			//se respalda la matriz

			comp_copy4(reserve, matrix);
			//se calcula el determinante general o coeficiente de la matriz

			det00 = determinent4(matrix);

			while (sp1<4)
			{
				detr[cont] = copy4(matrix, matrix1, sp1);
				comp_copy4(matrix, reserve);
				cont++;
				sp1++;
			}
			cont = 0;
			while (cont<4)
			{
				std::cout << "x" << teen << " = " << (mpreal)detr[cont] << " /" << (mpreal)det00 << " [" << (mpreal)detr[cont] / (mpreal)det00 << "]" << std::endl;
				cont++;
				teen++;
			}
			mpreal x1 = ((mpreal)detr[0] / (mpreal)det00);
			mpreal x2 = ((mpreal)detr[1] / (mpreal)det00);
			mpreal x3 = ((mpreal)detr[2] / (mpreal)det00);
			mpreal x4 = ((mpreal)detr[3] / (mpreal)det00);
			mpreal err_total = 0;
			cont = 0;
			/*for (int k = 0; k<4; k++) {
				for (int l = 0; l<4; l++) {
					std::cout << matrix[k][l] << " ";
				}
				std::cout << "\n";
			}
			for (int l = 0; l<4; l++) {
				std::cout << matrix1[0][l] << "\n";
			}*/
			while (cont<4)
			{
				//calculo de errores
				/*std::cout << "test x1 " << matrix[cont][0] * x1 << "\n";
				std::cout << "test x2 " << matrix[cont][1] * x2 << "\n";
				std::cout << "test x3 " << matrix[cont][2] * x3 << "\n";
				std::cout << "test x4 " << matrix[cont][3] * x4 << "\n";
				std::cout << "test sum" << matrix[cont][0] * x1 + matrix[cont][1] * x2 + matrix[cont][2] * x3 + matrix[cont][3] * x4 << "\n";*/
				mpreal resultado = matrix[cont][0] * x1 + matrix[cont][1] * x2 + matrix[cont][2] * x3 + matrix[cont][3] * x4;
				/*std::cout << "Resultado " << resultado << "\n";*/
				mpreal error = 100 - ((matrix1[cont][0] / resultado) * 100);
				std::cout << "\n[ESPERADO / OBTENIDO] para Ecuacion " << cont + 1 << "\n";
				std::cout << "[" << (mpreal)matrix1[cont][0] << " / " << (mpreal)resultado << "] - ";
				std::cout << "Error Relativo del " << abs((mpreal)error) << "%\n\n";
				err_total += abs((mpreal)error);
				cont++;
			}
			std::cout << "Error Promedio del " << (float)err_total / 4 << "%\n";

			std::cout << "Finalizado, presione X para terminar\n";
			std::cin >> in;
			if (in == 'x' || in == 'X')
				return 1;
		}
	}
	std::cout.flush();
	return 0;
}
void input(mpreal array[][3], mpreal array1[][1])
{
	//traspasa los valores de la matriz X Y Z
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			array[i][j] = (mpreal)matrixChi[i][j];
		}
	}
	//traspasa los resultados de cada ecuacion X + Y + Z = RESULTADO
	for (int i = 0; i < 3; i++) {
		array1[i][0] = (mpreal)matrixChi[i][3];
	}
}
mpreal determinent(mpreal array[][3])
{
	int rows = 1, col = 1;
	int z = 0;
	temp = 0;
	int cont = 1;
	int x = 0;
	while (x<3)
	{
		//por cada columna los determinantes de 2x2 y los multiplica por l columna nueva dando el 3x3
		temp = temp + cont*(array[0][x] * calculate(array, rows, col, z));
		col = col * 0;
		z = z + cont;
		cont = cont*-1;
		x++;
	}
	if(cram == 3)
	std::cout << "\nDeterminante de la matrix 3x3 es " << temp << "\n\n";
	return temp;
}
mpreal calculate(mpreal array[][3], int a, int b, int c)
{
	//calcula los determinantes de 2x2
	temp1 = (array[a][b] * array[a + 1][b + 1 + c]) - (array[a + 1][b] * array[a][b + 1 + c]);
	return temp1;
}
mpreal copy(mpreal array[][3], mpreal array1[][1], int a)
{
	//traspasa los valores de la matriz de resultados a la matriz de valores
	int col = 0;
	temp = 0;
	while (col<3)
	{
		array[col][a] = array1[col][0];
		col++;
	}
	int i = 0, j = 0;
	while (i<3)
	{
		j = 0;
		while (j<3)
		{
			std::cout << array[i][j] << "  ";
			j++;
		}
		std::cout << std::endl;
		i++;
	}
	temp = determinent(array);
	return temp;
}
void comp_copy(mpreal array[][3], mpreal array1[][3])
{
	//traspasa los valores de la matriz de resultados a la matriz de valores
	int rows = 0, col = 0;
	while (rows<3)
	{
		col = 0;
		while (col<3)
		{
			array[rows][col] = array1[rows][col];
			col++;
		}
		rows++;
	}
}
mpreal determinent4(mpreal array[][4]) {
	int i, j, k;
	aux = 0;
	std::vector<std::vector<mpreal> > matrix(3, std::vector<mpreal>(3));
	mpreal matrixaux[3][3];
	for (i = 0; i<4; i++) {
		matrix = newdet(array, i);
		for (j = 0; j<3; j++) {
			for (k = 0; k<3; k++) {
				matrixaux[j][k] = matrix[j][k];
			}
		}
		aux = aux + pow(-1.0, (mpreal)i)*array[0][i] * (determinent(matrixaux));
	}
	std::cout << "\nDeterminante de la matrix 4x4 es " << aux << "\n\n";
	return aux;
}
std::vector<std::vector<mpreal> > newdet(mpreal array[][4], int col) {
	std::vector<std::vector<mpreal> > matrix(3, std::vector<mpreal>(3));
	int cont = 0, i, j;
	for (i = 1; i<4; i++) {
		for (j = 0; j<4; j++) {
			if (j != col) {
				matrix[i - 1][cont] = array[i][j];
				cont++;
			}
		}
		cont = 0;
	}
	return matrix;
}
void comp_copy4(mpreal array[][4], mpreal array1[][4])
{
	//traspasa los valores de la matriz de resultados a la matriz de valores
	int rows = 0, col = 0;
	while (rows<4)
	{
		col = 0;
		while (col<4)
		{
			array[rows][col] = array1[rows][col];
			col++;
		}
		rows++;
	}
}
mpreal copy4(mpreal array[][4], mpreal array1[][1], int a)
{
	//traspasa los valores de la matriz de resultados a la matriz de valores
	int col = 0;
	temp = 0;
	while (col<4)
	{
		array[col][a] = array1[col][0];
		col++;
	}
	int i = 0, j = 0;
	while (i<4)
	{
		j = 0;
		while (j<4)
		{
			std::cout << array[i][j] << "  ";
			j++;
		}
		std::cout << std::endl;
		i++;
	}
	std::cout << "\n";
	temp = determinent4(array);

	return temp;
}
void input4(mpreal array[][4], mpreal array1[][1])
{
	//traspasa los valores de la matriz X Y Z
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			array[i][j] = (mpreal)matrixChi[i][j];
		}
	}
	//traspasa los resultados de cada ecuacion X + Y + Z = RESULTADO
	for (int i = 0; i < 4; i++) {
		array1[i][0] = (mpreal)matrixChi[i][4];
	}
}