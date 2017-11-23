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

typedef long big_int;

int dimension, ciclo;
long double temp, temp1, det00;

std::vector < std::vector <long double> > matrixChi;
std::vector < std::vector <long double> > matrixChi2;


void input(long double array[][3], long double array1[][1]);
long double determinent(long double array[][3]);
long double calculate(long double array[][3], int a, int b, int c);
long double copy(long double array[][3], long double array1[][1], int a);
void comp_copy(long double array[][3], long double array1[][3]);

int main()
{
	ciclo = 1;

	std::cout << "Ingrese dimensiones del problema,\n";
	std::cin >> dimension;
	//generando CHI
	matrixChi.resize(dimension);
	for (int i = 0; i < dimension; i++)
	{
		matrixChi[i].resize(dimension+1);
	}
	matrixChi2.resize(dimension - 1);
	for (int i = 0; i < dimension - 1; i++)
	{
		matrixChi2[i].resize(dimension);
	}
	//numeros random
	for (int i = 0; i < dimension; i++)
		for (int j = 0; j < dimension+1; j++)
			matrixChi[i][j] = rand() % 19 + (-9);

	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension+1; j++)
			if (j == dimension) std::cout << " | " << (long double)(matrixChi[i][j]) << " ";
			else std::cout << matrixChi[i][j] << " ";
		std::cout << "\n";
	}


	while (dimension > 3) {
		std::cout << "Generando Matrix de nivel " << dimension << "\n";	


		//sacando pivote
		std::cout << "Pivote es " << matrixChi[0][0] << "\n";
		for (int i = 0; i < dimension - 1; i++)
			for (int j = 0; j < dimension; j++)
				matrixChi2[i][j] = (long double)(matrixChi[0][0] * matrixChi[i + 1][j + 1] - matrixChi[i + 1][0] * matrixChi[0][j + 1]);

		std::cout << "Condensacion a ciclo " << ciclo << "\n";
		ciclo++;
		for (int i = 0; i < dimension-1; i++) {
			for (int j = 0; j < dimension; j++)
				if(j == dimension -1 ) std::cout <<" | " << (long double)(matrixChi2[i][j]) << " ";
				else std::cout << (long double)(matrixChi2[i][j]) <<" ";
			std::cout << "\n";
		}

		//generando CHI
		dimension--;
		matrixChi.resize(dimension);
		for (int i = 0; i < dimension; i++)
		{
			matrixChi[i].resize(dimension+1);
		}

		//poner datos de chi2 en chi
		for (int i = 0; i < dimension; i++)
			for (int j = 0; j < dimension+1; j++)
				matrixChi[i][j] = (long double)matrixChi2[i][j];

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
		//long double permite valores de 1.8 × 10^308
		long double matrix[3][3];
		long double matrix1[3][1];
		long double reserve[3][3];
		long double detr[3], sp1 = 0, teen = 1;
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
			std::cout << "x" << teen << " = " << (long double)detr[cont] << " /" << (long double)det00 << " ["<< (long double)detr[cont]/ (long double)det00 << "]" << std::endl;
			cont++;
			teen++;
		}
		long double x1 = ((long double)detr[0] / (long double)det00);
		long double x2 = ((long double)detr[1] / (long double)det00);
		long double x3 = ((long double)detr[2] / (long double)det00);
		long double err_total=0;
		cont = 0;
		while (cont<3)
		{
			//calculo de errores
			long double resultado = matrix[cont][0] * x1 + matrix[cont][1] * x2 + matrix[cont][2] * x3;
			long double error = 100-((matrix1[cont][0] / resultado) * 100);
			std::cout << "[ESPERADO / OBTENIDO] para Ecuacion " << cont+1 << "\n";
			std::cout << "[" << (long double)matrix1[cont][0] << " / " << (long double)resultado << "] - ";
			std::cout << "Error Relativo del "<< abs((long double)error) <<"%\n\n";
			err_total += abs((long double)error);
			cont++;
		}
		std::cout << "Error Promedio del " << (float)err_total/3 << "%\n";

		std::cout << "Finalizado, presione X para terminar\n";
		std::cin >> in;
		if (in == 'x' || in == 'X')
			return 1;
	}
	std::cout.flush();
	return 0;
}
void input(long double array[][3], long double array1[][1])
{
	//traspasa los valores de la matriz X Y Z
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			array[i][j] = matrixChi[i][j];
		}
	}
	//traspasa los resultados de cada ecuacion X + Y + Z = RESULTADO
	for (int i = 0; i < 3; i++) {
			array1[i][0] = matrixChi[i][3];
	}
}
long double determinent(long double array[][3])
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
	std::cout << "\nDeterminante de la matrix es " << temp << "\n\n";
	return temp;
}
long double calculate(long double array[][3], int a, int b, int c)
{
	//calcula los determinantes de 2x2
	temp1 = (array[a][b] * array[a + 1][b + 1 + c]) - (array[a + 1][b] * array[a][b + 1 + c]);
	return temp1;
}
long double copy(long double array[][3], long double array1[][1], int a)
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
void comp_copy(long double array[][3], long double array1[][3])
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
