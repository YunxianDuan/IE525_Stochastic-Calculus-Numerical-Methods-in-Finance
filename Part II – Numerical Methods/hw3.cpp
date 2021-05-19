//  hw3.cpp
//  Created by DYX on 2021/5/2.

#include <iostream>
#include <cmath>
#include "F:/MSFE/semester 1/IE 523 Financial computing/NEWMAT/newmat11/newmatap.h"
#include "F:/MSFE/semester 1/IE 523 Financial computing/NEWMAT/newmat11/newmat.h"
#include "F:/MSFE/semester 1/IE 523 Financial computing/NEWMAT/newmat11/newmatio.h"

double U = 120;
double L = 80;
double x_max = log(U);
double x_min = log(L);
double M = 799;
double delta_x = (x_max - x_min) / (M + 1);

double T = 0.25;
// double delta_t = 0.001;
// int N = T / delta_t;

double r = 0.01;
double q = 0.015;
double sigma = 0.2;
double mu = r - q - 0.5 * pow(sigma, 2);
double K = 100;

double a = pow(sigma, 2) / (2 * pow(delta_x, 2)) - mu / (2 * delta_x);
double b = r + pow(sigma, 2) / pow(delta_x, 2);
double c = pow(sigma, 2) / (2 * pow(delta_x, 2)) + mu / (2 * delta_x);


double max(float a, float b)
{
	return (b < a) ? a : b;
}

Matrix matrix_x()
{
	Matrix x(M, 1);
	x = 0.0;
	for (int i = 1; i <= M; i++)
	{
		x(i, 1) = x_min + i * delta_x;
	}
	return x;
}

Matrix matrix_f0()
{
	Matrix f0(M, 1);
	f0 = 0.0;
	Matrix x = matrix_x();
	for (int i = 1; i <= M; i++)
	{
		f0(i, 1) = max(exp(x(i, 1)) - K, 0);
	}
	return f0;
}

Matrix matrix_A()
{
	Matrix A(M, M);
	A = 0.0;
	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			if (i == (j - 1))
				A(i, j) = -c;
			if (i == j)
				A(i, j) = b;
			if (i == (j + 1))
				A(i, j) = -a;
		}
	}
	return A;
}

Matrix implicit_euler(int N)
{
	double delta_t = T / N;

	Matrix A = matrix_A();
	Matrix f0 = matrix_f0();
	Matrix fn(M, 1);
	IdentityMatrix I(M);
	Matrix Inverse = (I + delta_t * A);

	Matrix ai(M, 1), bi(M, 1), ci(M, 1), ci_p(M - 1, 1), di_p(M, 1);
	ai(1, 1) = 0, bi(1, 1) = Inverse(1, 1), ci(1, 1) = Inverse(1, 2);
	ai(M, 1) = Inverse(M, M - 1), bi(M, 1) = Inverse(M, M), ci(M, 1) = 0;
	for (int i = 2; i < M; i++)
	{
		ai(i, 1) = Inverse(i, i - 1);
		bi(i, 1) = Inverse(i, i);
		ci(i, 1) = Inverse(i, i + 1);
	}

	fn = f0;
	ci_p(1, 1) = ci(1, 1) / bi(1, 1);
	for (int i = 2; i <M; i++) 
	{
		ci_p(i, 1) = ci(i, 1) / (bi(i, 1) - ai(i, 1) * ci_p(i - 1, 1));
	}

	for (int i = 1; i <= N; i++)
	{
		di_p(1, 1) = fn(1, 1) / bi(1, 1);
		for (int j = 2; j <= M; j++)
		{
			di_p(j, 1) = (fn(j, 1) - ai(j, 1) * di_p(j - 1, 1)) / (bi(j, 1) - ai(j, 1) * ci_p(j - 1, 1));
		}

		fn(M, 1) = di_p(M, 1);
		for (int j = (M - 1); j >=1; j--)
		{
			fn(j, 1) = di_p(j, 1) - ci_p(j, 1) * fn(j + 1, 1);
		}
	}
	return fn;
}

Matrix crank_nicolson(int N)
{
	double delta_t = T / N;
	
	Matrix A = matrix_A();
	Matrix f0 = matrix_f0();
	Matrix fn(M, 1);
	IdentityMatrix I(M);
	Matrix Inverse_plus = (I + 0.5 * delta_t * A);
	Matrix Inverse_minus = (I - 0.5 * delta_t * A);
	// Matrix Inverse = Inverse_plus * Inverse_minus;

	Matrix ai(M, 1), bi(M, 1), ci(M, 1), ci_p(M - 1, 1), di_p(M, 1);
	ai(1, 1) = 0, bi(1, 1) = Inverse_plus(1, 1), ci(1, 1) = Inverse_plus(1, 2);
	ai(M, 1) = Inverse_plus(M, M - 1), bi(M, 1) = Inverse_plus(M, M), ci(M, 1) = 0;
	for (int i = 2; i < M; i++)
	{
		ai(i, 1) = Inverse_plus(i, i - 1);
		bi(i, 1) = Inverse_plus(i, i);
		ci(i, 1) = Inverse_plus(i, i + 1);
	}

	fn = f0;
	ci_p(1, 1) = ci(1, 1) / bi(1, 1);
	for (int i = 2; i < M; i++)
	{
		ci_p(i, 1) = ci(i, 1) / (bi(i, 1) - ai(i, 1) * ci_p(i - 1, 1));
	}

	for (int i = 1; i <= N; i++)
	{
		fn = Inverse_minus * fn;
		di_p(1, 1) = fn(1, 1) / bi(1, 1);
		for (int j = 2; j <= M; j++)
		{
			di_p(j, 1) = (fn(j, 1) - ai(j, 1) * di_p(j - 1, 1)) / (bi(j, 1) - ai(j, 1) * ci_p(j - 1, 1));
		}

		fn(M, 1) = di_p(M, 1);
		for (int j = (M - 1); j >= 1; j--)
		{
			fn(j, 1) = di_p(j, 1) - ci_p(j, 1) * fn(j + 1, 1);
		}

	}
	return fn;
}

Matrix extrapolation()
{
	double N1 = 1, N2 = 2, N3 = 3, N4 = 4, N5 = 5, N6 = 6, N7 = 7;

	Matrix fn_11 = implicit_euler(N1);

	Matrix fn_21 = implicit_euler(N2);
	Matrix fn_22 = fn_21 + (fn_21 - fn_11) / ((N2 / N1) - 1);

	Matrix fn_31 = implicit_euler(N3);
	Matrix fn_32 = fn_31 + (fn_31 - fn_21) / ((N3 / N2) - 1);
	Matrix fn_33 = fn_32 + (fn_32 - fn_22) / ((N3 / N1) - 1);

	Matrix fn_41 = implicit_euler(N4);
	Matrix fn_42 = fn_41 + (fn_41 - fn_31) / ((N4 / N3) - 1);
    Matrix fn_43 = fn_42 + (fn_42 - fn_32) / ((N4 / N2) - 1);
	Matrix fn_44 = fn_43 + (fn_43 - fn_33) / ((N4 / N1) - 1);

	Matrix fn_51 = implicit_euler(N5);
	Matrix fn_52 = fn_51 + (fn_51 - fn_41) / ((N5 / N4) - 1);
	Matrix fn_53 = fn_52 + (fn_52 - fn_42) / ((N5 / N3) - 1); 
	Matrix fn_54 = fn_53 + (fn_53 - fn_43) / ((N5 / N2) - 1); 
	Matrix fn_55 = fn_54 + (fn_54 - fn_44) / ((N5 / N1) - 1);

	Matrix fn_61 = implicit_euler(N6);
	Matrix fn_62 = fn_61 + (fn_61 - fn_51) / ((N6 / N5) - 1);
	Matrix fn_63 = fn_62 + (fn_62 - fn_52) / ((N6 / N4) - 1); 
	Matrix fn_64 = fn_63 + (fn_63 - fn_53) / ((N6 / N3) - 1); 
	Matrix fn_65 = fn_64 + (fn_64 - fn_54) / ((N6 / N2) - 1);
	Matrix fn_66 = fn_65 + (fn_65 - fn_55) / ((N6 / N1) - 1);

	Matrix fn_71 = implicit_euler(N7);
	Matrix fn_72 = fn_71 + (fn_71 - fn_61) / ((N7 / N6) - 1); 
	Matrix fn_73 = fn_72 + (fn_72 - fn_62) / ((N7 / N5) - 1); 
	Matrix fn_74 = fn_73 + (fn_73 - fn_63) / ((N7 / N4) - 1); 
	Matrix fn_75 = fn_74 + (fn_74 - fn_64) / ((N7 / N3) - 1); 
	Matrix fn_76 = fn_75 + (fn_75 - fn_65) / ((N7 / N2) - 1); 
	Matrix fn_77 = fn_76 + (fn_76 - fn_66) / ((N7 / N1) - 1);

	return fn_77;
}

double error_n(int N1, int N2)
{
	double error = 0.0;
	Matrix CN1 = crank_nicolson(N1);
	Matrix CN2 = crank_nicolson(N2);

	for (int i = 1; i <= M; i++)
	{
		if (abs(CN1(i, 1) - CN2(i, 1)) > error)
		{
			error = abs(CN1(i, 1) - CN2(i, 1));
		}
	}
	return error;
}

double error_m(Matrix a, Matrix b)
{
	double error = 0.0;
	Matrix matrix_1 = a;
	Matrix matrix_2 = b;

	for (int i = 1; i <= M; i++)
	{
		if (abs(matrix_1(i, 1) - matrix_2(i, 1)) > error)
		{
			error = abs(matrix_1(i, 1) - matrix_2(i, 1));
		}
	}
	return error;
}


int main(int argc, char* argv[])
{
	int N = 950;
	// double delta_t = T / N;
	Matrix IE = implicit_euler(N);
	Matrix CN = crank_nicolson(N);
	cout << IE(441, 1) << endl;
	cout << CN << endl;
	cout << "max error: " << error_n(357, 950) << endl;
	cout << "max error: " << error_m(extrapolation(), crank_nicolson(N)) << endl;

}