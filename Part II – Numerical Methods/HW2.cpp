#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <random>
#include <math.h>
#include <time.h>
using namespace std;

//unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
//std::default_random_engine generator(seed);


double max(float a, float b) {
	return (b < a) ? a : b;
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z) * (z / 2.0));
	double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
}

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& q,
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + (r - q) * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	// double d2 = (log(S / K) + (r - q) * time) / (sigma * time_sqrt) - 0.5 * sigma * time_sqrt;
	return S * exp(-q * time) * N(d1); //- K * exp(-r * time) * N(d2);
}

double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

double calculate_St(double S, double sigma, double r, double q, double T)
{
	return S * exp((r - q - 0.5 * pow(sigma, 2)) * T + sigma * pow(T, 0.5) * gaussrand());
}

double St_or_0(double St, double K)
{
	if (St > K)
		return St;
	else
		return 0;
}

double calculate_bstar(double S, double sigma, double r, double q, double T, double K)
{
	double X[1000];
	double Y[1000];
	double X_sum = 0; 
	double Y_sum = 0;

	for (int i = 0; i < 1000; i++)
	{
		X[i] = calculate_St(S, sigma, r, q, T);
		Y[i] = St_or_0(calculate_St(S, sigma, r, q, T), K);
		X_sum += X[i];
		Y_sum += Y[i];
	}

	double X_mean = X_sum / 1000;
	double Y_mean = Y_sum / 1000;
	double cov_XY = 0;
	double var_X = 0;

	for (int i = 0; i < 1000; i++)
	{
		cov_XY += (X[i] - X_mean) * (Y[i] - Y_mean);
		var_X += pow((X[i] - X_mean), 2);
	}

	return cov_XY / var_X;
}

double calculate_Xt(double sigma, double mu_hat, double T)
{
	return mu_hat + sigma * pow(T, 0.5) * gaussrand();
}

double calculate_fg(double sigma, double r, double q, double T, double mu_hat, double x)
{
	double mu = (r - q - 0.5 * pow(sigma, 2)) * T;
	return exp((pow(mu, 2) - pow(mu_hat, 2) + 2 * (mu_hat - mu) * x) / (-2 * pow(sigma, 2) * T));
}

double ln_or_0(double x, double K, double S0)
{
	if (x > log(K / S0))
		return 1;
	else
		return 0;
}

int main()
{
	double price_sum, price, ci;
	price_sum = 0;

	double risk_free_rate = 0.005;
	double strike_price = 2200;
	double expiration_time = 1.0 / 12;
	double initial_stock_price = 2000;
	double volatility = 0.3;
	double dividend_yield = 0.02;
	int sample_size = 1000;
	double mu_hat = 0.13;

	double Yib = 0, Yib_2 = 0, se;
	double b_star;
	double Xi = 0, Yi = 0;
	double control;
	//double Xt, fg, St, se;
	//double X;
	//double x_bar = 0, x_bar2 = 0;

	double start, end, time;
	srand(5);

	start = clock();

	// b_star = calculate_bstar(initial_stock_price, volatility, risk_free_rate, dividend_yield, expiration_time, strike_price);

	for (int i = 1; i <= sample_size; i++)
	{
		b_star = calculate_bstar(initial_stock_price, volatility, risk_free_rate, dividend_yield, expiration_time, strike_price);
		Xi = calculate_St(initial_stock_price, volatility, risk_free_rate, dividend_yield, expiration_time);
		Yi = St_or_0(Xi, strike_price); 
		control = Yi + b_star * (initial_stock_price * exp((risk_free_rate - dividend_yield) * expiration_time) - Xi);
		Yib = (1.0 - 1.0 / i) * Yib + (1.0 / i) * control;
		Yib_2 = (1.0 - 1.0 / i) * Yib_2 + (1.0 / i) * pow(control, 2);

		//Xt = calculate_Xt(volatility, mu_hat, expiration_time);
		//fg = calculate_fg(volatility, risk_free_rate, dividend_yield, expiration_time, mu_hat, Xt);
		//St = initial_stock_price * exp(Xt);
		//X = St * ln_or_0(Xt, strike_price, initial_stock_price) * fg;
		//x_bar = (1.0 - 1.0 / i) * x_bar + (1.0 / i) * X;
		//x_bar2 = (1.0 - 1.0 / i) * x_bar2 + (1.0 / i) * pow(X, 2);

	}

	end = clock();

	Yib = exp(-risk_free_rate * expiration_time) * Yib;
	Yib_2 = exp(-2 * risk_free_rate * expiration_time) * Yib_2;

	price = Yib;
	se = pow((1.0 / (sample_size - 1.0)) * (Yib_2 - pow(Yib, 2)), 0.5);

	//x_bar = exp(-risk_free_rate * expiration_time) * x_bar;
	//x_bar2 = exp(-2 * risk_free_rate * expiration_time) * x_bar2;

	//price = x_bar;
	//se = pow((1.0 / (sample_size - 1.0)) * (x_bar2 - pow(x_bar, 2)), 0.5);


	time = double(end - start) / CLOCKS_PER_SEC;

	double price_BS = option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, dividend_yield, volatility, expiration_time);

	cout << "Call Price according to Black-Scholes = " << price_BS << endl;
	cout << "Sample size = " << sample_size << endl;
	cout << "The estimated price = " << price << endl;
	cout << "The estimated standard error = " << se << endl;
	cout << "The computational time in seconds = " << time << "s" << endl;
	cout << "The efficiency measure = " << pow(se, 2) * time << endl;
}
