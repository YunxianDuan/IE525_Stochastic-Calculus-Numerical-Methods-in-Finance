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
	double d2 = (log(S / K) + (r - q) * time) / (sigma * time_sqrt) - 0.5 * sigma * time_sqrt;
	return S * exp(-q * time) * N(d1) - K * exp(-r * time) * N(d2);
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

double calculate_St(double S, double sigma, double r, double q, double T, double a)
{
	if (a == 1)
	{
		return S * exp((r - q - 0.5 * pow(sigma, 2)) * T + sigma * pow(T, 0.5) * gaussrand());
	}
	else
	{
		return S * exp((r - q - 0.5 * pow(sigma, 2)) * T - sigma * pow(T, 0.5) * gaussrand());
	}
	// return S * exp((r - q - 0.5 * pow(sigma, 2)) * T + sigma * pow(T, 0.5) * gaussrand());
}

int main()
{
	double price_sum, price, ci;
	price_sum = 0;

	double risk_free_rate = 0.003866;
	double strike_price = 1870;
	double expiration_time = 1.0 / 52;
	double initial_stock_price = 1868.99;
	double volatility = 0.2979;
	double dividend_yield = 0.0232;
	int sample_size = 40000000;

	double x_bar = 0, y_bar = 0, se;
	double upper, lower, z_score = 1.96;

	double start, end, time;
	srand(5);

	start = clock();

	for (int i = 1; i <= sample_size; i++)
	{
		// ci = exp(-risk_free_rate * expiration_time) * (max(calculate_St(initial_stock_price, volatility, risk_free_rate, dividend_yield, expiration_time, 1) - strike_price, 0)); // the standard approach
		ci = 0.5 * exp(-risk_free_rate * expiration_time) * (max(calculate_St(initial_stock_price, volatility, risk_free_rate, dividend_yield, expiration_time, 1) - strike_price, 0) 
			+ (max(calculate_St(initial_stock_price, volatility, risk_free_rate, dividend_yield, expiration_time, 0) - strike_price, 0))); // the antithetic approach
		x_bar = (1.0 - 1.0 / i) * x_bar + (1.0 / i) * ci;
		y_bar = (1.0 - 1.0 / i) * y_bar + (1.0 / i) * pow(ci, 2);
		price_sum += ci;
	}

	end = clock();

	price = price_sum / sample_size;
	se = pow((1.0 / (sample_size - 1.0)) * (y_bar - pow(x_bar, 2)), 0.5);
	upper = x_bar + z_score * se;
	lower = x_bar - z_score * se;
	time = double(end - start) / CLOCKS_PER_SEC;

	double price_BS = option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, dividend_yield, volatility, expiration_time);

	cout << "Call Price according to Black-Scholes = " << price_BS << endl;
	cout << "Sample size = " << sample_size << endl;
	cout << "The estimated price = " << price << endl;
	cout << "The estimated standard error = " << se << endl;
	cout << "The 95% confidence interval = [" << lower << ", " << upper << "]" << endl;
	cout << "The computational time in seconds = " << time << "s" << endl;
	cout << "The efficiency measure = " << pow(se, 2) * time << endl;


}