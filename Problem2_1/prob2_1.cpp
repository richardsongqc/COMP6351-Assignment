#include <stdio.h>
#include "..\\cmn\\ludcom.h"
#include <fstream>
#include <iomanip>

void main()
{
	int N = 0, i = 0, k = 0;
	const double PI = 3.1415926535897932384626433832795;
	const double PI2 = pow(PI,2);

	std::ofstream ofile;
	ofile.open("prob2_1.csv");

	for (N = 5; N < 1000000; N *=2 )
	{
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), u(nDim), y(N+1), f(N + 1);

		// Prepare the vectors of the discretized boundary value problem
		double d = 2 * PI / N;
		for (i = 0; i < N+1; i++)
		{
			y[i] = sin(i*d);
			f[i] = -5 * y[i];
			//printf("f[%d] = %11.8f\n", i, f[i]);
		}

		y[N] = 0;
		f[N] = 0;
		double d2 = pow(N, 2) / PI2;
		double d1 = -2 * d2 - (double)10/(double)12;
		double d0 = d2 - (double)1 /(double)12;
		for (i = 0; i < nDim; i++)
		{
			l2[i] = d1;

			if (i > 0)
			{
				l3[i - 1] = d0;
			}

			if (i + 1 < nDim)
			{
				l1[i + 1] = d0;
			}

			b[i] = (f[i] + 10 * f[i + 1] + f[i + 2]) / 12;
			
			//printf("b[%d] = (%11.8f + 10 * %11.8f + %11.8f)/12 = %11.8f\n", i, f[i], f[i + 1], f[i + 2], b[i]);
		}

		// Apply LU-Decomp to solve it
		LUDecompTridiagnoal(l1, l2, l3, b, u);

		// Compute the maximum error
		double dblMaxError = 0;

		//std::vector<double> vecError;
		//vecError.clear();

		for (i = 0; i < nDim; i++)
		{
			//y[i+1] = sin((i)*d);

			double temp = abs(y[i + 1]- u[i]);///*b[i] -f[i+1]*//*u[i] - y[i]*/);
			//vecError.push_back(temp);
			//printf("u[%6d] = %11.8f\t y[%6d] = %11.8f\n", i, u[i], i + 1, y[i+1]);
			if (temp > dblMaxError)
			{
				dblMaxError = temp;
			}
		}

		// Output the results to one excel file, later we can analyze this file.
		ofile << std::setw(6) << N << "," << std::setw(25) << dblMaxError << std::endl;
		
		printf("%7d) %-15.11f\n", N, dblMaxError);
	}

	ofile.close();

	printf("\n");
}