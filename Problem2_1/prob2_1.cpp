#include <stdio.h>
#include "..\\cmn\\ludcom.h"
#include <fstream>
#include <iomanip>

void main()
{
	int N = 0, i = 0, k = 0;
	const double PI = 3.1415926535897932384626433832795;
	const double PI2 = 9.8696044010893586188344909998762;

	std::ofstream ofile;
	ofile.open("prob2_1.csv");

	for (N = 10; N < 1000000; N *= 2)
	{
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), u(nDim), g(nDim), y(nDim), f(N + 1);

		// Prepare the vectors of the discretized boundary value problem
		double d = 2 * PI / N;
		for (i = 0; i < N+1; i++)
		{
			
			f[i] = -5 * sin(i*d);
			//printf("%f\n", f[i]);
		}

		f[N] = 0;
		double d2 = pow(N, 2) / PI2;
		double d1 = -2 * d2 - 10 / 12;
		double d0 = d2 - 1 / 12;
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
			
			//printf("%11.8f + 10 * %11.8f + %11.8f = b[%d] = %11.8f\n", f[i], f[i + 1], f[i + 2], i, b[i]);
		}

		// Apply LU-Decomp to solve it
		LUDecompTridiagnoal(l1, l2, l3, b, u, g);

		// Compute the maximum error
		double dblMaxError = 0;

		for (i = 0; i < nDim; i++)
		{
			y[i] = sin((i + 1)*d);
			//double dbl = 0;
			//if (i == 0)
			//{
			//	dbl = (10 * x[i] + x[i + 1]) / 12;
			//}
			//else if (i + 2 <= nDim)
			//{
			//	dbl = (x[i-1] + 10 * x[i] + x[i + 1]) / 12;
			//}
			//else if (i + 1 == nDim)
			//{
			//	dbl = (x[i-1] + 10 * x[i] ) / 12;
			//}
			//
			//y[i] = dbl;
			
			double temp = abs(b[i] -f[i+1]/*u[i] - y[i]*/);
			if (temp > dblMaxError)
			{
				dblMaxError = temp;
			}
		}


		// Output the results to one excel file, later we can analyze this file.
		ofile << std::setw(6) << N << "," << std::setw(20) << dblMaxError << std::endl;
		
		printf("%d) %f\n", N, dblMaxError);
	}

	ofile.close();

	printf("\n");
}