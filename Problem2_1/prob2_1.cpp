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
	ofile.open("d:\\prob2_1.csv");

	for (N = 5; N < 1000000; N++)
	{
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), x(nDim), g(nDim), f(N + 1);

		// Prepare the vectors of the discretized boundary value problem
		for (i = 0; i < N; i++)
		{
			f[i] = -5 * sin(2 * i * PI / N);
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
		}

		// Apply LU-Decomp to solve it
		LUDecompTridiagnoal(l1, l2, l3, b, x, g);

		//printf("%d) ", N);
		//for (i = 0; i < nDim; i++)
		//{
		//	printf("%10f, ", x[i]);
		//}
		//printf("\n");

		// Compute the maximum error
		double dblMaxError = 0;

		for (i = 0; i < nDim; i++)
		{
			double y = 0;
			if (i == 0)
			{
				y = ( 10 * x[0] + x[1] ) / 12;
			}
			else if (i + 2 < nDim)
			{
				y = ( x[i] + 10 * x[i + 1] + x[i + 2] ) / 12;
			}
			else if (i +1 == nDim)
			{
				y = (x[i - 1] + 10 * x[i]) / 12;
			}

			double temp = abs( x[i] - y );
			if (temp > dblMaxError)
			{
				dblMaxError = temp;
			}
		}

		// Output the results to one excel file, later we can analyze this file.
		ofile << std::setw(6) << N << "," << std::setw(15) << dblMaxError << std::endl;
		
		printf("%d) \n", N);
	}

	ofile.close();

	printf("\n");
}