#include <stdio.h>
#include "..\\cmn\\ludcom.h"
#include <fstream>
#include <iomanip>

void main()
{
	int N = 0, i = 0;
	//const double PI = 3.1415926535897932384626433832795;
	//const double PI2 = 9.8696044010893586188344909998762;

	//std::ofstream ofile;
	//ofile.open("d:\\prob2_2.csv");

	for (N = 5; N < 1000000; N++)
	{
		int k = 0;
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), deltaU(nDim), g(nDim), u(N+1);

		// Prepare the vectors of the discretized boundary value problem
		double dbl = pow(N, 2);
		for (i = 0; i < N; i++)
		{
			u[i] = (double)i / (double)N;
		}

		do
		{
			u[N] = 0;

			for (i = 0; i < nDim; i++)
			{
				double u0 = 0;  
				double u1 = 0;  
				double u2 = 0;  

				if (i == 0)
				{
					u0 = u[i];
					u1 = u[i + 1];
					u2 = u[i + 2];  
				}
				else if (i + 2 <= nDim)
				{
					u0 = u[i];  
					u1 = u[i+ 1];
					u2 = u[i + 2];
				}
				else if (i + 1 == nDim)
				{
					u0 = u[nDim-1];  
					u1 = u[nDim];
					u2 = u[N];
				}
			           
				b[i] = -dbl* (u0 - 2 * u1 + u2) - (exp(u0) + 10 * exp(u1) + exp(u2)) / 12;

				double d0 = dbl + exp(u0) / 12;
				double d1 = -2*dbl + 10*exp(u1)/12;
				double d2 = dbl+exp(u2)/12;

				l2[i] = d1;

				if (i > 0)
				{
					l3[i - 1] = d2;
				}

				if (i + 1 < nDim)
				{
					l1[i + 1] = d0;
				}
			}

			// Apply LU-Decomp to solve it
			LUDecompTridiagnoal(l1, l2, l3, b, deltaU, g);

			//printf("\nk = %d\n", k);
			for (i = 0; i < nDim; i++)
			{
				u[i+1] = deltaU[i] + u[i+1];
				//printf("deltaU[%d]=%15.11f\tu[%d]=%12.11f\n", i, deltaU[i], i, u[i] );
			}

			if (deltaU < 0.00000000001)
			{
				printf("\nk = %d\n", k);
				break;
			}

			printf("\n");

			k++;
		} 
		while (k < 1000  );

		printf("%d) \n", N);
	}

	//ofile.close();

	printf("\n");
}