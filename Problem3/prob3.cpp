#include <stdio.h>
#include "..\\cmn\\ludcom.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

class MyException
{
};


double GetIntegral(CVector v, int N)
{
	double dbl = 0;
	double dblMeshSize = (double)1.0 / (double)N;
	for (int i = 0; i < v.GetSize(); i++)
	{
		dbl += v[i] * dblMeshSize;
	}

	return dbl;

}

void main()
{
	int N = 40, i = 0;
	CVector u(N + 1), bu(N+1), x(N+1), y(N+1);
	double dblLamda = 0;

	std::ofstream ofile;
	ofile.open("prob3.csv");

	ofile << ",";

	for (i = 0; i < N; i++)
	{
		ofile << (double)i / (double)N << ",";
	}

	ofile << "1" << std::endl;

	// Prepare the vectors of the discretized boundary value problem
	double dbl = pow(N, 2);

	x[0] = 0.0000000000	 ;
	x[1] = 0.0343082496	 ;
	x[2] = 0.0686164992	 ;
	x[3] = 0.1029247489	 ;
	x[4] = 0.1372329985	 ;
	x[5] = 0.1630599908	 ;
	x[6] = 0.1888869832	 ;
	x[7] = 0.2147139755	 ;
	x[8] = 0.2405409678	 ;
	x[9] = 0.2606269044	 ;
	x[10] = 0.2807128409 ;
	x[11] = 0.3007987775 ;
	x[12] = 0.3208847140 ;
	x[13] = 0.3435678045 ;
	x[14] = 0.3662508950 ;
	x[15] = 0.3889339855 ;
	x[16] = 0.4116170759 ;
	x[17] = 0.4311921132 ;
	x[18] = 0.4507671504 ;
	x[19] = 0.4703421876 ;
	x[20] = 0.4899172249 ;
	x[21] = 0.5095079643 ;
	x[22] = 0.5290987038 ;
	x[23] = 0.5486894433 ;
	x[24] = 0.5682801828 ;
	x[25] = 0.5880030858 ;
	x[26] = 0.6077259888 ;
	x[27] = 0.6274488918 ;
	x[28] = 0.6471717948 ;
	x[29] = 0.6704400385 ;
	x[30] = 0.6937082823 ;
	x[31] = 0.7169765260 ;
	x[32] = 0.7402447698 ;
	x[33] = 0.7722219345 ;
	x[34] = 0.8041990992 ;
	x[35] = 0.8361762639 ;
	x[36] = 0.8681534286 ;
	x[37] = 0.9011150715 ;
	x[38] = 0.9340767143 ;
	x[39] = 0.9670383572 ;
	x[40] = 1.0000000000 ;

	y[0 ] = 0.0000000000000;
	y[1 ] = 0.2823891844000;
	y[2 ] = 0.5616385780200;
	y[3 ] = 0.8367352184100;
	y[4 ] = 1.1063627597000;
	y[5 ] = 1.3047245315000;
	y[6 ] = 1.4981548973000;
	y[7 ] = 1.6856021186000;
	y[8 ] = 1.8658348095000;
	y[9 ] = 2.0001139772000;
	y[10] = 2.1284244625000;
	y[11] = 2.2499505967000;
	y[12] = 2.3638180656000;
	y[13] = 2.4820436474000;
	y[14] = 2.5879537854000;
	y[15] = 2.6801791150000;
	y[16] = 2.7574051180000;
	y[17] = 2.8110651920000;
	y[18] = 2.8519926112000;
	y[19] = 2.8796584800000;
	y[20] = 2.8936924388000;
	y[21] = 2.8938958517000;
	y[22] = 2.8802500390000;
	y[23] = 2.8529420916000;
	y[24] = 2.8123377648000;
	y[25] = 2.7585628902000;
	y[26] = 2.6925236742000;
	y[27] = 2.6150004774000;
	y[28] = 2.5268465560000;
	y[29] = 2.4103991867000;
	y[30] = 2.2818895009000;
	y[31] = 2.1427664169000;
	y[32] = 1.9944047168000;
	y[33] = 1.7777560044000;
	y[34] = 1.5489690335000;
	y[35] = 1.3105196672000;
	y[36] = 1.0644562450000;
	y[37] = 0.8045977215000;
	y[38] = 0.5398528351600;
	y[39] = 0.2713590271900;
	y[40] = 0.0000000000000;

	u[0] = 0;
	for (i = 0; i < N; i++)
	{
		double k = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
		double dblX = (double)(i + 1) / (double)N;
		if (dblX < x[i + 1])
		{
			u[i + 1] = y[i] + k*(dblX - x[i]);
		}
		else if (dblX == x[i + 1])
		{
			u[i + 1] = y[i + 1];
		}
		else if (dblX > x[i + 1])
		{
			u[i + 1] = y[i] + k*(dblX - x[i]);
		}
	}







	double deltaLamda = 0.001;
	for (dblLamda = 2; dblLamda < 5; dblLamda += deltaLamda)
	{
		int k = 0;
		int nDim = N - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		CVector b(nDim), deltaU(nDim), g(nDim);

		printf("%.15f) \n", dblLamda);

		try
		{
			
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
						u1 = u[i + 1];
						u2 = u[i + 2];
					}
					else if (i + 1 == nDim)
					{
						u0 = u[nDim - 1];
						u1 = u[nDim];
						u2 = u[N];
					}

					b[i] = -dbl* (u0 - 2 * u1 + u2) - dblLamda*(exp(u0) + 10 * exp(u1) + exp(u2)) / ((double)12);

					double d0 = dbl + dblLamda*(exp(u0) / ((double)12));
					double d1 = -2 * dbl + 10 * dblLamda* (exp(u1) / ((double)12));
					double d2 = dbl + dblLamda* (exp(u2) / ((double)12));

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
					u[i + 1] = deltaU[i] + u[i + 1];
					if (_finite(u[i + 1]) == 0)
					{
						throw MyException();
					}
					//printf("deltaU[%d]=%15.11f\tu[%d]=%12.11f\n", i, deltaU[i], i, u[i+1] );
					// Output the results to one excel file, later we can analyze this file.
					//ofile << std::setw(25) << u[i + 1] << ",";
				}

				if (deltaU < 0.00000000001)
				{
					printf("\nk = %d\n", k);
					break;
				}

				//printf("\n");

				k++;
			} 
			while (k < 1000);
		}
		catch (MyException & e)
		{
			//cout << "*****************************************************************" << endl;
			//cout << "\t!!!Lamda = " << dblLamda << "!!!" << endl;
			//dblLamda -= deltaLamda;
			//deltaLamda /= (double)10;
			//u = bu;
			//if (deltaLamda < pow(0.1,15))
			//{
			//	printf("\n============================================================\nLamda = %.15f\n", dblLamda);
			//	break;
			//}
			continue;
		}

		printf("\nLamda = %.15f \t delta = %.15f\n", dblLamda, deltaLamda);
		ofile << std::setprecision(17) << dblLamda << ",";
		for (i = 0; i < N; i++)
		{
			ofile << std::setw(15) << u[i] << ",";
		}
		ofile << "0," << GetIntegral(u,N) << std::endl;

		bu = u;
	}

	ofile.close();

	printf("\n");
}