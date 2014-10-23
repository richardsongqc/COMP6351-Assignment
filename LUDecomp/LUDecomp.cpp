#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"
#include <fstream>
#include <iomanip>

// General LU-Decomposition 
void LUDecomp(
	CMatrix a, 
	CVector b, 
	CMatrix& l, 
	CMatrix& u, 
	CVector& x)
{
	CVector  z(b.GetSize());
	int i = 0, j = 0, k = 0, n = b.GetSize();

	// computations of L and U matrices
	for (i = 0; i < n; i++)
	{
		l[i][0] = a[i][0];
	}

	for (j = 1; j < n; j++)
	{
		u[0][j] = a[0][j] / l[0][0];
	}

	for (i = 0; i < n; i++)
	{
		u[i][i] = 1;
	}

	for (i = 1; i < n; i++)
	{
		for (j = 1; j < n; j++)
		{
			if (i >= j)
			{
				l[i][j] = a[i][j];
				for (k = 0; k <= j - 1; k++)
				{
					l[i][j] -= l[i][k] * u[k][j];
				}
			}
			else
			{
				u[i][j] = a[i][j];
				for (k = 0; k <= i - 1; k++)
				{
					u[i][j] -= l[i][k] * u[k][j];
				}
				u[i][j] /= l[i][i];
			}
		}
	}

	// solve lz=b by forward substitution
	z[0] = b[0] / l[0][0];
	for (i = 1; i < n; i++)
	{
		z[i] = b[i];
		for (j = 0; j <= i - 1; j++)
		{
			z[i] -= l[i][j] * z[j];
		}
		z[i] /= l[i][i];
	}

	//solve Ux=z by backward substitution 
	//x[n] = z[n];
	x = z;
	for (i = n - 1; i >= 0; i--)
	{
		x[i] = z[i];
		for (j = i + 1; j < n; j++)
		{
			x[i] -= u[i][j] * x[j];
		}
	}
}

// Tridiagnonal Systems (Matrix Form)
void LUDecompTridiagnoal(
	CMatrix a, 
	CVector b, 
	CVector& x, CVector& g)
{
	int i = 0, j = 0, k = 0, n = b.GetSize();
	CVector l1(b.GetSize());		//	vector a
	CVector l2(b.GetSize());		//	vector b
	CVector l3(b.GetSize());		//	vector c

	// Get the elements and put them into the vector a, b, c
	for (i = 0; i < n; i++)
	{
		if (i > 0)
		{
			l3[i-1] = a[i][i-1];
		}

		if (i + 1 < n)
		{
			l1[i + 1] = a[i + 1][i];
		}

		l2[i] = a[i][i];
	}

	CVector alpha(b.GetSize());
	CVector beta(b.GetSize());
	CVector gama(b.GetSize());
	

	// Apply the algorithm of Gauss elimination ( Page 37 )
	beta[0] = l2[0];
	g[0]    = b[0];

	for (i = 1; i < n; i++)
	{
		gama[i] = l1[i] / beta[i-1];
		beta[i] = l2[i] - gama[i] * l3[i-1];
		g[i] = b[i] - gama[i]*g[i-1];
	}

	// backsubstitution algorithm ( Page 38 )
	x[n - 1] = g[n - 1] / beta[n - 1];
	for (i = n - 2; i >= 0; i-- )
	{
		x[i] = (g[i] - l3[i] * x[i + 1]) / beta[i];
	}

	return;
}

// Tridiagnonal Systems (Vector Form)
void LUDecompTridiagnoal(
	CVector l1, CVector l2, CVector l3,							// vector a, b, c in Tridiagnonal matrix
	CVector b,													
	//CMatrix& l, CMatrix& u, 
	CVector& x, CVector& g)
{
	int i = 0, j = 0, k = 0, n = b.GetSize();

	CVector alpha(b.GetSize());
	CVector beta(b.GetSize());
	CVector gama(b.GetSize());


	// Apply the algorithm of Gauss elimination ( Page 37 )
	beta[0] = l2[0];
	g[0] = b[0];

	for (i = 1; i < n; i++)
	{
		gama[i] = l1[i] / beta[i - 1];
		beta[i] = l2[i] - gama[i] * l3[i - 1];
		g[i] = b[i] - gama[i] * g[i - 1];
	}

	//// computation of L matrix and U matrix
	//for (i = 0; i < n; i++)
	//{
	//	l[i][i] = 1;
	//	u[i][i] = beta[i];

	//	if (i + 1 < n)
	//	{
	//		l[i + 1][i] = gama[i + 1];
	//		u[i][i + 1] = l3[i];
	//	}
	//}

	// backsubstitution algorithm ( Page 38 )
	x[n - 1] = g[n - 1] / beta[n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		x[i] = (g[i] - l3[i] * x[i + 1]) / beta[i];
	}

	return;
}


void main()
{
	int i = 0, j = 0, k = 0;
	const double PI  = 3.1415926535897932384626433832795;
	const double PI2 = 9.8696044010893586188344909998762;

	std::ofstream ofile;              
	ofile.open("d:\\myfile.csv");     

	for (i = 10; i < 1000000; i*=2)
	{
		int nDim = i - 1;
		CVector l1(nDim);
		CVector l2(nDim);
		CVector l3(nDim);
		//CMatrix l(nDim, nDim);
		//CMatrix u(nDim, nDim);
		CVector b(nDim), x(nDim), g(nDim);

		// Prepare the vectors of the discretized boundary value problem
		double d2 = pow(i, 2) / PI2;
		double d1 = -1 - 2 * d2;
		for (j = 0; j < nDim; j++)
		{
			l2[j] = d1;
		

			if (j > 0)
			{
				l3[j - 1] = d2;
			}

			if (j + 1 < nDim)
			{
				l1[j + 1] = d2;
			}

			b[j] = -5 * sin(2 * (j+1) * PI / i);
		}

		// Apply LU-Decomp to solve it
		LUDecompTridiagnoal(l1, l2, l3, b, /*l, u, */x, g);

		// Compute the maximum error
		double dblMaxError = 0;

		for (j = 0; j < nDim; j++)
		{
			double temp = abs(x[j] - sin(2 * (j + 1) * PI / i));
			if (temp > dblMaxError)
			{
				dblMaxError = temp;
			}
		}
		
		double dbl = 4;
		dbl /= 3;
		dbl /= d2;
		// Output the results to one excel file, later we can analyze this file.
		ofile << std::setw(6) << i << "," << std::setw(15) << dblMaxError << "," << std::setw(15) << dbl << std::endl;
		printf("%d\n", i);
	}

	ofile.close();

	printf("\n");
}