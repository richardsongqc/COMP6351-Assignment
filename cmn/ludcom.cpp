#include "Matrix.h"
#include "Vector.h"


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