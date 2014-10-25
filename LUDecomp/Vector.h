#pragma once

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <vector>


class CVector 
{
private:
	int m_maxSize;
	std::vector<double> m_vec;

public:
	CVector();
	CVector(int MaxSize);
	CVector(const CVector& vec);
	~CVector();

	CVector& operator =(const CVector& vec);				// Assignment operator 
	// access to the elements of the CMatrix                     
	double& operator [](int nPos);
	const double& operator [](int nPos) const;

	int CheckSize(int nSize);
	// get Rows and Columns of CMatrix
	int GetSize() const;
};

#endif //__VECTOR_H__