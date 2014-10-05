#pragma once

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>

typedef std::vector<std::vector<double> >    MATRIXELEMENTTYPE;

class CMatrix
{
private:
	int m_row;
	int m_col;												// rows and columns of CMatrix
	MATRIXELEMENTTYPE   m_vec;        

public:
	// constructors
	CMatrix();
	CMatrix(int row, int col);						
	CMatrix(const CMatrix& mat);                   
	~CMatrix();												// Destructor        
	CMatrix& operator =(const CMatrix& mat);				// Assignment operator  
	CMatrix& operator +=(const CMatrix& mat);

	// access to the elements of the CMatrix                     
	std::vector<double>& operator [](int row);
	const std::vector<double>& operator [](int row) const;


	int CheckRowCol(int row, int col);
	// get Rows and Columns of CMatrix
	int GetRows() const;           
	int GetCols() const;           
};

#endif //__MATRIX_H__