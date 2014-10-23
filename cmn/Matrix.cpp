#include "Matrix.h"
#include <vector>

using namespace std;

CMatrix::CMatrix() :
	m_row(1),
	m_col(1)
{
	m_vec.assign(m_row, vector<double>(m_col));
}

// two constructors
	
CMatrix::CMatrix(int row, int col) :
	m_row(row),
	m_col(col)
{

	m_vec.assign(row, vector<double>(col));
}

	
CMatrix::CMatrix( const CMatrix& mat)
{
	m_row = mat.m_row;
	m_col = mat.m_col;
	m_vec.assign(mat.m_row, vector<double>(mat.m_col));
	for( int i = 0; i < m_row ; ++i )
	{
		for ( int j = 0; j < m_col; ++j )
		{
			m_vec[i][j] = mat.m_vec[i][j];
		}
	}
}

// mathematical operation
CMatrix& CMatrix::operator =(const CMatrix& mat)
{
	for ( int i = 0; i < m_row; ++i) 
	{
		for ( int j = 0; j < m_col; ++j) 
		{
			m_vec[i][j] = mat.m_vec[i][j];
		}
	}

	return *this;
}

// Destructor
CMatrix::~CMatrix()
{
	m_row = 0;
	m_col = 0;
	m_vec.clear();
}



CMatrix& CMatrix::operator +=(const CMatrix& mat)            
{                                                                   
	for (int i = 0; i < m_row; ++i)
	{
		for (int j = 0; j < m_col; ++j)
		{
			m_vec[i][j] = mat.m_vec[i][j];
		}
	}

	return *this;
}

// access to the elements of the CMatrix                     
std::vector<double>& CMatrix::operator [](int row)
{
	return m_vec[row];
}

const std::vector<double>& CMatrix::operator [](int row) const
{
	return m_vec[row];
}


int CMatrix::GetRows() const 
{ 
	return m_row; 
}

int CMatrix::GetCols() const 
{ 
	return m_col; 
}

int CMatrix::CheckRowCol(int row, int col)
{
	return ( ( row >= 0 ) && ( row < m_row ) && ( col >= 0 ) && ( col < m_col ) ) ? 1 : 0;
}