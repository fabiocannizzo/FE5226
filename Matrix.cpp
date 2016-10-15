#include <iostream>
#include <cstring>

using namespace std;

class Matrix
{
    double *m_data;
    unsigned m_nRows, m_nCols;

    // auxiliary function to convert from (row,col) coordinates to the
    // position in the array m_data
    unsigned index(unsigned r, unsigned c) const { return r*m_nCols+c; }

public:
   // default constructor
   Matrix() : m_data(NULL), m_nRows(0), m_nCols(0) {}

   // another constructor
   Matrix(unsigned nr, unsigned nc, double v )
        : m_data(new double[nr*nc])  // allocate memory
        , m_nRows(nr)
        , m_nCols(nc)
        {
            // set every element in the array to v
            // we do the loop using pointer arithmetic, but we could have
            // also used the usual square brackets
            unsigned nElementsToCopy = nr*nc;
            const double *pEnd = m_data+nElementsToCopy;
            for (double *p = m_data; p != pEnd; ++p)
                *p = v;
        }

    // copy constructor
    // we do not use the default copy constructor because we do not want
    // to copy the pointer. We want to create a new storage space and copy
    // in it the content of the other matrix
   Matrix(const Matrix& m)
        : m_data(new double[m.m_nRows*m.m_nCols])
        , m_nRows(m.m_nRows)
        , m_nCols(m.m_nCols)
        {
            // we use the crt function memcpy
            unsigned nBytesToCopy = m.m_nRows*m.m_nCols*sizeof(double);
            memcpy(m_data, m.m_data, nBytesToCopy);
        }

    // destructor
    ~Matrix()
    {
        delete [] m_data;  // release memory
    }

    // accessor methods
    unsigned nRows() const { return m_nRows; }
    unsigned nColumns() const { return m_nCols; }

    // non const (read-write) implementation of operator()
    double& operator()(unsigned r, unsigned c)
    {
        return m_data[index(r,c)];
    }

    // const version (read only) implementation of operator()
    const double& operator()(unsigned r, unsigned c) const
    {
        return m_data[index(r,c)];
    }

    // assignment operator
    Matrix& operator=(const Matrix& m)
    {
        // we check that the user is not trying to assign the matrix to itself
        if (this != &m) {
            unsigned nElements = m.m_nRows*m.m_nCols;
            delete [] m_data; // release memory (note that if nElements<currentsize, we could skip memory reallocation)
            m_data = NULL; // for safety, set the the pointer to NULL, in case something goes wrong in the allocation
            m_data = new double [nElements];  // allocate memory
            m_nRows = m.m_nRows;
            m_nCols = m.m_nCols;
            unsigned nBytesToCopy = nElements*sizeof(double);
            memcpy( m_data, m.m_data, nBytesToCopy);
        }
        return *this;
    }
};

int main()
{
    // create a 3 by 2 matrix and initialize all elements to 0.0
    Matrix m1(3, 2, 0.0);
    Matrix m2(m1);  // create copy of the matrix via copy constructor

    // requires methods nRows() and nColumns()
    for (unsigned i=0, n = min(m2.nRows(),m2.nColumns()); i < n; ++i)
       m2(i,i) = 1;   // requires operator()(int,int)

    Matrix m3;  // invoke default constructor, we set m_data=0
    m3 = m2;    // create copy of the matrix via assignment operator
    m3 = m3;    // is your assignment operator smart enough to handle this?

    const Matrix& cm = m3;
    for (unsigned i=0; i < cm.nRows(); ++i) {
        for (unsigned j=0; j < cm.nColumns(); ++j)
            cout << cm(i,j) << " "; // requires operator()(int,int) const
        cout << endl;
    }

    return 0;
}