/*
A simple matrix class
c++ code
Author: Jos de Jong, Nov 2007. Updated March 2010

With this class you can:
  - create a 2D matrix with custom size
  - get/set the cell values
  - use operators +, -, *, /
  - use functions Ones(), Zeros(), Eye(), Diag(), Det(), Inv(), Size(),
    Transpose(), FlipLR(), FlipUD(), Hcat(), Vcat(), ToGF2()
  - print the content of the matrix

Usage:
  you can create a matrix by:
    Matrix A;
    Matrix A = Matrix(rows, cols);
    Matrix A = B;

  you can get and set matrix elements by:
    A(2,3) = 5.6;    // set an element of Matix A
    value = A(3,1);   // get an element of Matrix A
    value = A.get(3,1); // get an element of a constant Matrix A
    A = B;        // copy content of Matrix B to Matrix A

  you can apply operations with matrices and doubles:
    A = B + C;
    A = B - C;
    A = -B;
    A = B * C;
    A = B / C;

  the following functions are available:
    A = Ones(rows, cols);
    A = Zeros(rows, cols);
    A = Eye(n);
    A = Diag(B);
    d = Det(A);
    A = Inv(B);
    A = Transpose(B);
    A = FlipLR(B);
    A = FlipUD(B);
    A = Hcat(B, C);
    A = Vcat(B, C);
    A = ToGF2(B);
    cols = A.GetCols();
    rows = A.GetRows();
    cols = Size(A, 1);
    rows = Size(A, 2);

  you can quick-print the content of a matrix in the console with:
    A.Print();
    A.PrintGF2();
*/

#include <cstdlib>
#include <cstdio>
#include <math.h>

#define PAUSE {printf("Press \"Enter\" to continue\n"); fflush(stdin); getchar(); fflush(stdin);}

// Declarations
class Matrix;
double Det(const Matrix& a);
Matrix Eye(const int n);
Matrix Diag(const Matrix& v);
Matrix Inv(const Matrix& a);
Matrix Transpose(const Matrix& a);
Matrix FlipLR(const Matrix& a);
Matrix FlipUD(const Matrix& a);
Matrix Hcat(const Matrix& a, const Matrix& b);
Matrix Vcat(const Matrix& a, const Matrix& b);
Matrix ToGF2(const Matrix& a);
Matrix Ones(const int rows, const int cols);
int Size(const Matrix& a, const int i);
Matrix Zeros(const int rows, const int cols);


/*
 * a simple exception class
 * you can create an exeption by entering
 *   throw Exception("...Error description...");
 * and get the error message from the data msg for displaying:
 *   err.msg
 */
class Exception
{
public:
  const char* msg;
  Exception(const char* arg)
   : msg(arg)
  {
  }
};



class Matrix
{
public:
  // constructor
  Matrix()
  {
    //printf("Executing constructor Matrix() ...\n");
    // create a Matrix object without content
    p = NULL;
    rows = 0;
    cols = 0;
  }

  // constructor
  Matrix(const int row_count, const int column_count)
  {
    // create a Matrix object with given number of rows and columns
    p = NULL;

    if (row_count > 0 && column_count > 0)
    {
      rows = row_count;
      cols = column_count;

      p = new double*[rows];
      for (int r = 0; r < rows; r++)
      {
        p[r] = new double[cols];

        // initially fill in zeros for all values in the matrix;
        for (int c = 0; c < cols; c++)
        {
          p[r][c] = 0;
        }
      }
    }
  }

  // assignment operator
  Matrix(const Matrix& a)
  {
    rows = a.rows;
    cols = a.cols;
    p = new double*[a.rows];
    for (int r = 0; r < a.rows; r++)
    {
      p[r] = new double[a.cols];

      // copy the values from the matrix a
      for (int c = 0; c < a.cols; c++)
      {
        p[r][c] = a.p[r][c];
      }
    }
  }

  // index operator. You can use this class like myMatrix(col, row)
  // the indexes are one-based, not zero based.
  double& operator()(const int r, const int c)
  {
    if (p != NULL && r > 0 && r <= rows && c > 0 && c <= cols)
    {
      return p[r-1][c-1];
    }
    else
    {
      throw Exception("Subscript out of range");
    }
  }

  // index operator. You can use this class like myMatrix.get(col, row)
  // the indexes are one-based, not zero based.
  // use this function get if you want to read from a const Matrix
  double get(const int r, const int c) const
  {
    if (p != NULL && r > 0 && r <= rows && c > 0 && c <= cols)
    {
      return p[r-1][c-1];
    }
    else
    {
      throw Exception("Subscript out of range");
    }
  }

  // assignment operator
  Matrix& operator= (const Matrix& a)
  {
    rows = a.rows;
    cols = a.cols;
    p = new double*[a.rows];
    for (int r = 0; r < a.rows; r++)
    {
      p[r] = new double[a.cols];

      // copy the values from the matrix a
      for (int c = 0; c < a.cols; c++)
      {
        p[r][c] = a.p[r][c];
      }
    }
    return *this;
  }

  // add a double value (elements wise)
  Matrix& Add(const double v)
  {
    for (int r = 0; r < rows; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        p[r][c] += v;
      }
    }
     return *this;
  }

  // subtract a double value (elements wise)
  Matrix& Subtract(const double v)
  {
    return Add(-v);
  }

  // multiply a double value (elements wise)
  Matrix& Multiply(const double v)
  {
    for (int r = 0; r < rows; r++)
    {
      for (int c = 0; c < cols; c++)
      {
        p[r][c] *= v;
      }
    }
     return *this;
  }

  // divide a double value (elements wise)
  Matrix& Divide(const double v)
  {
     return Multiply(1/v);
  }

  // addition of Matrix with Matrix
  friend Matrix operator+(const Matrix& a, const Matrix& b)
  {
    // check if the dimensions match
    if (a.rows == b.rows && a.cols == b.cols)
    {
      Matrix res(a.rows, a.cols);

      for (int r = 0; r < a.rows; r++)
      {
        for (int c = 0; c < a.cols; c++)
        {
          res.p[r][c] = a.p[r][c] + b.p[r][c];
        }
      }
      return res;
    }
    else
    {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // addition of Matrix with double
  friend Matrix operator+ (const Matrix& a, const double b)
  {
    Matrix res = a;
    res.Add(b);
    return res;
  }
  // addition of double with Matrix
  friend Matrix operator+ (const double b, const Matrix& a)
  {
    Matrix res = a;
    res.Add(b);
    return res;
  }

  // subtraction of Matrix with Matrix
  friend Matrix operator- (const Matrix& a, const Matrix& b)
  {
    // check if the dimensions match
    if (a.rows == b.rows && a.cols == b.cols)
    {
      Matrix res(a.rows, a.cols);

      for (int r = 0; r < a.rows; r++)
      {
        for (int c = 0; c < a.cols; c++)
        {
          res.p[r][c] = a.p[r][c] - b.p[r][c];
        }
      }
      return res;
    }
    else
    {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // subtraction of Matrix with double
  friend Matrix operator- (const Matrix& a, const double b)
  {
    Matrix res = a;
    res.Subtract(b);
    return res;
  }
  // subtraction of double with Matrix
  friend Matrix operator- (const double b, const Matrix& a)
  {
    Matrix res = -a;
    res.Add(b);
    return res;
  }

  // operator unary minus
  friend Matrix operator- (const Matrix& a)
  {
    Matrix res(a.rows, a.cols);

    for (int r = 0; r < a.rows; r++)
    {
      for (int c = 0; c < a.cols; c++)
      {
        res.p[r][c] = -a.p[r][c];
      }
    }

    return res;
  }

  // operator multiplication
  friend Matrix operator* (const Matrix& a, const Matrix& b)
  {
    // check if the dimensions match
    if (a.cols == b.rows)
    {
      Matrix res(a.rows, b.cols);

      for (int r = 0; r < a.rows; r++)
      {
        for (int c_res = 0; c_res < b.cols; c_res++)
        {
          for (int c = 0; c < a.cols; c++)
          {
            res.p[r][c_res] += a.p[r][c] * b.p[c][c_res];
          }
        }
      }
      return res;
    }
    else
    {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // multiplication of Matrix with double
  friend Matrix operator* (const Matrix& a, const double b)
  {
    Matrix res = a;
    res.Multiply(b);
    return res;
  }
  // multiplication of double with Matrix
  friend Matrix operator* (const double b, const Matrix& a)
  {
    Matrix res = a;
    res.Multiply(b);
    return res;
  }

  // division of Matrix with Matrix
  friend Matrix operator/ (const Matrix& a, const Matrix& b)
  {
    // check if the dimensions match: must be square and equal sizes
    if (a.rows == a.cols && a.rows == a.cols && b.rows == b.cols)
    {
      Matrix res(a.rows, a.cols);

      res = a * Inv(b);

      return res;
    }
    else
    {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // division of Matrix with double
  friend Matrix operator/ (const Matrix& a, const double b)
  {
    Matrix res = a;
    res.Divide(b);
    return res;
  }

  // division of double with Matrix
  friend Matrix operator/ (const double b, const Matrix& a)
  {
    Matrix b_matrix(1, 1);
    b_matrix(1,1) = b;

    Matrix res = b_matrix / a;
    return res;
  }

  /**
   * returns the minor from the given matrix where
   * the selected row and column are removed
   */
  Matrix Minor(const int row, const int col) const
  {
    Matrix res;
    if (row > 0 && row <= rows && col > 0 && col <= cols)
    {
      res = Matrix(rows - 1, cols - 1);

      // copy the content of the matrix to the minor, except the selected
      for (int r = 1; r <= (rows - (row >= rows)); r++)
      {
        for (int c = 1; c <= (cols - (col >= cols)); c++)
        {
          res(r - (r > row), c - (c > col)) = p[r-1][c-1];
        }
      }
    }
    else
    {
      throw Exception("Index for minor out of range");
    }

    return res;
  }

  /*
   * returns the size of the i-th dimension of the matrix.
   * i.e. for i=1 the function returns the number of rows,
   * and for i=2 the function returns the number of columns
   * else the function returns 0
   */
  int Size(const int i) const
  {
    if (i == 1)
    {
      return rows;
    }
    else if (i == 2)
    {
      return cols;
    }
    return 0;
  }

  // returns the number of rows
  int GetRows() const
  {
    return rows;
  }

  // returns the number of columns
  int GetCols() const
  {
    return cols;
  }

  // print the contents of the matrix
  void Print() const
  {
    if (p != NULL)
    {
      printf("[");
      for (int r = 0; r < rows; r++)
      {
        if (r > 0)
        {
          printf(" ");
        }
        for (int c = 0; c < cols-1; c++)
        {
          printf("%.2f, ", p[r][c]);
        }
        if (r < rows-1)
        {
          printf("%.2f;\n", p[r][cols-1]);
        }
        else
        {
          printf("%.2f]\n", p[r][cols-1]);
        }
      }
    }
    else
    {
      // matrix is empty
      printf("[ ]\n");
    }
  }

  // print the contents of the matrix
  void PrintGF2() const
  {
    if (p != NULL)
    {
      printf("[");
      for (int r = 0; r < rows; r++)
      {
        if (r > 0)
        {
          printf(" ");
        }
        for (int c = 0; c < cols-1; c++)
        {
          printf("%d, ", int(p[r][c]));
        }
        if (r < rows-1)
        {
          printf("%d;\n", int(p[r][cols-1]));
        }
        else
        {
          printf("%d]\n", int(p[r][cols-1]));
        }
      }
    }
    else
    {
      // matrix is empty
      printf("[ ]\n");
    }
  }

public:
  // destructor
  ~Matrix()
  {
    // clean up allocated memory
    for (int r = 0; r < rows; r++)
    {
      delete p[r];
    }
    delete p;
    p = NULL;
  }

private:
  int rows;
  int cols;
  double** p;     // pointer to a matrix with doubles
};

/*
 * i.e. for i=1 the function returns the number of rows,
 * and for i=2 the function returns the number of columns
 * else the function returns 0
 */
int Size(const Matrix& a, const int i)
{
  return a.Size(i);
}


/**
 * returns a matrix with size cols x rows with ones as values
 */
Matrix Ones(const int rows, const int cols)
{
  Matrix res = Matrix(rows, cols);

  for (int r = 1; r <= rows; r++)
  {
    for (int c = 1; c <= cols; c++)
    {
      res(r, c) = 1;
    }
  }
  return res;
}

/**
 * returns a matrix with size cols x rows with zeros as values
 */
Matrix Zeros(const int rows, const int cols)
{
  return Matrix(rows, cols);
}


/**
 * returns an identity matrix with size n x n
 * @param  v a vector
 * @return a diagonal matrix with ones on the diagonal
 */
Matrix Eye(const int n)
{
  Matrix res = Matrix(n, n);
  for (int i = 1; i <= n; i++)
  {
    res(i, i) = 1;
  }
  return res;
}

/**
 * returns a diagonal matrix
 * @param  v a vector
 * @return a diagonal matrix with the given vector v on the diagonal
 */
Matrix Diag(const Matrix& v)
{
  Matrix res;
  if (v.GetCols() == 1)
  {
    // the given matrix is a vector n x 1
    int rows = v.GetRows();
    res = Matrix(rows, rows);

    // copy the values of the vector to the matrix
    for (int r=1; r <= rows; r++)
    {
      res(r, r) = v.get(r, 1);
    }
  }
  else if (v.GetRows() == 1)
  {
    // the given matrix is a vector 1 x n
    int cols = v.GetCols();
    res = Matrix(cols, cols);

    // copy the values of the vector to the matrix
    for (int c=1; c <= cols; c++)
    {
      res(c, c) = v.get(1, c);
    }
  }
  else
  {
    throw Exception("Parameter for diag must be a vector");
  }
  return res;
}

/*
 * returns the determinant of Matrix a
 * Note: This implementation users the "minors" method
 * which is inefficient for large matrices. The method
 * requires a whopping N!/2 (factorial(N) / 2) computations
 * of 2-by-2 sub-determinants. As a ersult the function
 * now thorws on any input Matrix larger than 10-by-10.
 * TODO: Implement using a more efficient method such as
 * LU decomposition.
 */
double Det(const Matrix& a)
{
  double d = 0;    // value of the determinant
  int rows = a.GetRows();
  int cols = a.GetRows();

  if (rows == cols)
  {
    if (rows > 10)
    {
      throw Exception("Det(): Max size (10-by-10) exceeded.");
    }
    // this is a square matrix
    if (rows == 1)
    {
      // this is a 1 x 1 matrix
      d = a.get(1, 1);
    }
    else if (rows == 2)
    {
      // this is a 2 x 2 matrix
      // the determinant of [a11,a12;a21,a22] is det = a11*a22-a21*a12
      d = a.get(1, 1) * a.get(2, 2) - a.get(2, 1) * a.get(1, 2);
    }
    else
    {
      // this is a matrix of 3 x 3 or larger
      for (int c = 1; c <= cols; c++)
      {
        Matrix M = a.Minor(1, c);
        //d += pow(-1, 1+c) * a(1, c) * Det(M);
        d += (c%2 + c%2 - 1) * a.get(1, c) * Det(M); // faster than with pow()
      }
    }
  }
  else
  {
    throw Exception("Matrix must be square");
  }
  return d;
}

// swap two values
void Swap(double& a, double& b)
{
  double temp = a;
  a = b;
  b = temp;
}

/*
 * returns the inverse of Matrix a
 */
Matrix Inv(const Matrix& a)
{
  Matrix res;
  double d = 0;    // value of the determinant
  int rows = a.GetRows();
  int cols = a.GetRows();

  if (rows == cols)
  {
    // this is a square matrix
    if (rows == 1)
    {
      // this is a 1 x 1 matrix
      res = Matrix(rows, cols);
      res(1, 1) = 1 / a.get(1, 1);
    }
    else if (rows == 2)
    {
      // this is a 2 x 2 matrix
      // d = determinant of a 2-by-2 matrix
      d = a.get(1, 1) * a.get(2, 2) - a.get(2, 1) * a.get(1, 2);
      if (d == 0)
      {
        throw Exception("Determinant of matrix is zero");
      }
      res = Matrix(rows, cols);
      res(1, 1) = a.get(2, 2);
      res(1, 2) = -a.get(1, 2);
      res(2, 1) = -a.get(2, 1);
      res(2, 2) = a.get(1, 1);
      res = (1/d) * res;
    }
    else
    {
      // this is a matrix of 3 x 3 or larger
      // calculate inverse using gauss-jordan elimination
      //   http://mathworld.wolfram.com/MatrixInverse.html
      //   http://math.uww.edu/~mcfarlat/inverse.htm
      res = Eye(rows);   // a diagonal matrix with ones at the diagonal
      Matrix ai = a;    // make a copy of Matrix a

      for (int c = 1; c <= cols; c++)
      {
        // element (c, c) should be non zero. if not, swap content
        // of lower rows
        int r;
        for (r = c; r <= rows && ai(r, c) == 0; r++)
        {
        }
        if (r != c)
        {
          // swap rows
          for (int s = 1; s <= cols; s++)
          {
            Swap(ai(c, s), ai(r, s));
            Swap(res(c, s), res(r, s));
          }
        }

        // eliminate non-zero values on the other rows at column c
        for (int r = 1; r <= rows; r++)
        {
          if (ai(c, c) == 0)
          {
            throw Exception("Matrix is singular.");
          }
          if(r != c)
          {
            // eleminate value at column c and row r
            if (ai(r, c) != 0)
            {
              double f = - ai(r, c) / ai(c, c);

              // add (f * row c) to row r to eleminate the value
              // at column c
              for (int s = 1; s <= cols; s++)
              {
                ai(r, s) += f * ai(c, s);
                res(r, s) += f * res(c, s);
              }
            }
          }
          else
          {
            // make value at (c, c) one,
            // divide each value on row r with the value at ai(c,c)
            double f = ai(c, c);
            for (int s = 1; s <= cols; s++)
            {
              ai(r, s) /= f;
              res(r, s) /= f;
            }
          }
        }
      }
    }
  }
  else
  {
    if (rows != cols)
    {
      throw Exception("Matrix must be square");
    }
  }
  return res;
}

Matrix Transpose(const Matrix& a)
{
  int rows = a.Size(1);
  int cols = a.Size(2);
  Matrix res = Matrix(cols, rows);
  for (int r=1;  r <= rows; r++)
  {
    for (int c=1;  c <= cols; c++)
    {
      res(c, r) = a.get(r, c);
    }
  }
  return res;
}

Matrix FlipLR(const Matrix& a)
{
  int rows = a.Size(1);
  int cols = a.Size(2);
  Matrix res = Matrix(rows, cols);
  for (int r=1;  r <= rows; r++)
  {
    for (int c=1;  c <= cols; c++)
    {
      res(r, cols - c + 1) = a.get(r, c);
    }
  }
  return res;
}

Matrix FlipUD(const Matrix& a)
{
  int rows = a.Size(1);
  int cols = a.Size(2);
  Matrix res = Matrix(rows, cols);
  for (int r=1;  r <= rows; r++)
  {
    for (int c=1;  c <= cols; c++)
    {
      res(rows - r + 1, c) = a.get(r, c);
    }
  }
  return res;
}

Matrix Hcat(const Matrix& a, const Matrix& b)
{
  int rows_a = a.Size(1);
  int cols_a = a.Size(2);
  int rows_b = b.Size(1);
  int cols_b = b.Size(2);
  if (rows_a != rows_b) {
    throw Exception("Dimensions of arrays being concatenated are not consistent.");
  }
  Matrix res = Matrix(rows_a, cols_a + cols_b);
  for (int r=1;  r <= rows_a; r++)
  {
    for (int c=1;  c <= cols_a; c++)
    {
      res(r, c) = a.get(r, c);
    }
    for (int c=1;  c <= cols_b; c++)
    {
      res(r, cols_a + c) = b.get(r, c);
    }
  }
  return res;
}

Matrix Vcat(const Matrix& a, const Matrix& b)
{
  int rows_a = a.Size(1);
  int cols_a = a.Size(2);
  int rows_b = b.Size(1);
  int cols_b = b.Size(2);
  if (cols_a != cols_b) {
    throw Exception("Dimensions of arrays being concatenated are not consistent.");
  }
  Matrix res = Matrix(rows_a + rows_b, cols_a);
  for (int c=1;  c <= cols_a; c++)
  {
    for (int r=1;  r <= rows_a; r++)
    {
      res(r, c) = a.get(r, c);
    }
    for (int r=1;  r <= rows_b; r++)
    {
      res(rows_a + r, c) = b.get(r, c);
    }
  }
  return res;
}

Matrix ToGF2(const Matrix& a)
{
  int rows = a.Size(1);
  int cols = a.Size(2);
  Matrix res = Matrix(rows, cols);
  for (int r=1;  r <= rows; r++)
  {
    for (int c=1;  c <= cols; c++)
    {
      res(r, c) = double(int(abs(round(a.get(r, c)))) % 2);
    }
  }
  return res;
}
