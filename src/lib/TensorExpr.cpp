#include "TensorExpr.hpp"

namespace Sundance {

Expr contractQuadraticForm(const Array<Array<Expr> >& Q, 
  const Expr& x, const Expr& y)
{
  Expr rtn = 0.0;
  
  TEUCHOS_TEST_FOR_EXCEPT(Q.size() != x.size());
  for (int i=0; i<x.size(); i++) {
    TEUCHOS_TEST_FOR_EXCEPT(Q[i].size() != y.size());
    for (int j=0; j<y.size(); j++) {
      rtn = rtn + Q[i][j]*x[i]*y[j];
    }
  }

  return rtn;
}

Expr contractCubicForm(const Array<Array<Array<Expr> > > & T, const Expr& x, 
  const Expr& y, const Expr& z)
{
  
  Expr rtn = 0.0;
  
  TEUCHOS_TEST_FOR_EXCEPT(T.size() != x.size());
  for (int i=0; i<x.size(); i++) {
    TEUCHOS_TEST_FOR_EXCEPT(T[i].size() != y.size());
    for (int j=0; j<y.size(); j++) {
      TEUCHOS_TEST_FOR_EXCEPT(T[i][j].size() != z.size());
      for (int k=0; k<z.size(); k++) {
        rtn = rtn + T[i][j][k]*x[i]*y[j]*z[k];
      }
    }
  }

  return rtn;
}

Array<Array<Expr> > makeZeroMatrix(int m, int n)
{
  Array<Array<Expr> > rtn(m);
  for (int i=0; i<m; i++) {
    rtn[i] = Array<Expr>(n);
    for (int j=0; j<n; j++) rtn[i][j]=0.0;
  }

  return rtn;
}

Array<Array<Expr> > diagonalMatrix(const Expr& d)
{
  int N = d.size();
  Array<Array<Expr> > rtn(N);
  for (int i=0; i<N; i++) {
    rtn[i] = Array<Expr>(N);
    for (int j=0; j<N; j++) {
      if (i == j) rtn[i][j] = d[i];
      else rtn[i][j] = 0.0;
    }
  }

  return rtn;
}

Array<Array<Array<Expr> > > makeZeroCube(int m, int n, int ell)
{
  Array<Array<Array<Expr> > > rtn(m);
  for (int i=0; i<m; i++) {
    rtn[i] = makeZeroMatrix(n, ell);
  }
  return rtn;
}

Expr transpose(const Expr& x) {
  
  int m = x.size();
  int n = -1;
  for (int i=0; i<m; i++) {
    int k = x[i].size();
    if (n == -1) {
      n = k;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(n != k, RuntimeError, "irregular list structure");     
  }   
        
  Array<Array<Expr> > tmp = makeZeroMatrix(n, m);
  Array<Expr> rtn(n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      tmp[i][j] = x[j][i];
    }
    rtn[i] = new ListExpr(tmp[i]);
  }

  return new ListExpr(rtn);
  
}   


}
