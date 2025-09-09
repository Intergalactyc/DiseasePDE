//#include "/home/intergalactyc/Code/TTUTrilinos/BUILDS/OPT/include/Sundance.hpp"
#include "Sundance.hpp"

#ifndef tensorexprheader
#define tensorexprheader

namespace Sundance {

  /* Do the contraction Q_{i,j} x_i y_j */
  Expr contractQuadraticForm(const Array<Array<Expr> > & Q, 
    const Expr& x, const Expr& y) ;
  
  /* Do the contraction T_{i,j} x_i y_j z_k */
  Expr contractCubicForm(const Array<Array<Array<Expr> > > & T, 
    const Expr& x, const Expr& y, const Expr& z);

  /* Construct a matrix of all zeros. */
  Array<Array<Expr> > makeZeroMatrix(int m, int n);

  /* Construct a tensor cube of all zeros. */
  Array<Array<Array<Expr> > > makeZeroCube(int m, int n, int ell);


  /* Transpose an order 2 tensor */
  Expr transpose(const Expr& x);

  /* Construct a diagonal matrix based on a list of diagonal entries. */
  Array<Array<Expr> > diagonalMatrix(const Expr& d);
}

#endif
