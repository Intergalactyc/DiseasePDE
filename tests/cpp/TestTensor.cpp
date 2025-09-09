#include "TensorExpr.hpp"

using namespace Sundance;

int main(int argc, char** argv) {

  Sundance::init(&argc, &argv);

  try {

    BasisFamily P1 = new Lagrange(1);

    Expr v1 = new TestFunction(P1, "v1");
    Expr v2 = new TestFunction(P1, "v2");
    Expr v3 = new TestFunction(P1, "v3");
    Expr V = List(v1, v2, v3);
    Expr u1 = new UnknownFunction(P1, "u1");
    Expr u2 = new UnknownFunction(P1, "u2");
    Expr u3 = new UnknownFunction(P1, "u3");
    Expr U = List(u1, u2, u3);

    Array<Array<Expr> > A = makeZeroMatrix(3,3);
    A[0][0]=1.0;
    A[0][2]=1.0;
    A[1][0]=0.5;
    A[1][1]=2.0;
    A[2][0]=0.25;
    A[2][2]=3.0;

    Array<Array<Array<Expr> > > T = makeZeroCube(3,3,3);
    T[0][1][2] = 6.0;
    T[1][2][0] = 0.3;
    T[2][0][1] = 0.6; 
  
    Expr w1 = contractQuadraticForm(A, V, U);
    Expr w2 = contractCubicForm(T, V, U, U);

    Expr w = w1 + w2;

    Out::root() << "weak form = " << w << endl;

    Expr d = List(1.0, 2.0, 3.0);
    
    Array<Array<Expr> > D = diagonalMatrix(d);
    Out::root() << "diagonal matrix D= " << D << endl;

    Expr grad = gradient(2);

    Out::root() << "computing gradient " << endl;
    Expr gradV = transpose(outerProduct(grad, V));
    Expr gradU = transpose(outerProduct(grad, U));
    Out::root() << "grad v = " << gradV << endl;
    Expr w3 = contractQuadraticForm(D, gradV, gradU);
    Out::root() << "wf with gradient contraction" << endl << w3 << endl;


  }
  catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
}
