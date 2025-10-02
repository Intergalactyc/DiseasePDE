#include "Sundance.hpp"
#include "SIR.hpp"

Expr grad = gradient(2);

Expr f_weak(Expr U, Expr UHat, Expr x, Expr y, Expr t, ModelParams p, bool isTest)
{
  Expr S = U[0];
  Expr I = U[1];
  Expr R = U[2];
  Expr SHat = UHat[0];
  Expr IHat = UHat[1];
  Expr RHat = UHat[2];
  return p.D_S*(grad*SHat)*(grad*S)
    + p.D_I*(grad*IHat)*(grad*I)
    + p.D_R*(grad*RHat)*(grad*R)
    + p.mu*SHat*S + p.beta*SHat*I*S - p.ell*SHat*R - p.Lambda*SHat
    + (p.mu+p.w+p.gamm)*IHat*I - p.beta*IHat*(I*S)
    + (p.mu+p.ell)*RHat*R - p.gamm*RHat*I
    - isTest*resid(x,y,t,p)*(UHat);
}

struct State {
  Expr U;
  Expr UHat;
  Expr UPrev;
  Expr t;
  Expr tPrev;
  ModelParams p;
  Expr x;
  Expr y;
  double dt;
  bool isTest;
};

Expr fwe(State s)
{
  return s.UHat*(s.U-s.UPrev)
    + s.dt * f_weak(s.UPrev, s.UHat, s.x, s.y, s.t, s.p, s.isTest);
}

Expr bwe(State s)
{
  return s.UHat*(s.U-s.UPrev)
    + s.dt * f_weak(s.U, s.UHat, s.x, s.y, s.t, s.p, s.isTest);
}

Expr itr(State s)
{
  return s.UHat*(s.U-s.UPrev)
    + s.dt * 0.5 * (f_weak(s.UPrev, s.UHat, s.x, s.y, s.tPrev, s.p, s.isTest) + f_weak(s.U, s.UHat, s.x, s.y, s.t, s.p, s.isTest));
}

Expr rk4(State s)
{
  // this won't work anyway b/c the u_i are arguments of weak forms, not actual U values
  Expr k1 = f_weak(s.U, s.UHat, s.x, s.y, s.tPrev, s.p, s.isTest);
  Expr u1 = s.UHat * (s.U-s.UPrev)
    + 0.5 * s.dt * k1;
  Expr k2 = f_weak(u1, s.UHat, s.x, s.y, s.tPrev+0.5*s.dt, s.p, s.isTest);
  Expr u2 = s.UHat * (s.U-s.UPrev)
    + 0.5 * s.dt * k2;
  Expr k3 = f_weak(u2, s.UHat, s.x, s.y, s.tPrev+0.5*s.dt, s.p, s.isTest);
  Expr u3 = s.UHat * (s.U-s.UPrev)
    + 0.5 * s.dt * k3;
  Expr k4 = f_weak(u3, s.UHat, s.x, s.y, s.t, s.p, s.isTest);
  return s.UHat*(s.U-s.UPrev)
    + s.dt * (1./6.) * (k1 + k2 + k3 + k4);
}

int main(int argc, char** argv)
{
  try
  {
    /* Run tests? */
    bool isTest = false;

    /* Define number of time steps and final time */
    int nSteps = 64;
    double T_final = 2.0;
    int nx = 32;

    std::string paramFile = "base-params.xml";
    std::string solverFile = "playa-newton-amesos.xml";
    std::string outputLocation = "Results/SIR";
    std::string method = "itr";

    /* Handle command-line options */
    Sundance::setOption("paramFile", paramFile, "XML file containing parameters");
    Sundance::setOption("solver", solverFile, "Name of XML file for solver");
    Sundance::setOption("out", outputLocation, "Prefix for VTU file output");
    Sundance::setOption("nx", nx, "Number of elements in each spatial direction");
    Sundance::setOption("nt", nSteps, "Number of timesteps");
    Sundance::setOption("tf", T_final, "Final time");
    Sundance::setOption("method", method, "Numerical method");
    Sundance::setOption("test", "release", isTest, "Flag to run test");

    ModelParams p;

    if(paramFile != "") {
      // Attempt to parse parameter XML file using Teuchos
      try{
        FileInputSource paramInput = FileInputSource(paramFile);
        XMLObject paramXML = paramInput.getObject();

        // Use parameter values from the XML file
        p.D_S = xml_parameter(paramXML, "Diffusion", "susceptible", p.D_S);
        p.D_I = xml_parameter(paramXML, "Diffusion", "infectious", p.D_I);
        p.D_R = xml_parameter(paramXML, "Diffusion", "recovered", p.D_R);
        p.mu = xml_parameter(paramXML, "Population", "base_mortality", p.mu);
        p.Lambda = xml_parameter(paramXML, "Population", "spontaneous_creation", p.Lambda);
        p.beta = xml_parameter(paramXML, "Disease", "infection_rate", p.beta);
        p.gamm = xml_parameter(paramXML, "Disease", "recovery_rate", p.gamm);
        p.w = xml_parameter(paramXML, "Disease", "excess_mortality", p.w);
        p.ell = xml_parameter(paramXML, "Disease", "immunity_loss", p.ell);

        Out::root() << "Parameter input file parsed successfully\n";
      } catch(...) {
        Out::root() << "Error parsing parameter input file\n";
      }
    }

    /* Initialize */
    Sundance::init(&argc, &argv);

    double dt = T_final/((double) nSteps);

    std::map<std::string, std::string> methodName{{"bwe","Backward Euler"},{"itr","Implicit Trapezoidal"}};
    Out::root() << "Running simulation.\n";
    Out::root() << "Method selected: " << methodName[method] << "\n";
    Out::root() << "Time steps: " << nSteps << " (to time " << T_final << ", dt=" << dt << ")\n";
    Out::root() << "Spatial elements: " << nx << "\n";

    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    /* Create/load a mesh. */
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, 
        0.0, 1.0, nx, meshType);
    Mesh mesh = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
      
    /* Create unknown and test functions, discretized using first-order
     * Lagrange interpolants */
    BasisFamily bas = new Lagrange(1);
    Expr S = new UnknownFunction(bas, "S");
    Expr I = new UnknownFunction(bas, "I");
    Expr R = new UnknownFunction(bas, "R");
    Expr SHat = new TestFunction(bas, "SHat");
    Expr IHat = new TestFunction(bas, "IHat");
    Expr RHat = new TestFunction(bas, "RHat");
    Expr U = List(S, I, R);
    Expr UHat = List(SHat, IHat, RHat);

    /* Create coordinate function */
    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    
    /* Initial profile*/
    Expr UStart;
    if (isTest){
      UStart = List(cos(pi*x)*cos(3*pi*y),
        4*cos(pi*x)*cos(2*pi*y),
        0.0);
    } else{
      double epsilon = 0.05;
      double delta = 0.2;
      UStart = List(1.0,
        4*epsilon/(pi*delta*delta)*exp(-(x*x+y*y)/delta/delta),
        0.0);
    }
    Expr SStart = UStart[0];
    Expr IStart = UStart[1];
    Expr RStart = UStart[2];

    /* Represent the time variable as a parameter expression, NOT as
     * a double variable. The reason is that we need to be able to update
     * the time value without rebuilding expressions. */
    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);

    /* Project onto P1 basis to form UPrev */
    DiscreteSpace discSpace(mesh, List(bas, bas, bas), vecType);
    L2Projector projector(discSpace, UStart);
    Expr UPrev = projector.project();

    Expr SPrev = UPrev[0];
    Expr IPrev = UPrev[1];
    Expr RPrev = UPrev[2];

    // Current Newton approximation
    Expr UNewt = copyDiscreteFunction(UPrev, "SIRNewt");

    QuadratureFamily quad = new GaussianQuadrature(4);

    /* Define the weak form, semidiscretized in time */
    Expr weak;
    State s;
    s.U = U;
    s.UHat = UHat;
    s.UPrev = UPrev;
    s.t = t;
    s.tPrev = tPrev;
    s.p = p;
    s.x = x;
    s.y = y;
    s.dt = dt;
    s.isTest = isTest;
    if(method == "bwe"){
      // Backward Euler
      weak = bwe(s);
    } else if (method == "itr") {
      // Implicit Trapezoidal
      weak = itr(s);
    } else if (method == "rk4") {
      // 4th-order Runge-Kutta
      weak = rk4(s);
    } else {
      throw std::invalid_argument("solution method '" + method + "' not recognized.");
    }
    Expr eqn = Integral(interior, weak, quad);

    /* There are no Dirichlet BCs */
    Expr bc;

    /* We can now set up the nonlinear problem! */
    NonlinearProblem prob(mesh, eqn, bc, UHat, U, UNewt, vecType);

    NonlinearSolver<double> solver 
      = NonlinearSolverBuilder::createSolver(solverFile);

    /* Write the initial conditions */
    {
      FieldWriter writer = new VTKWriter(outputLocation + "-0"); 
      writer.addMesh(mesh);
      writer.addField("S", new ExprFieldWrapper(UPrev[0]));
      writer.addField("I", new ExprFieldWrapper(UPrev[1]));
      writer.addField("R", new ExprFieldWrapper(UPrev[2]));
      writer.write();
    }

    /* Loop over timesteps */
    double maxErr = 0.0;
    for (int i=0; i<nSteps; i++)
    {
      Out::root() << "timestep #" << i << endl;
      t.setParameterValue((i+1)*dt);
      tPrev.setParameterValue(i*dt);

      SolverState<double> state = prob.solve(solver);

      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
        std::runtime_error,
        "Nonlinear solve failed to converge: message=" << state.finalMsg());
          
      updateDiscreteFunction(UNewt, UPrev);
      if (isTest) {
        double err = L2Norm(mesh,interior,UNewt-UExact(x,y,t),quad);
        if (err > maxErr){
          maxErr = err;
        }
      }
      
      FieldWriter writer = new VTKWriter(outputLocation + "-" 
        + Teuchos::toString(i+1));
      writer.addMesh(mesh);
      writer.addField("S", new ExprFieldWrapper(UPrev[0]));
      writer.addField("I", new ExprFieldWrapper(UPrev[1]));
      writer.addField("R", new ExprFieldWrapper(UPrev[2]));
      writer.write();
    }
    
    if (isTest) {
      double tol = 0.01;
      Sundance::passFailTest(maxErr, tol);
    }

  }
  catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}

// xml_parameter functions to attempt discovering an attribute value in the parameter file
template <class T>
T xml_parameter(XMLObject sourceObj, std::string param_type, std::string attribute, T default_value){
  int child_id = sourceObj.findFirstChild(param_type.append("Parameters"));
  if (child_id != -1){
    XMLObject childObj = sourceObj.getChild(child_id);
    if (childObj.hasAttribute(attribute)){
      return T(std::stod(childObj.getAttribute(attribute)));
    }
  }
  return default_value;
}

// Exact solution
Expr UExact(const Expr& x, const Expr& y, const Expr& t){
  return List(cos(pi*x)*cos(3*pi*y)*(1 + sin(t)/2.),
        cos(pi*x)*cos(2*pi*y)*(2 + cos(2*t) + exp(-t)),
        cos(2*pi*x)*cos(pi*y)*(1 - exp(-t)));
}

// Residual
Expr resid(const Expr& x, const Expr& y, const Expr& t, const ModelParams p){
  return List(-p.Lambda + 
        (cos(t)*cos(pi*x)*cos(3*pi*y))/
          2. - p.ell*cos(2*pi*x)*cos(pi*y)*
          (1 - exp(-t)) + 
          p.mu*cos(pi*x)*cos(3*pi*y)*
          (1 + sin(t)/2.) + 
          10*p.D_S*cos(pi*x)*cos(3*pi*y)*
          pow(pi,2)*(1 + sin(t)/2.) + 
          p.beta*cos(2*pi*y)*cos(3*pi*y)*
          (2 + cos(2*t) + exp(-t))*
          pow(cos(pi*x),2)*
          (1 + sin(t)/2.),
        (p.gamm + p.mu + p.w)*cos(pi*x)*
          cos(2*pi*y)*
          (2 + cos(2*t) + exp(-t)) + 
          5*p.D_I*cos(pi*x)*cos(2*pi*y)*
          (2 + cos(2*t) + exp(-t))*
          pow(pi,2) - 
          p.beta*cos(2*pi*y)*cos(3*pi*y)*
          (2 + cos(2*t) + exp(-t))*
          pow(cos(pi*x),2)*
          (1 + sin(t)/2.) + 
          cos(pi*x)*cos(2*pi*y)*
          (-exp(-t) - 2*sin(2*t)),
        (p.ell + p.mu)*cos(2*pi*x)*cos(pi*y)*
          (1 - exp(-t)) + 
          cos(2*pi*x)*cos(pi*y)*exp(-t) - 
          p.gamm*cos(pi*x)*cos(2*pi*y)*
          (2 + cos(2*t) + exp(-t)) + 
          5*p.D_R*cos(2*pi*x)*cos(pi*y)*
          (1 - exp(-t))*pow(pi,2));
}
