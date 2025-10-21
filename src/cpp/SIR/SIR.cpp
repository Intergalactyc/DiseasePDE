#include "Sundance.hpp"
#include "SIR.hpp"

#include "exodusII.h"

const int nComp = 3;

Expr grad = gradient(2);

Expr f_weak(Expr U, Expr UHat, ModelParams p)
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
    + (p.mu+p.ell)*RHat*R - p.gamm*RHat*I;
}

Expr fwe(Expr U, Expr UHat, Expr UPrev, ModelParams p, double dt)
{
  return UHat*(U-UPrev)
    + dt * f_weak(UPrev, UHat, p);
}

Expr bwe(Expr U, Expr UHat, Expr UPrev, ModelParams p, double dt)
{
  return UHat*(U-UPrev)
    + dt * f_weak(U, UHat, p);
}

Expr itr(Expr U, Expr UHat, Expr UPrev, ModelParams p, double dt)
{
  return UHat*(U-UPrev)
    + dt * 0.5 * (f_weak(U, UHat, p) + f_weak(UPrev, UHat, p));
}

int main(int argc, char** argv)
{
  try
  {
    /* Define number of time steps and final time */
    int nSteps = 64;
    double T_final = 2.0; // Rescale this to make sense with data timescale (data is daily with T=k being k days since start; won't see much change in 1 day step)

    std::string meshFile = "";
    std::string paramFile = "base-params.xml"; // Want to update to make more sense (especially update recovery rate to match our 8 day assumption in data gen)
    std::string solverFile = "playa-newton-amesos.xml";
    std::string outputLocation = "../../../data_products/simulated/SIR";
    std::string method = "itr";

    /* Handle command-line options */
    Sundance::setOption("meshFile", meshFile, "Mesh file");
    Sundance::setOption("paramFile", paramFile, "XML file containing parameters");
    Sundance::setOption("solver", solverFile, "Name of XML file for solver");
    Sundance::setOption("out", outputLocation, "Prefix for VTU file output");
    Sundance::setOption("nt", nSteps, "Number of timesteps");
    Sundance::setOption("tf", T_final, "Final time");
    Sundance::setOption("method", method, "Numerical method");

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
    ex_opts(EX_VERBOSE); // For debugging purposes (q - is this built in to Sundance::init or anything?)

    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    /* Create/load a mesh. */
    MeshType meshType = new BasicSimplicialMeshType();
    if (meshFile == ""){
      meshFile = "/home/intergalactyc/Code/DiseasePDE/data_products/datameshes/SIR_100";
      //"../../../data_products/datameshes/SIR_0";
    }
    
    Out::root() << "Loading mesh: " << meshFile << ".exo\n";

    MeshSource reader = new ExodusMeshReader(meshFile, meshType); // This exo file has both mesh and data (for all timesteps) - what to do different?
    Mesh mesh = reader.getMesh();
    RCP<Array<Array<double> > > nodeAttrValues;
    RCP<Array<Array<double> > > elemAttrValues;
    reader.getAttributes(nodeAttrValues, elemAttrValues);
    Out::root() << "Mesh loaded successfully!\n";

    // Get values with *nodeAttrValues()

    /* Represent the time variable as a parameter expression, NOT as
     * a double variable. The reason is that we need to be able to update
     * the time value without rebuilding expressions. */
    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);
    double dt = T_final/((double) nSteps);

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
      
    /* Create unknown and test functions, discretized using first-order
     * Lagrange interpolants */
    BasisFamily bas = new Lagrange(0); // Here we use P1; element data needs >=P0, while nodal data needs >=P1
    Expr S = new UnknownFunction(bas, "S");
    Expr I = new UnknownFunction(bas, "I");
    Expr R = new UnknownFunction(bas, "R");
    Expr SHat = new TestFunction(bas, "SHat");
    Expr IHat = new TestFunction(bas, "IHat");
    Expr RHat = new TestFunction(bas, "RHat");
    Expr U = List(S, I, R);
    Expr UHat = List(SHat, IHat, RHat);
    
    /* Create the discrete space */
    DiscreteSpace discSpace(mesh, List(bas, bas, bas), vecType);

    /* Use the data to form a DiscreteFunction representing the initial profile */
    Out::root() << "Forming initial conditions DiscreteFunction with loaded data\n";

    Array<Array<double> > data = *elemAttrValues();

    int nCells = mesh.numCells(2); // 0 for nodes, 1 for edges, 2 for elements

    // Validate that the size of the data is compatible with the mesh
    TEUCHOS_TEST_FOR_EXCEPTION(
      data.size() != nComp,
      std::runtime_error,
      "Data size incompatible with problem (" << data.size() << " components found instead of " << nComp << ")\n"
    );
    for (int i = 0; i < nComp; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        data[i].size() != nCells,
        std::runtime_error,
        "Data size incompatible with mesh (in component " << i << ", found " << data[i].size() << " elements instead of " << nCells << ")"
      );
    }

    Expr UStart = new DiscreteFunction(discSpace, 0.0, "UStart");

    Vector<double> vec = DiscreteFunction::discFunc(UStart)->getVector();
    const RCP<DOFMapBase>& dofMap = DiscreteFunction::discFunc(UStart)->map();
    for (int i = 0; i < nComp; i++) {
      for (int j = 0; j < nCells; j++) {
        Array<int> dofs;
        dofMap->getDOFsForCell(2, j, i, dofs);
        Out::root() << dofs << "\n";
        int dof = dofs[0];
        vec[dof] = data[i][j];
      }
    }
    DiscreteFunction::discFunc(UStart)->setVector(vec);

    Out::root() << UStart << "\n";

    Expr SStart = UStart[0];
    Expr IStart = UStart[1];
    Expr RStart = UStart[2];

    Out::root() << "Initial conditions set!\n";

    Out::root() << "Basis: " << discSpace.basis() << " (size " << discSpace.basis().size() << ")\n"; //Tmp
    Out::root() << "UStart size: " << UStart.size() << "\n"; //Tmp
    Out::root() << "UStart[i] size: " << UStart[0].size() << "\n"; //Tmp

    /* Project onto the P1 basis to form UPrev */
    L2Projector projector(discSpace, UStart);
    Expr UPrev = projector.project();
    Expr SPrev = UPrev[0];
    Expr IPrev = UPrev[1];
    Expr RPrev = UPrev[2];

    // Current Newton approximation
    Expr UNewt = copyDiscreteFunction(UPrev, "SIRNewt");

    /* Use 4th order Gaussian quadrature */
    QuadratureFamily quad = new GaussianQuadrature(4);

    /* Define the weak form, semidiscretized in time, based on solution method */
    Expr weak;
    if (method == "fwe") {
      // Forward Euler
      weak = fwe(U, UHat, UPrev, p, dt);
    } else if (method == "bwe") {
      // Backward Euler
      weak = bwe(U, UHat, UPrev, p, dt);
    } else if (method == "itr") {
      // Implicit Trapezoidal
      weak = itr(U, UHat, UPrev, p, dt);
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

    std::map<std::string, std::string> methodName{{"bwe","Backward Euler"},{"itr","Implicit Trapezoidal"}};
    Out::root() << "Running simulation.\n";
    Out::root() << "Method selected: " << methodName[method] << "\n";
    Out::root() << "Time steps: " << nSteps << " (to time " << T_final << ", dt=" << dt << ")\n";

    /* Write the initial conditions */
    {
      FieldWriter writer = new ExodusWriter(outputLocation + "-0");
      writer.addMesh(mesh);
      writer.addField("S", new ExprFieldWrapper(UPrev[0]));
      writer.addField("I", new ExprFieldWrapper(UPrev[1]));
      writer.addField("R", new ExprFieldWrapper(UPrev[2]));
      writer.write();
    }

    /* Loop over timesteps */
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
      
      FieldWriter writer = new ExodusWriter(outputLocation + "-" 
        + Teuchos::toString(i+1));
      writer.addMesh(mesh);
      writer.addField("S", new ExprFieldWrapper(UPrev[0]));
      writer.addField("I", new ExprFieldWrapper(UPrev[1]));
      writer.addField("R", new ExprFieldWrapper(UPrev[2]));
      writer.write();
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
