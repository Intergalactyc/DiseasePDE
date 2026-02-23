#include "Sundance.hpp"
#include "SIR.hpp"
#include <unordered_map>
#include "exodusII.h"

/* Define a global (spherical, surface-constrained) coordinate system */
class GlobalCoordinateSystem : public CoordinateSystemBase
{
public:
  GlobalCoordinateSystem(double R, bool isRadians)
    : longitude_(new CoordExpr(0)),
    latitude_(new CoordExpr(1)),
    R_(R),
    conv_(isRadians ? 1 : deg2rad)
  {}
  // for our purposes any R is okay (constant factor doesn't affect inner product use for projection), and isRadians=false

  Expr jacobian() const {
    return R_ * R_ * cos(conv_*latitude_);
  }

  bool supportsMeshDimension(int dim) const
  {
    return dim == 2;
  }

  std::ostream& print(std::ostream& os) const
  {
    os << "Global (Radius " << R_ << ")";
    return os;
  }

  virtual RCP<CoordinateSystemBase> getRcp()
  {
    return rcp(this);
  }

private:
  Expr longitude_;
  Expr latitude_;
  Expr R_; 
  Expr conv_;
};

// Theta is colatitude (adjust) // phi is -longitude
// correct for degrees, and when considering values with /km^2 units consider scale factor
// can use grad[0] and grad[1] for x/y components
Expr grad = gradient(2);
Expr d_dtheta = grad[0];
Expr d_dphi = -grad[1];

/* RHS of the model weak form  */
Expr f_weak(Expr U, Expr UHat, ModelParams p, Expr longitude, Expr latitude)
{
  Expr S = U[0];
  Expr I = U[1];
  Expr R = U[2];
  Expr SHat = UHat[0];
  Expr IHat = UHat[1];
  Expr RHat = UHat[2];
  // Given lat/lon are in degrees: phi = deg2rad * longitude, theta = deg2rad * (90 - latitude) are az/pol angles
  Expr sine_theta = cos(deg2rad * latitude);
  // deg2rad / rad2deg scale factors cancel out between the d_dtheta and dtheta (d_dphi and dphi) terms 
  return - p.D_S * ((d_dtheta*S)*(d_dtheta*SHat)*sine_theta + (d_dphi*S)*(d_dphi*SHat)/sine_theta)
    - p.D_I * ((d_dtheta*I)*(d_dtheta*IHat)*sine_theta + (d_dphi*I)*(d_dphi*IHat)/sine_theta)
    - p.D_R * ((d_dtheta*R)*(d_dtheta*RHat)*sine_theta + (d_dphi*R)*(d_dphi*RHat)/sine_theta)
    + p.mu*SHat*S + p.beta*SHat*I*S - p.ell*SHat*R - p.Lambda*SHat
    + (p.mu+p.w+p.gamm)*IHat*I - p.beta*IHat*(I*S)
    + (p.mu+p.ell)*RHat*R - p.gamm*RHat*I;
}

/* Full weak form integrand for Forward Euler */
Expr fwe(Expr U, Expr UHat, Expr UPrev, ModelParams p, double dt, Expr longitude, Expr latitude)
{
  return UHat*(U-UPrev)
    + dt * f_weak(UPrev, UHat, p, longitude, latitude);
}

/* Full weak form integrand for Backward Euler */
Expr bwe(Expr U, Expr UHat, Expr UPrev, ModelParams p, double dt, Expr longitude, Expr latitude)
{
  return UHat*(U-UPrev)
    + dt * f_weak(U, UHat, p, longitude, latitude);
}

/* Full weak form integrand for Implicit Trapezoidal Rule */
Expr itr(Expr U, Expr UHat, Expr UPrev, ModelParams p, double dt, Expr longitude, Expr latitude)
{
  return UHat*(U-UPrev)
    + dt * 0.5 * (f_weak(U, UHat, p, longitude, latitude) + f_weak(UPrev, UHat, p, longitude, latitude));
}

/* Helper function for reading strings from Exodus files */
string clean_exo_string(const std::string& s) {
    size_t end = s.find('\0');
    std::string trimmed = s.substr(0, end);
    while (!trimmed.empty() && trimmed.back() == ' ')
        trimmed.pop_back();
    return trimmed;
}

/* Component map for ensuring that components stay in the right order between input, solving, and output */
std::unordered_map<int, int> get_component_map(Array<string> attrNames) {
  std::unordered_map<int, int> comps;
  Array<string> all( {"S", "I", "R"} );
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (all[i] == clean_exo_string(attrNames[j])){
        comps[i] = j;
        break;
      }
    }
  }
  return comps;
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
    std::string solverFile = "playa-newton-amesos.xml";// "aztec-ml.xml"; // "playa-newton-amesos.xml"; // Sparse direct solver
    std::string outputLocation = "../../../data_products/simulated/";
    std::string outputPrefix = "L2SIR";
    std::string method = "itr";

    /* Handle command-line options */
    Sundance::setOption("meshFile", meshFile, "Mesh file");
    Sundance::setOption("paramFile", paramFile, "XML file containing parameters");
    Sundance::setOption("solver", solverFile, "Name of XML file for solver");
    Sundance::setOption("out", outputLocation, "Location of file output directory");
    Sundance::setOption("prefix", outputPrefix, "Prefix of output file names (default is 'SIR')");
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

    /* Load the mesh into memory */
    MeshSource reader = new ExodusMeshReader(meshFile, meshType); // This exo file has both mesh and data (for all timesteps) - what to do different?
    Mesh mesh = reader.getMesh();
    RCP<Array<Array<double> > > nodeAttrValues;
    RCP<Array<Array<double> > > elemAttrValues;
    reader.getAttributes(nodeAttrValues, elemAttrValues);
    RCP<Array<string>> nodeAttrNames;
    RCP<Array<string>> elemAttrNames;
    reader.getAttributeNames(nodeAttrNames, elemAttrNames);
    std::unordered_map<int, int> componentMap = get_component_map(*elemAttrNames());
    Out::root() << "Mesh loaded successfully!\n";

    /* Represent the time variable as a parameter expression, NOT as
     * a double variable. The reason is that we need to be able to update
     * the time value without rebuilding expressions. */
    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);
    double dt = T_final/((double) nSteps);

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
      
    /* Cell-wise (P0) basis and discrete space, for loaded data*/
    BasisFamily p0 = new Lagrange(0);
    DiscreteSpace discCellSpace(mesh, List(p0, p0, p0), vecType);
    
    /* Use the data to form a DiscreteFunction representing the initial profile */
    Out::root() << "Forming initial conditions DiscreteFunction with loaded data\n";

    Array<Array<double> > data = *elemAttrValues();

    int nCells = mesh.numCells(2); // 0 for nodes, 1 for edges, 2 for elements (cells)

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

    Expr cellData = new DiscreteFunction(discCellSpace, 0.0, "UStart");

    Vector<double> vec = DiscreteFunction::discFunc(cellData)->getVector();
    const RCP<DOFMapBase>& dofMap = DiscreteFunction::discFunc(cellData)->map();
    for (int i = 0; i < nComp; i++) {
      int c = componentMap[i];
      for (int j = 0; j < nCells; j++) {
        Array<int> dofs;
        dofMap->getDOFsForCell(2, j, i, dofs);
        int dof = dofs[0];
        vec[dof] = data[c][j];
      }
    }
    DiscreteFunction::discFunc(cellData)->setVector(vec);

    /* Node-wise (P1: first order Lagrange interpolants) basis, for solving */
    BasisFamily p1 = new Lagrange(1);
    DiscreteSpace discNodeSpace(mesh, List(p1, p1, p1), vecType);

    /* Define global coordinate system so projectors account for Jacobian */
    CoordinateSystem coordSystem = new GlobalCoordinateSystem(R, false);

    /* Define projector from P0->P1 and project our initial data onto the P1 basis */
    L2Projector nodeProjector(discNodeSpace, coordSystem, cellData);
    Expr UStart = nodeProjector.project();
    Expr SStart = UStart[0];
    Expr IStart = UStart[1];
    Expr RStart = UStart[2];

    Out::root() << "Initial conditions set!\n";

    /* Create unknown and test functions in node space*/
    Expr S = new UnknownFunction(p1, "S");
    Expr I = new UnknownFunction(p1, "I");
    Expr R = new UnknownFunction(p1, "R");
    Expr SHat = new TestFunction(p1, "SHat");
    Expr IHat = new TestFunction(p1, "IHat");
    Expr RHat = new TestFunction(p1, "RHat");
    Expr U = List(S, I, R);
    Expr UHat = List(SHat, IHat, RHat);

    /* Form UPrev */
    Expr UPrev = copyDiscreteFunction(UStart);
    Expr SPrev = UPrev[0];
    Expr IPrev = UPrev[1];
    Expr RPrev = UPrev[2];

    /* Current Newton approximation */
    Expr UNewt = copyDiscreteFunction(UPrev, "SIRNewt");

    /* Use 4th order Gaussian quadrature */
    QuadratureFamily quad = new GaussianQuadrature(4);

    /* Create latitude/longitude coordinate expressions */
    Expr longitude = new CoordExpr(0); // x coordinates in mesh are longitudes in degrees
    Expr latitude = new CoordExpr(1); // y coordinates in mesh are latitudes in degrees

    /* Define the weak form, semidiscretized in time, based on solution method */
    Expr weak;
    if (method == "fwe") {
      // Forward Euler
      weak = fwe(U, UHat, UPrev, p, dt, longitude, latitude);
    } else if (method == "bwe") {
      // Backward Euler
      weak = bwe(U, UHat, UPrev, p, dt, longitude, latitude);
    } else if (method == "itr") {
      // Implicit Trapezoidal
      weak = itr(U, UHat, UPrev, p, dt, longitude, latitude);
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

    std::map<std::string, std::string> methodName{{"bwe","Backward Euler"},{"itr","Implicit Trapezoidal"},{"fwe","Forward Euler"}};
    Out::root() << "Running simulation.\n";
    Out::root() << "Method selected: " << methodName[method] << "\n";
    Out::root() << "Time steps: " << nSteps << " (to time " << T_final << ", dt=" << dt << ")\n";
    Out::root() << "Mesh size: " << nCells << " cells\n";

    /* Write the initial conditions (from the P0 data) */
    {
      FieldWriter writer = new ExodusWriter(outputLocation + outputPrefix + "_0.0");

      writer.addMesh(mesh);
      writer.addField("S0", new ExprFieldWrapper(cellData[0]));
      writer.addField("I0", new ExprFieldWrapper(cellData[1]));
      writer.addField("R0", new ExprFieldWrapper(cellData[2]));

      writer.addField("S", new ExprFieldWrapper(UStart[0]));
      writer.addField("I", new ExprFieldWrapper(UStart[1]));
      writer.addField("R", new ExprFieldWrapper(UStart[2]));
      writer.write();
    }

    exit(0);

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
      
      // TODO: try to determine how to write the timestep # / associated time value to the Exodus file
      // It seems that in write_exo_mesh (called at some point in the ExodusWriter write),
      // writing of data at multiple timesteps is supported - but how are they to be specified?
      // And how would a different timestep be specified to start at? (It seems that the times will always start at 0, based on 773 & 775...?)
      // /home/intergalactyc/Code/TTUTrilinos/packages/seacas/libraries/exodus/cbind/test/create_mesh.c
      
      /* Project (P1) solution at this step back onto P0 cell-space before saving */
      L2Projector cellProjector(discCellSpace, UPrev);
      Expr cellU = cellProjector.project();
      
      /* Write this step's output to an Exodus file */
      FieldWriter writer = new ExodusWriter(outputLocation + outputPrefix + "_" + Teuchos::toString((i+1)*dt));
      writer.addMesh(mesh);
      writer.addField("S", new ExprFieldWrapper(cellU[0]));
      writer.addField("I", new ExprFieldWrapper(cellU[1]));
      writer.addField("R", new ExprFieldWrapper(cellU[2]));
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
