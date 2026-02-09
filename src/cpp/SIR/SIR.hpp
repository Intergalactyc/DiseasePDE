#include "Sundance.hpp"

class ModelParams{
  public:
    /* Diffusion rates */
    double D_S=0.005;
    double D_I=0.002;
    double D_R=0.005;
    /* Parameters */
    double ell = 0.4;    // Immunity loss
    double mu = 0.0;      // Base mortality rate
    double beta = 3.0;    // Infection rate
    double w = 0.0;       // Excess mortality rate
    double gamm = 0.6;   // Recovery rate
    double Lambda = 0.0;  // Spontaneous creation
};

template <class T>
T xml_parameter(XMLObject sourceObj, std::string param_type, std::string attribute, T default_value);

Expr UExact(const Expr& x, const Expr& y, const Expr& t);
Expr resid(const Expr& x, const Expr& y, const Expr& t, const ModelParams p);

const double pi = 4.0*atan(1.0);
const double deg2rad = pi/180.;
const double rad2deg = 180./pi;

const double R = 1.; // radius of Earth, arbitrary (only for GlobalCoordinateSystem definition - cancels out everywhere it's used by us)

const int nComp = 3;
