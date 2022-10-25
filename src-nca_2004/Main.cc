#include <sys/stat.h>
#include <sys/types.h>
#include "parser.h"
#include "timer.h"
#include "util.h"
#include "complex.h"
#include "mesh.h"
#include "function.h"
#include "Common.h"

using namespace std;

namespace Parameters{
  StringPar   Sig        ("Sig",  "sig.inp", "\t\t# The name of the input Self-energies");
  StringPar   Ac         ("Ac",    "Ac.inp", "\t\t# The name of the input bath function");
  StringPar   out        ("out",        ".", "\t\t# The name of the output directory");
  StringPar   gloc       ("gloc","gloc.out", "\t# The name of the output Green's function");
  StringPar   sig        ("sig",  "sig.out", "\t\t# The name of the output Self-energy");
  StringPar   cix        ("cix",  "cix.dat", "\t\t# The name of the input file containing information about bands and their degeneracy");
  Par<double> Ed         ("Ed",         -2., "\t\t# Energy level");
  Par<double> U          ("U",           4., "\t\t\t# Coulomb repulsion");
  Par<double> T          ("T",          0.2, "\t\t# Temperature ");
  Par<double> Q0         ("Q",            1, "\t\t\t# Default average Q in grand-canonical ansamble");
  Par<int>    prt        ("prt",          1, "\t\t# Weather to print intermediate results");
  Par<double> alpha      ("alpha",      0.5, "\t\t# The fraction of the new self energy to be used in the next iteration");
  Par<double> max_diff   ("max_diff",  1e-6, "\t# Criterium to finish the procedure");
  Par<int>    max_steps  ("max_steps",  300, "\t# Criterium to finish the procedure");
  Par<double> StartLambda("StartLambda", -1, "\t# Where to start searching for the lambda0");
  Par<double> EndLambda  ("EndLambda",    1, "\t\t# Where to stop searching for the lambda0");
  Par<double> dLambda    ("dLambda",    0.1, "\t\t# Step in searching for the lambda");
  Par<int>    followPeak ("followPeak",  -3, "\t# Wheather to determin zero frequency from the diverging pseudo-particle (-2: lamdba0==StartLambda, -1: Q==Q0, 0: follow b, 1: follow f, 2: foolow a)");
  Par<double> MissDopSt  ("MissDopSt", -100, "\t# Missing doping due to projection and finite mesh starting at MissDopSt");
}

int common::Na;
double common::U;
double common::T;
int common::baths;
function1D<double> common::Ed;
function1D<int> common::Ns;
function2D<int> common::Ms;
function1D<int> common::Mtot;
function1D<int> common::deg;
function2D<int> common::ncab;	  // index for hole diagrams
function2D<int> common::ncaf;	  // index for particle diagrams
function2D<double> common::prefactb; // prefactor for hole digrams
function2D<double> common::prefactf; // prefactor for particle diagrams
function2D<double> common::prefactG; // prefactor to calculate local Green's function
function1D<string> common::Eds;
function1D<double> common::nalpha;
function1D<double> common::miss_nd;
double common::beta;
double common::lambda0;
double common::Q;
double common::Q0;
double common::nd;
string common::outdir;
int common::totDeg;

int main (int argc, char *argv[], char *env[]){
  using namespace Parameters;
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  DynArray arr(80, &Sig, &Ac, &cix, &out, &gloc, &sig, &Ed, &U, &T, &Q0, &alpha, &max_diff, &max_steps,
	       &StartLambda, &EndLambda, &dLambda,  &followPeak, &prt, &MissDopSt, NULL);
  
  if (argc<2) {
    arr.printl(clog);
    return 0;
  }
  arr.ParsComLine(argc, argv);
  arr.print (clog);

  common::ParsInputFile(cix);
  
  common::SetParameters(Ed,U,T,Q0,out);
  RememberParams (argc, argv);
  Physical ph(common::Na, common::baths);
  Auxiliary au(common::Na, common::baths);
  
  if (!au.ReadSelfEnergy(Sig,Ed,T,U)) exit(1);
  if (!ph.ReadBathFunction(Ac)) exit(1);

  au.SetUpAverageAc(ph.omega(), ph.momega(), ph.Ac0(), ph.fe());

  if (followPeak==-3){
    StartLambda = StartLambda+au.minEnergy;
    EndLambda = EndLambda+au.minEnergy;
    common::Q0 = common::totDeg;
    followPeak=-1;
  }

  int n=0; double diff=1;
  while (diff>max_diff && n++<max_steps){
    clog<<"  "<<n<<".) KramarsKronig"<<endl;
    au.KramarsKronig();
    clog<<"  "<<n<<".) DeterminSpectralFunctions ";
    
    common::nd = au.DeterminSpectralFunctions(StartLambda,EndLambda,dLambda,followPeak);
    clog<<"  "<<n<<".) nd:  "<<common::nd<<endl; 
    au.PrintNorm(clog);
    if (prt) au.Print(n);

    clog<<"  "<<n<<".) SelfConsistentStep"<<endl;
    au.SetSignToZero();
    au.CalcSigmab();
    au.CalcSigmaf();
    
    clog<<"  "<<n<<".) Difference between steps: \t\t"<<COLOR(GREEN,(diff=au.Difference()))<<endl;
    au.DeterminSelfEnergies(Parameters::alpha);
  }
  ph.CalculateA00(au.omega(), au._Gp(), au._Gm());

  ph.DeterminG00(1.0);
  ph.KramarsKronig();
  ph.CalcSelfEnergy();

  ph.MissingDoping(MissDopSt);
  
  au.Print(0);
  ofstream outg((common::outdir+"/"+static_cast<string>(gloc)).c_str());
  print(outg, ph.omd, ph.G00, 25);
  ofstream outs((common::outdir+"/"+static_cast<string>(sig)).c_str());
  print(outs, ph.omd, ph.Sig, 25);
  return 0;
}
