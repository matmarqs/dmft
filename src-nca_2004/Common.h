#ifndef _COMMON_
#define _COMMON_
#include "zeroin.h"
#include "average.h"
#ifdef _STRSTREAM
#include <strstream>
#endif

using namespace std;

int Binomial(int n, int m)
{
  int Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  int r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return r/Mf;
}

//Common constants and variables
class common{
public:
  static double U;
  static double T;
  static int baths;
  static int Na;
  static function1D<int> Ns;
  static function2D<int> Ms;
  static function1D<int> Mtot;
  static function1D<int> deg;
  static function2D<int> ncab;	   // index for hole diagrams
  static function2D<int> ncaf;	   // index for particle diagrams
  static function2D<double> prefactb; // prefactor for hole digrams
  static function2D<double> prefactf; // prefactor for particle diagrams
  static function2D<double> prefactG; // prefactor to calculate local Green's function
  static function1D<double> Ed;
  static function1D<double> nalpha;
  static function1D<double> miss_nd;
  static double beta;
  static double delta;
  static double Q;
  static double Q0;
  static double nd;
  static double nd0;
  static double lambda0;
  static string outdir;
  static int totDeg;
  static function1D<string> Eds;
  static void SetParameters(Par<double>& Ed_, double U_, double T_, double Q0_, const string& outdir_)
  {
    Ed.resize(baths);
    int i=0;
    while (Ed_.next() && i<baths) {
      Ed[i] = Ed_;
      i++;
    }
    for (int j=i; j<baths; j++) Ed[j]=Ed[i-1];
    T = T_; 
    U = U_; 
    beta=1/T_;
    Q0 = Q0_;
    outdir = outdir_;
    Eds.resize(baths);
    for (int i=0; i<baths; i++){
      stringstream t; t<<"E"<<i;
      Eds[i] = t.str();
    }
    nalpha.resize(baths);
    miss_nd.resize(baths);
    for (int i=0; i<baths; i++) miss_nd[i]=0;
  }
  static void ParsInputFile(const string& filename);
  static void PrintParsedData(ostream& stream);
  static ostream& printHead(ostream& stream);
};

// Auxiliary self-energies and spectral functions
class Auxiliary{
  const int Na, baths;
  mesh1D om;
  function1D<double> fe;
  function1D<double> fedh;
  function1D<double> logo;
  function2D<double> Sigt;
  function2D<double> Sigtn;
  function2D<dcomplex> Sigc;
  function2D<double> Gt;
  function2D<double> Gp;
  function2D<double> Gm;
  vector<function2D<double> > aAc;
  function1D<double> Acx;
  AvFun<double> aF;
  function1D<double> Energy;
public:
  Auxiliary (int Na_, int baths_) : Na(Na_), baths(baths_), aAc(2*baths){};
  bool ReadSelfEnergy(const string& filename, const Par<double>& Ed, const Par<double>& T, const Par<double>& U);
  void KramarsKronig();
  double DeterminSpectralFunctions(double StartLambda, double EndLambda, double dLamdba, int followPeak);
  void PrintOutMeanQ(double StartLambda, double EndLambda);
  void PrintNorm(ostream& stream);
  void Print(int l, string dir);
  void Printn(int l);
  void SetSignToZero(){Sigtn=0.0;}
  
  void SetUpAverageAc(const mesh1D& omd, const mesh1D& momd, const function2D<double>& Ack, const function1D<double>& fed);
  void CalcSigmab();
  void CalcSigmaf();
  
  double Difference();
  void DeterminSelfEnergies(double alpha);
  const mesh1D& omega() const {return om;}
  double ferm(int i) const {return fe[i];}
  const function2D<double>& _Gp() const {return Gp;}
  const function2D<double>& _Gm() const {return Gm;}
  
  void PrintSign();
  
  double Q(double lambda);
  double operator()(double lambda);
  double minEnergy;
private:
  void Print_aAc(int l);
  void Print_Qux(int l);
  void Print_Sign(int l, int st, int en);
  void PrintOutMeanQ(int M, double StartLambda, double EndLambda);
};

// Physical electron spectral function and suscpetibility
// Physical observables
class Physical{
public:
  const int Na, baths;
  mesh1D omd;
  function2D<dcomplex> G00;
  function2D<double> A00;
  function2D<dcomplex> Sig;
private:
  mesh1D momd;
  function1D<double> fed;
  function1D<double> logod;
  function2D<double> Ac;
  vector<AvFun<double> > aF;
  function2D<double> Gtx;
  function2D<double> Cmp;
  function1D<double> tG;
  function1D<bool> Pexists;
public:
  Physical(int Na_, int baths_);
  bool ReadBathFunction(const string& filename);
  void CalculateA00(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm);
  void KramarsKronig();
  void DeterminG00(double alpha);
  double Difference();
  void Print(int l, string dir);
  
  const mesh1D& omega() const {return omd;}
  const mesh1D& momega() const {return momd;}
  const function1D<double>& fe() const {return fed;}
  const function2D<double>& Ac0() const {return Ac;}
  
  void PrintA00(ostream& out);
  void CalcSelfEnergy();
  void MissingDoping(double start);
private:
  //  void SetUpProducts(const mesh1D& om, const function2D<double>& Gp);
  void CalculateProducts(double u, double fu, const mesh1D& om, const function2D<double>& Gm);
  
  bool ReadBeginning(const string& filename, istream& input, int& n, int& m, bool& begincomment, double& center);
};

void AverageFunction(const mesh1D& omx, double u, const mesh1D& eps, AvFun<double>& aF, function<double>& aAc)
{
  apar ap;
  cintpar pi;
  tint position = omx.InitInterpLeft();
  InterpLeft(eps[0]-u, omx, position, pi);
  aF.InterpolateFirst(pi);
  InterpLeft(eps[1]-u, omx, position, pi);
  ap.SetUpCsFirst(u, eps);
  aAc[0] = aF.InterpolateNext(pi, ap) * eps.Dh(0);
  for (int j=1; j<eps.size()-1; j++){
    InterpLeft(eps[j+1]-u, omx, position, pi);
    ap.SetUpCs(u, j, eps, omx.Dh(pi.i));
    aAc[j] = aF.InterpolateNext(pi, ap) * eps.Dh(j);
  }
  ap.SetUpCsLast(u, eps);
  aAc[eps.size()-1] = aF.InterpolateLast(ap) * eps.Dh(eps.size()-1);
}

inline double product(const double* A, const double* G, int size)
{
  double sum = 0;
  for (int i=0; i<size; i++) sum += A[i]*G[i];
  return sum;
}

void Auxiliary::SetUpAverageAc(const mesh1D& omd, const mesh1D& momd, const function2D<double>& Ack, const function1D<double>& fed)
{
  int m = om.find_(0.0)+1;
  
  Acx.resize(omd.size());

  for (int b=0; b<baths; b++){
    aAc[b].resize(om.size(),om.size());

    for (int i=0; i<omd.size(); i++) Acx[i] = Ack[b][i]*(1-fed[i]);
    aF.SetUp(Acx,omd);
    for (int i=0; i<m; i++) AverageFunction(omd,om[i],om,aF,aAc[b][i]);

    for (int i=0; i<omd.size(); i++) Acx[i] = Ack[b][i]*fed[i];
    aF.SetUp(Acx,omd);
    for (int i=m; i<om.size(); i++) AverageFunction(omd,om[i],om,aF,aAc[b][i]);
  
    aAc[baths+b].resize(om.size(),om.size());
  
    for (int i=0; i<momd.size(); i++) Acx[momd.size()-i-1] = Ack[b][i]*fed[i];
    aF.SetUp(Acx,momd);
    for (int i=0; i<m; i++) AverageFunction(momd,om[i],om,aF,aAc[baths+b][i]); 

    for (int i=0; i<momd.size(); i++) Acx[momd.size()-i-1] = Ack[b][i]*(1-fed[i]);
    aF.SetUp(Acx,momd);
    for (int i=m; i<om.size(); i++) AverageFunction(momd,om[i],om,aF,aAc[baths+b][i]);
  }
}

void Auxiliary::CalcSigmab()
{
  int m = om.find_(0.0)+1;
  
  for (int j=0; j<Na; j++){
    for (int i=0; i<om.size(); i++){
      double sum=0;
      for (int b=0; b<baths; b++){
	int l = common::ncab[j][b];
	if (l<0) continue;
	double prf = common::prefactb[j][b];
	if (i<m)
	  sum += prf * product(aAc[b][i].MemPt(), Gm[l].MemPt(), om.size())/fe[i];
	else
	  sum += prf * product(aAc[b][i].MemPt(), Gp[l].MemPt(), om.size())/(1-fe[i]);
      }
      Sigtn[j][i] += sum;
    }
  }
}

void Auxiliary::CalcSigmaf()
{
  int m = om.find_(0.0)+1;
  
  for (int j=0; j<Na; j++){
    for (int i=0; i<om.size(); i++){
      double sum=0;
      for (int b=0; b<baths; b++){
	int l = common::ncaf[j][b];
	if (l<0) continue;
	double prf = common::prefactf[j][b];
	if (i<m)
	  sum += prf * product(aAc[baths+b][i].MemPt(), Gm[l].MemPt(), om.size())/fe[i];
	else
	  sum += prf * product(aAc[baths+b][i].MemPt(), Gp[l].MemPt(), om.size())/(1-fe[i]);
      }
      Sigtn[j][i] += sum;
    }
  }
}

inline ostream& common::printHead(ostream& stream)
{
  stream<<"# ";
  for (int i=0; i<baths; i++)stream<<" "<<Eds[i]<<"="<<Ed[i];
  stream<<" T="<<T<<" nd="<<nd<<" U="<<U<<" lambda0="<<lambda0;
  for (int i=0; i<baths; i++) stream<<" nb"<<i<<"="<<nalpha[i];
  double tmiss_nd=0;
  for (int i=0; i<baths; i++) {stream<<" md"<<i<<"="<<miss_nd[i]; tmiss_nd+=miss_nd[i];}
  stream<<" mdt="<<tmiss_nd;
  for (int i=0; i<baths; i++) stream<<" Ns"<<i<<"="<<Ns[i];
  return stream;
}

void RememberParams (int argc, char *argv[]){
  ofstream param ((common::outdir+"/history.nca").c_str(), ios::app);
  if (!param) cerr<<" Didn't suceeded to open params file!"<<(common::outdir+"/history.nca")<<endl;
  for (int i=0; i<argc; i++) param << argv[i] << " ";
  param << endl;
}

template <class T>
bool ReadValue(T& a, const std::string& variable, const std::string& str){
  std::string::size_type pos = str.find(variable);
  if (pos < std::string::npos){
    std::string::size_type poseq = str.find("=",pos);
    if (poseq<std::string::npos){
      std::istringstream streambuff(std::string(str,poseq+1));
      streambuff >> a;
    }
    return true;
  }
  return false;
}

bool Auxiliary::ReadSelfEnergy(const string& filename, const Par<double>& Ed, const Par<double>& T, const Par<double>& U){
  ifstream inputf(filename.c_str());
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);
  
  if (!input) {
    cerr << "Can't open input file: " << filename << endl;
    return false;
  }
  // Is the input file started with comment?
  bool begincomment = false;
  int n = 0;
  string str;
  const double SpecNumber = -100000; 
  double T_ = SpecNumber, U_ = SpecNumber;
  function1D<double> Ed_(baths);
  Ed_ = SpecNumber;
  double center = 0;
  getline(input,str);
  if (str.find('#')<string::npos){
    begincomment = true;
    for (int i=0; i<baths; i++) ReadValue(Ed_[i], common::Eds[i], str);
    ReadValue(T_, "T", str);
    ReadValue(U_, "U", str);
    if (!ReadValue(center, "peakposition", str)) center=0;
  } else n++;
  
  if (!Ed.IsSet() && Ed_[0]!=SpecNumber) for (int i=0; i<baths; i++) common::Ed[i] = Ed_[i];
  if (!T.IsSet()  && T_!=SpecNumber) common::T = T_;
  if (!U.IsSet() && U_!=SpecNumber) common::U = U_;
  common::beta = 1./common::T;

  Energy.resize(Na);
  minEnergy=0;
  // Calculates auxiliary Energies
  for (int i=0; i<Na; i++){
    Energy[i] = 0;
    for (int j=0; j<baths; j++)  Energy[i] += common::Ed[j]*common::Ms[i][j];
    Energy[i] += 0.5*common::Mtot[i]*(common::Mtot[i]-1)*common::U;
    if (Energy[i]<minEnergy) minEnergy = Energy[i];
  }

  clog<<"************* Parameters ****************"<<endl;
  clog<<"  		U  = "<<common::U<<endl;
  for (int i=0; i<baths; i++)
    clog<<"  		Ed"<<i<<" = "<<common::Ed[i]<<endl;
  clog<<"  		T  = "<<common::T<<endl;
  for (int i=0; i<baths; i++)
    clog<<"  		N"<<i<<" = "<<common::Ns[i]<<endl;
  for (int i=0; i<Na; i++){
    clog<<"	    state"<<i<<" = ";
    for (int j=0; j<baths; j++) clog<<setw(2)<<common::Ms[i][j];
    clog<<" with Energy"<<i<<"  = "<<Energy[i]<<endl;
  }
  clog<<"*****************************************"<<endl;
  
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for Sigm" << endl;
    return false;
  }
  getline(input,str);  n++;
#ifdef _STRSTREAM
  strstream oneline;
  oneline << str <<ends;
#else
  istringstream oneline(str);
#endif
  int m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;
  
  clog << filename << ": Number of entries: "<< n <<endl;
  clog << filename << ": Number of columns: "<< m <<endl;
  clog << filename << ": Peak-position "<< center <<endl;

  if (m<Na+1){
    cerr<<"ERROR: Not enough columns is input Sigma file. Exiting!"<<endl;
    return false;
  }
  inputf.seekg(0,ios::beg);
  clog<<"Premaknil na "<< inputf.tellg()<<endl;
  if (begincomment) inputf.ignore(1000,'\n');
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  
  om.resize(n);
  Sigt.resize(Na,n);
  
  int l=0;
  double omega;
  while (inputf>>omega && l<n){
    om[l] = omega;
    for (int i=0; i<Na; i++) {inputf>>Sigt(i,l); Sigt(i,l)*=-1;}
    getline(inputf, str);
    l++;
  }
  inputf.close();
  if (l<n) cerr<<"Something wrong by reading file "<<filename<<endl;
  om.SetUp(center);

  Sigc.resize(Na,om.size());
  Sigtn.resize(Na,om.size());
  Gt.resize(Na,om.size());
  Gp.resize(Na,om.size());
  Gm.resize(Na,om.size());
  fe.CalcFermOnMesh(common::beta, om);
  logo.CalcLogOnMesh(om);
  fedh.resize(om.size());
  for (int i=0; i<om.size(); i++) fedh[i] = fe[i]*om.Dh(i);
  return true;
}

void Auxiliary::KramarsKronig()
{
  for (int l=0; l<Na; l++){
    for (int i=0; i<om.size(); i++) Sigc(l,i).imag() = Sigt(l,i)*(1-fe[i]);
    Sigc[l].KramarsKronig(om, logo);
  }
}

double Lambda(double E, const function<dcomplex>& Sigc, const function<double>& Sigx, const mesh1D& om)
{
  // looking for lambda such that \widetilde{G} has maximum at zero frequency.
  // Sufficient condition is that the derivative of 1/\widetilde{G} is zero at zero frequency.
  // One gets a quadratic equation for lambda and thus two roots. Then one chooses the root that maximizes \widetilde{G}.
  // If no root exists, than we take lambda that minimizes linear coeficient in the expansion of 1/\widetilde{G}.
  // The latter equation is linear and one always gets unique solution.
  intpar p = om.Interp(0.0); int i=p.i;
  dcomplex cs = -E-Sigc(p);
  dcomplex ds = (Sigc[i+1]-Sigc[i])*om.Delta(i);
  double cr = cs.real();
  double ci = cs.imag();
  double dcr = 1-ds.real();
  double dci = -ds.imag();
  double dSigx   = (Sigx[i+1]-Sigx[i])*om.Delta(i);
  double x = Sigx[i]/dSigx;
  double determinant2 = x*(x*dcr*dcr+2*ci*dci)-ci*ci;
  // Minimum can not be at zero. Try to find lambda that minimizes the linear coefficient in the expansion of 1/G
  // If 1/G = a + b omega + c omega^2 +... and the below determinant is smaller than zero, coefficient b can not be
  // set to zero. Than return lambda that gives the smallest b.
  if (determinant2<=0) return dcr*x-cr;
  double d2 = -sqrt(determinant2);
  double d1 = -cr + dcr*x;
  double v1 = 1/(sqr(ci)+sqr(cr+d1+d2));
  double v2 = 1/(sqr(ci)+sqr(cr+d1-d2));
  cout<<"Lambda="<<d1+d2<<" "<<d1-d2<<" "<<v1<<" "<<v2<<endl;
  if (fabs(v1)>fabs(v2)) return d1+d2;
  else return d1-d2;
}

double Auxiliary::Q(double lambda)
{
  double sumQ=0;
  for (int j=0; j<Na; j++){
    double mune = -Energy[j]+lambda;
    double sum=0;
    for (int i=0; i<om.size(); i++)
      sum -= Sigt(j,i)*fedh[i]/(sqr(om[i]+mune-Sigc(j,i).real())+sqr(Sigc(j,i).imag()));
    sumQ += sum*common::deg[j];
  }
  return (sumQ/M_PI);
}

inline double Auxiliary::operator()(double lambda)
{
  double Q_ = Q(lambda);
  return Q_-common::Q0;
}

void Auxiliary::PrintOutMeanQ(double StartLambda, double EndLambda)
{
  double a0 = StartLambda;
  int M = 100;
  double da0 = (EndLambda-StartLambda)/M;
  cout.precision(16);
  for (int i=0; i<M; i++){
    cout << a0 << setw(25) << operator()(a0) << endl;
    a0 += da0;
  }
}

double Auxiliary::DeterminSpectralFunctions(double StartLambda, double EndLambda, double dLambda, int followPeak)
{
  double lambda0;
  if (followPeak>=0 && followPeak<Na)
    lambda0 = Lambda(Energy[followPeak], Sigc[followPeak], Sigt[followPeak], om);
  else if (followPeak==-2){
    lambda0 = minEnergy;
  }else{
    double a0 = StartLambda, b0 = 0;
    int sign=0, nn=0;
    while (!sign && nn++<100){
      double pQ = operator()(a0);
      while (!sign && a0<=b0) {
	double sQ = operator()(a0+dLambda);
	sign = pQ*sQ<0;
	pQ = sQ;
	if (!sign) a0 += dLambda;
      }
      if (!sign) dLambda /= 2.0;
    }
    
    if (nn>=100) {
      cerr << "Can't find root for <Q>" << endl;
      PrintOutMeanQ(StartLambda, EndLambda);
      exit(1);
    }
    
    // loking for zero (lambda0)
    lambda0 = zeroin(a0, a0+dLambda, *this, 1e-15*common::Q0);
  }

  common::lambda0 = lambda0;
  clog << setprecision(16) << "; lambda = "<<lambda0<<endl;

  double sumQ = 0, sumnd=0;
  function1D<double> dQ(Na);
  for (int j=0; j<Na; j++){
    double mune = -Energy[j]+lambda0;
    dQ[j]=0;
    for (int i=0; i<om.size(); i++){
      Gt(j,i) = Sigt(j,i)/(sqr(om[i]+mune-Sigc(j,i).real())+sqr(Sigc(j,i).imag()));
      Gm(j,i) = fe[i]*Gt(j,i);
      Gp(j,i) = (1-fe[i])*Gt(j,i);
      dQ[j] -= Gt(j,i)*fedh[i];
    }
    dQ[j] *= common::deg[j]/M_PI;
    sumQ += dQ[j];
    sumnd += dQ[j]*common::Mtot[j];
  }
  clog<<"       Q = "<<sumQ<<endl;
  for (int j=0; j<Na; j++)
    clog<<setprecision(16)<<"       n"<<j<<"="<<dQ[j]/sumQ<<endl;

  for (int b=0; b<baths; b++){
    common::nalpha[b]=0;
    for (int j=0; j<Na; j++) common::nalpha[b] += dQ[j]*common::Ms[j][b];
    common::nalpha[b]/=sumQ;
  }
  common::Q = sumQ;
  
  //  if (fabs(sumQ-common::Q0)>1e-10) cerr<<"Something wrong with Q "<<sumQ<<"!"<<endl;
  return sumnd/sumQ;
}

void Auxiliary::Print(int l, string dir="")
{
  string filename;
  if (l<0) filename = common::outdir+"/Sigma"+dir;
  else filename = NameOfFile(common::outdir+"/Sigma", l);
  ofstream out1(filename.c_str());  out1.precision(16);
  common::printHead(out1)<<" peakposition="<<om.dcenter()<<endl;
  for (int i=0; i<om.size(); i++){
    out1<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out1<<setw(25)<<-Sigt(j,i);
    out1<<endl;
  }
  if (l<0) filename = common::outdir+"/Spec"+dir;
  else filename = NameOfFile(common::outdir+dir+"/Spec", l);
  ofstream out2(filename.c_str());  out2.precision(16);
  common::printHead(out2)<<" peakposition="<<om.dcenter()<<endl;
  for (int i=0; i<om.size(); i++){
    out2<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out2<<setw(25)<<-Gt(j,i);
    out2<<endl;
  }
}

void Auxiliary::Printn(int l)
{
  string filename;
  filename = NameOfFile(common::outdir+"/nSigma", l);
  ofstream out1(filename.c_str());  out1.precision(16);
  common::printHead(out1)<<" peakposition="<<om.dcenter()<<endl;
  for (int i=0; i<om.size(); i++){
    out1<<setw(25)<<om[i];
    for (int j=0; j<Na; j++) out1<<setw(25)<<-Sigtn(j,i);
    out1<<endl;
  }
}

Physical::Physical(int Na_, int baths_) : Na(Na_), baths(baths_), aF(Na)
{
  Pexists.resize(Na);
  for (int i=0; i<Na; i++){
    Pexists[i]=false;
    for (int b=0; b<baths; b++) if (common::ncab[i][b]>0) {Pexists[i]=true; break;}
  }
}

bool Physical::ReadBeginning(const string& filename, istream& input, int& n, int& m, bool& begincomment, double& center)
{
  if (!input) {
    cerr << "Can't open input file: " << filename << endl;
    return false;
  }
  // Is the input file started with comment?
  begincomment = false;
  n = 0;
  string str;
  getline(input,str);
  if (str.find('#')<string::npos){
    begincomment = true;
    if (!ReadValue(center, "peakposition", str)) center=0;
  } else n++;
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for Sigm" << endl;
    return false;
  }
  getline(input,str);  n++;
  stringstream oneline;
  oneline << str << ends;
  m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;

  clog << filename << ": Number of entries: "<< n <<endl;
  clog << filename << ": Number of columns: "<< m <<endl;
  clog << filename << ": Peak-position "<< center <<endl;
  
  input.seekg(0, ios::beg);
  input.clear();
  if (begincomment) getline(input, str);
  
  return true;
}

bool Physical::ReadBathFunction(const string& filename)
{
  ifstream inputf(filename.c_str());
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);

  if (!input) {
    cerr << "Can't open input file: " << filename << endl;
    return false;
  }
  // Is the input file started with comment?
  bool begincomment = false;
  int n = 0;
  string str;
  double center=0;
  getline(input,str);
  if (str.find('#')<string::npos){
    begincomment = true;
    if (!ReadValue(center, "peakposition", str)) center=0;
  } else n++;
  // Computes the number of columns in file
  if (!input) {
    cerr << "ERROR: Wrong file format for Ac" << endl;
    return false;
  }
  getline(input,str);  n++;
#ifdef _STRSTREAM
  strstream oneline;
  oneline << str <<ends;
#else
  istringstream oneline(str);
#endif
  int m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;

  clog << filename << ": Number of entries: "<< n <<endl;
  clog << filename << ": Number of columns: "<< m <<endl;
  clog << filename << ": Peak-position "<< center <<endl;

  if (m<baths+1){
    cerr<<"ERROR: Not enough columns in bath input file! Exiting..."<<endl;
    return false;
  }
  inputf.seekg(0, ios::beg);
  clog<<"Premaknil na "<< inputf.tellg()<<endl;
  if (begincomment) inputf.ignore(1000,'\n');
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  
  omd.resize(n);
  momd.resize(n);
  G00.resize(baths,n);
  A00.resize(baths,n);
  Sig.resize(baths,n);
  Ac.resize(baths,n);
  
  int l=0;
  double omega;
  while (inputf>>omega && l<n){
    omd[l] = omega;
    for (int j=0; j<baths; j++) inputf>>Ac(j,l);
    getline(inputf, str);
    momd[n-l-1] = -omd[l];
    l++;
  }
  inputf.close();
  if (l<n) cerr<<"Something wrong by reading file "<<filename<<endl;
  omd.SetUp(center);
  momd.SetUp(-center);

  fed.CalcFermOnMesh(common::beta, omd);
  logod.CalcLogOnMesh(omd);
  
  return true;
}

void Physical::CalculateProducts(double u, double fu, const mesh1D& om, const function2D<double>& Gm)
{
  apar ap;
  cintpar pi;
  tint position = om.InitInterpLeft();
  InterpLeft(om[0]-u, om, position, pi);
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].InterpolateFirst(pi);
  InterpLeft(om[1]-u, om, position, pi);
  ap.SetUpCsFirst(u, om);
  for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,0) = aF[i].InterpolateNext(pi, ap) * om.Dh(0);
  for (int j=1; j<om.size()-1; j++){
    InterpLeft(om[j+1]-u, om, position, pi);
    ap.SetUpCs(u, j, om, om.Dh(pi.i+1));
    for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,j) = aF[i].InterpolateNext(pi, ap) * om.Dh(j);
  }
  ap.SetUpCsLast(u, om);
  for (int i=0; i<Na; i++) if (Pexists[i]) Gtx(i,om.size()-1) = aF[i].InterpolateLast(ap) * om.Dh(om.size()-1);

  Cmp.resize(Na,Na);
  for (int i=0; i<Na; i++){
    for (int b=0; b<baths; b++){
      int l=common::ncab(i,b);
      if (l>=0)
	Cmp(i,l) = product(Gtx[i].MemPt(),Gm[l].MemPt(),om.size())/fu;
    }
  }
}

void Physical::CalculateA00(const mesh1D& omega, const function2D<double>& Gp, const function2D<double>& Gm)
{
  int m = omd.find_(0.0)+1;

  Gtx.resize(Na, omega.size());
  
  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].SetUp(Gp[i],omega);
  for (int i=0; i<m; i++){
    CalculateProducts(omd[i], fed[i], omega, Gm);
    for (int b=0; b<baths; b++){
      double sum=0;
      for (int j=0; j<Na; j++)
	if (common::ncab[j][b]>=0)
	  sum += (common::prefactG[j][b]*Cmp(j,common::ncab[j][b]))/(M_PI*M_PI)/common::Q;
      A00(b,i) = sum;
    }
  }

  for (int i=0; i<Na; i++) if (Pexists[i]) aF[i].SetUp(Gm[i],omega);
  for (int i=m; i<omd.size(); i++){
    CalculateProducts(omd[i], (1-fed[i]), omega, Gp);
    for (int b=0; b<baths; b++){
      double sum=0;
      for (int j=0; j<Na; j++)
	if (common::ncab[j][b]>=0)
	  sum += (common::prefactG[j][b]*Cmp(j,common::ncab[j][b]))/(M_PI*M_PI)/common::Q;
      A00(b,i) = sum;
    }
  }
}

//  void Physical::CalcExternHilbert(const string& executable, const string& DOSfile)
//  {
//    {
//      ofstream A0temp("/tmp/A0temp"); A0temp.precision(16);
//      print(A0temp, omd, A00, 25);
//    }
//    stringstream command;
//    command<<executable<<" -i /tmp/A0temp "<<DOSfile<<" > /tmp/Actemp"<<ends;
//    int ret = system(command.str().c_str());
//    if (ret) cerr<<"system command returned "<<ret<<endl;
  
//    int l=0, mdata=7;
//    vector<double> data(mdata);
//    ifstream Actemp("/tmp/Actemp");
//    while(Actemp && l<omd.size()){
//      for (int j=0; j<mdata; j++) Actemp>>data[j];
//      Actemp.ignore(500, '\n');
//      if (fabs(omd[l]-data[0])>1e-5) cerr<<"Mesh has changed after external Hilbert!"<<endl;
//      Ac[l] = data[2];
//      Sig[l] = dcomplex(data[3],data[4]);
//      G00[l] = dcomplex(data[5],data[6]);
//      l++;
//    }
//  }

inline void Physical::KramarsKronig()
{
  for (int b=0; b<baths; b++) G00[b].KramarsKronig(omd,logod);
}

void Physical::CalcSelfEnergy()
{
  for (int b=0; b<baths; b++){
    for (int i=0; i<omd.size(); i++){
      double Deltar = ::KramarsKronig(Ac[b], omd, omd[i], i, Ac[b][i]);
      dcomplex Delta(-M_PI*Deltar,-M_PI*Ac[b][i]);
      Sig[b][i] = omd[i]-common::Ed[b]-Delta-1/G00[b][i];
      if (Sig[b][i].imag()>0) Sig[b][i].imag()=0.0;
    }
  }
}

void Physical::Print(int n, string dir="")
{
  string filename;
  if (n<0) filename = common::outdir+"/A00"+dir;
  else filename = common::outdir+NameOfFile("/A00",n,3);
  ofstream out(filename.c_str()); out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out <<setw(25)<<omd[i];
    for (int b=0; b<baths; b++)
      out<<setw(25)<<A00[b][i]<<setw(25)<<G00[b][i]<<setw(25)<<-Sig[b][i];
    out<<endl;
  }
}

void Auxiliary::DeterminSelfEnergies(double alpha){
  double beta=1-alpha;
  for (int j=0; j<Na; j++)
    for (int i=0; i<om.size(); i++)
      Sigt(j,i) = beta*Sigt(j,i)+alpha*Sigtn(j,i);
}

void Physical::DeterminG00(double alpha)
{
  double beta=1-alpha;
  double alphapi=-alpha*M_PI;
  for (int b=0; b<baths; b++){
    for (int j=0; j<omd.size(); j++)
      G00[b][j].imag()=beta*G00[b][j].imag()+alphapi*A00[b][j];
  }
}

void Auxiliary::PrintNorm(ostream& stream)
{
  stream<<"    Norm of Spectral functions: "<<endl<<"   ";
  stream.setf(ios::fixed);
  for (int i=0; i<Na; i++){
    double sum=0;
    for (int j=0; j<om.size(); j++)
      sum += Gp(i,j)*om.Dh(j);

    sum/=-M_PI;    
    double norm0=1;
    stream<<setprecision(4)<<" ";
    
    if (fabs(sum-norm0)<1e-2)
      stream<<COLOR(GREEN,setw(2)<<i<<":"<<setw(8)<<sum)<<" ";
    else if (fabs(sum-norm0)<1e-1)
      stream<<COLOR(YELLOW,setw(2)<<i<<":"<<setw(8)<<sum)<<" ";
    else 
      stream<<COLOR(PURPLE,setw(2)<<i<<":"<<setw(8)<<sum)<<" ";
    if ((i+1)%6==0) stream<<endl<<"   ";
  }
  stream<<endl;
  for (int b=0; b<baths; b++){
    stream<<setprecision(4)<<" "<<COLOR(BLUE,setw(2)<<b<<":"<<setw(8)<<common::nalpha[b])<<" ";
  }
  stream<<endl;
  stream.unsetf(ios::fixed);
}

void Physical::PrintA00(ostream& out)
{
  out.precision(16);
  common::printHead(out)<<" peakposition=" << omd.dcenter()<<endl;
  for (int i=0; i<omd.size(); i++){
    out<<setw(25)<<omd[i];
    for (int b=0; b<baths; b++)
      out<<setw(25)<<A00[i];
    out<<endl;
  }
}

//  void Physical::ReadA00(istream& inp)
//  {
//    inp.ignore(200,'\n');
//    double omega;
//    for (int i=0; i<omd.size(); i++){
//      inp>>omega;
//      if (fabs(omega-omd[i])>1e-5) cerr<<"Mesh has changed when the spectral function was corrected!"<<endl;
//      if (!inp) cerr<<"Something wrong in reading changed spectral function!"<<endl;
//      inp>>A00[i];
//      inp.ignore(200,'\n');
//    }
//  }

double Auxiliary::Difference(){
  double diff=0, norm=0;
  for (int j=0; j<Na; j++){
    for (int i=0; i<om.size(); i++){
      diff += fabs(Sigtn(j,i)-Sigt(j,i));
      norm += 0.5*fabs(Sigtn(j,i)+Sigtn(j,i));
    }
  }
  return diff/norm;
}

//  double Physical::Difference(){
//    double diff=0, norm=0;
//    for (int j=0; j<omd.size(); j++){
//      diff += fabs(G00[j].imag()+M_PI*A00[j]);
//      norm += 0.5*fabs(G00[j].imag()-M_PI*A00[j]);
//    }
//    return diff/norm;
//  }

/******************* Used only for debugging **********************/
void Auxiliary::PrintSign()
{
  for (int i=0; i<Na; i++){
    ofstream out(NameOfFile("Sign",i,2).c_str());
    out.precision(16);
    for (int j=0; j<om.size(); j++)
      out<<setw(25)<<om[j]<<setw(25)<<-Sigtn[i][j]<<endl;
  }
}

void Auxiliary::Print_aAc(int l)
{
  for (int i=0; i<aAc[0].size_N(); i++){
    ofstream out(NameOfFile_("aAc",l,i,1,3).c_str());
    out.precision(16);
    for (int j=0; j<aAc[0].size_Nd(); j++){
      out<<setw(25)<<om[j]<<setw(25)<<aAc[0][i][j]/om.Dh(j)<<endl;
    }
  }
}

/******************* New things ******************************/
void common::ParsInputFile(const string& filename)
{
  ifstream input(filename.c_str());
  string line;
  getline(input,line);
  input>>baths;
  Ns.resize(baths);
  for (int i=0; i<baths; i++) input>>Ns[i];
  input>>Na;
  getline(input,line); getline(input,line);
  if (!input){ cerr<<filename<<" file not recognized. Error in first 3 lines!"<<endl; exit(1);}
  deg.resize(Na);
  Ms.resize(Na,baths);
  Mtot.resize(Na);
  ncab.resize(Na, baths);
  ncaf.resize(Na, baths);
  prefactb.resize(Na, baths);
  prefactf.resize(Na, baths);
  prefactG.resize(Na, baths);
  for (int i=0; i<Na; i++){
    getline(input, line);
    if (!input){ cerr<<filename<<" file not recognized. Error in line number "<<i+3<<endl; exit(1);}
    stringstream thisline(line);
    int lc;
    thisline>>lc;
    for (int j=0; j<baths; j++) thisline>>Ms[i][j];
    thisline>>Mtot[i]>>deg[i];
    for (int j=0; j<baths; j++) thisline>>ncab[i][j];
    for (int j=0; j<baths; j++) thisline>>prefactb[i][j];
    for (int j=0; j<baths; j++) thisline>>ncaf[i][j];
    for (int j=0; j<baths; j++) thisline>>prefactf[i][j];
    for (int j=0; j<baths; j++) thisline>>prefactG[i][j];
    if (!input){ cerr<<filename<<" file not recognized. Error in line number "<<i+3<<endl; exit(1);}
  }
  PrintParsedData(cout);
  totDeg = 0;
  for (int i=0; i<Na; i++) totDeg += deg[i];
}

void common::PrintParsedData(ostream& stream)
{
  stream<<baths<<" ";
  for (int i=0; i<baths; i++) stream<<Ns[i]<<" ";
  stream<<Na<<endl;
  for (int i=0; i<Na; i++){
    for (int j=0; j<baths; j++) stream<<setw(2)<<Ms[i][j];
    stream<<setw(4)<<Mtot[i]<<setw(4)<<deg[i];
    for (int j=0; j<baths; j++) stream<<setw(4)<<ncab[i][j];
    for (int j=0; j<baths; j++) stream<<setw(4)<<prefactb[i][j];
    for (int j=0; j<baths; j++) stream<<setw(4)<<ncaf[i][j];
    for (int j=0; j<baths; j++) stream<<setw(4)<<prefactf[i][j];
    for (int j=0; j<baths; j++) stream<<setw(4)<<prefactG[i][j];
    stream<<endl;
  }
}

void print(std::ostream& stream, const mesh1D& om, const function2D<dcomplex>& f, int width=20)
{
  if (om.size()!=f.size_Nd()) std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++){
    stream <<std::setw(width)<<om[i];
    for (int j=0; j<f.size_N(); j++) stream<<std::setw(width)<<f(j,i);
    stream<<std::endl;
  }
}

void Physical::MissingDoping(double start)
{
  cout<<"Missing doping : ";
  for (int b=0; b<baths; b++){
    double sum = 0;
    for (int i=0; i<omd.size(); i++) {
      if (omd[i]>start) sum += G00[b][i].imag()*fed[i]*omd.Dh(i);
    }
    sum *= -common::Ns[b]/M_PI;
    common::miss_nd[b] = common::nalpha[b]-sum;
    cout<<b<<" : "<<common::miss_nd[b]<<" ";
  }
  cout<<endl;
}

#endif
