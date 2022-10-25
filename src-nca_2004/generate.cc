#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>

using namespace std;
typedef vector<int>::size_type vint;

// Crucial function that recursively generates the local base for any given number of baths ans its degeneracy
void generateBase(int b, vector<int>& state, vector<vector<int> >& base, int baths, const vector<int>& Ns) 
{ 
  for(int i = 0; i<=Ns[b]; i++){
    state[b] = i;
    if (b<baths-1) 
      generateBase(b+1, state,base,baths,Ns);
    else{      
      base.push_back(state);
    }
  }
}

// This comparisson is needed to sort the base by the total number of electrons in the
// local state. One could sort them also inside the subspace of the same number of particles.
// Not implementet - do not have a good idea what would be a good choice.
bool lessthan(const vector<int>& a, const vector<int>& b){
  int tota = 0; int totb = 0;
  for (vint i=0; i<a.size(); i++) tota += a[i];
  for (vint i=0; i<b.size(); i++) totb += b[i];
  return tota<totb;
}

// Just old good !
int Binomial(int n, int m)
{
  double Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  double r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return static_cast<int>(r/Mf);
}

// For parsing an entity in the parsInputChoice.
// Given a string like 10-12 return par of numbers
// 10 and 12.
pair<int,int> parsRange(const string& str){
  string::size_type pos;
  if ((pos=str.find('-'))<string::npos){
    int start = atoi(str.substr(0,pos).c_str());
    int end   = atoi(str.substr(pos+1).c_str());
    return make_pair(start,end);
  } else{
    int start = atoi(str.c_str());
    return make_pair(start,start);
  }
}

// Used to parse the input choice when the user cuts the base.
// Gets a string from the input (like 1-10,12,13) and returns
// a list of int numbers (in this 1,2,....10,12,13).
void parsInputChoice(string& str, list<int>& small, int base_size)
{
  list<pair<int,int> > keep;
  string::size_type pos;
  while((pos=str.find(','))<string::npos){
    keep.push_back(parsRange(str.substr(0,pos)));
    str.erase(0,pos+1);
  }
  
  keep.push_back(parsRange(str));
  for (list<pair<int,int> >::iterator i=keep.begin(); i!=keep.end(); i++)
    for (int j=(i->first); j<=(i->second); j++)
      if (j<base_size && j>=0) small.push_back(j);   
}

// Given the base and bath degeneracy calculates index to the base, total number of particles
// in the state, degeneracy of each state and NCA diagrams.
void SetUp(const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
	   vector<int>& Mtot, vector<int>& degeneracy){
  for(vint i=0; i<base.size(); i++){
    int nt=0; int dg=1;
    for (vint j=0; j<base[i].size(); j++){
      nt += base[i][j];
      dg *= Binomial(Ns[j],base[i][j]);
    }
    Mtot[i] = nt;
    degeneracy[i] = dg;
    index[base[i]]=i+1;
  }

}

void CalcNCADiag(int baths, const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
		 const vector<int>& degeneracy, vector<vector<int> >& ncab, vector<vector<int> >& ncaf)
{
  for (vint i=0; i<base.size(); i++){
    for (int b=0; b<baths; b++){
      vector<int> st = base[i];
      st[b]++;
      int indb = index[st]-1;
      ncab[i].push_back(indb);
      st[b]--; st[b]--;
      int indf = index[st]-1;
      ncaf[i].push_back(indf);
    }
  }
}

// Prints out the information calculated above
void Print(ostream& stream, int baths, const vector<vector<int> >& base, const vector<int>& Mtot, const vector<int>& degeneracy)
{
  stream<<"#"<<setw(2)<<"i"<<"  "<<setw(baths*2)<<"state"<<"   "<<setw(6)<<"Mtot"<<setw(6)<<"Deg"<<endl;
  for(vint i=0; i<base.size(); i++){
    stream<<setw(3)<<i<<"  ";
    for (vint j=0; j<base[i].size(); j++) stream<<setw(2)<<base[i][j];
    stream<<"   "<<setw(4)<<Mtot[i]<<setw(7)<<degeneracy[i]<<endl;
  }
}

void Print(ostream& stream, int baths, const vector<vector<int> >& base, const vector<int>& Ns, const vector<int>& Mtot,
	   const vector<int>& degeneracy, const vector<vector<int> >& ncab, const vector<vector<int> >& ncaf){
  stream<<"#"<<setw(2)<<"i"<<"  "<<setw(baths*2)<<"state"<<"  "<<setw(4)<<"Mtot"<<setw(7)<<"Deg"<<"  ";
  stream<<setw(5*baths)<<"ncab"<<setw(5*baths)<<"prefactb"<<setw(5*baths)<<"ncaf"<<setw(5*baths)<<"prefactf";
  stream<<setw(5*baths)<<"prefactG_loc"<<endl;
  for(vint i=0; i<base.size(); i++){
    stream<<setw(3)<<i<<"  ";
    for (int j=0; j<baths; j++) stream<<setw(2)<<base[i][j];
    stream<<"  "<<setw(4)<<Mtot[i]<<setw(7)<<degeneracy[i]<<"  ";
    for (int j=0; j<baths; j++) stream<<setw(5)<<ncab[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<Ns[j]-base[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<ncaf[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<base[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<(Ns[j]>0 ? degeneracy[i]*(Ns[j]-base[i][j])/Ns[j]: 0);
    stream<<endl;
  }
}
 
int main()
{
  int baths;
  cout<<" Give number of baths ... "<<endl;
  cin>>baths;
  if (!cin) {cerr<<"The value entered for number of baths is not valid\n";exit(1);}
    
  vector<int> Ns(baths);

  for (int i=0; i<baths; i++){
    cout<<" Give degeneracy of bath number "<<i<<" ... "<<endl;
    cin>>Ns[i];
    if (!cin) {cerr<<"The value entered for degeneracy of bath is not valid\n";exit(1);}
  }
  
  vector<vector<int> > base;
  vector<int> state(baths);
  generateBase(0,state,base,baths,Ns);
  sort(base.begin(),base.end(),lessthan);
  
  map<vector<int>,int> index;
  vector<int> Mtot(base.size());
  vector<int> degeneracy(base.size());
  
  SetUp(base, Ns, index, Mtot, degeneracy);
  Print(cout, baths, base, Mtot, degeneracy);
  
  cout<<" Printed are all possible states in this case. Which states would you like to keep?"<<endl;
  cout<<"  Enter your choice in format [xx-yy]+[,xx]+ (Example: 0-12,15,20)"<<endl;
  string str;
  cin>>str;

  list<int> small;
  parsInputChoice(str, small, base.size());
  
  vector<vector<int> > small_base(small.size());
  map<vector<int>,int> small_index;
  vector<int> small_Mtot(small_base.size());
  vector<int> small_degeneracy(small_base.size());
  vector<vector<int> > small_ncab(small_base.size());
  vector<vector<int> > small_ncaf(small_base.size());

  int l=0;
  for (list<int>::iterator i=small.begin(); i!=small.end(); i++,l++) small_base[l] = base[*i];
  
  SetUp(small_base, Ns, small_index, small_Mtot, small_degeneracy);
  CalcNCADiag(baths, small_base, Ns, small_index, small_degeneracy, small_ncab, small_ncaf);
  Print(cout, baths, small_base, Ns, small_Mtot, small_degeneracy, small_ncab,  small_ncaf);

  cout<<"Give output filename?"<<endl;
  string filename;
  cin>>filename;
  ofstream out(filename.c_str());
  out<<"# Input file for NCA impurity solver. Do not change this file if you are not absolutely aware of what you are doing!"<<endl;
  out<<baths<<" ";
  for (int i=0; i<baths; i++) out<<Ns[i]<<" ";
  out<<small_base.size()<<"  # Number of baths it's degeneracy an number of all local states"<<endl;
  Print(out, baths, small_base, Ns, small_Mtot, small_degeneracy, small_ncab,  small_ncaf);
  cout<<"Output written on "<<filename<<endl;
  return 0;
}
