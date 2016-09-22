#ifndef	__IFFS__
#define __IFFS__

#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
using namespace std;
using namespace Rcpp;

class IFFS {
  
public:
  IFFS(SEXP, SEXP, SEXP, SEXP, SEXP, int, int, float, float, float, float, float, float, float, float);
  
  
  /* Main methods */
  void ReadInput(SEXP);
  void ReadInput_SDB(SEXP);
  void ReadInput_ChipSet(SEXP);
  
  void Run();
  void RUN_IFFS(int);
  void SFS(int, set<int>&, set<int>&);
  double J(set<int>, int);
  double COD(set<int>, int);
  double COD_pen_rara(set<int>, int);
  double COD_pen_zero(set<int>, int);
  double IM(set<int>, int);
  double IM_pen_rara(set<int>, int);
  double IM_pen_zero(set<int>, int);
  double IM_Tsallis(set<int>, int);
  float SDB(set<int>, int);
  double TsallisEntropy(vector<vector<int> >, int, long double);
  double ConditionalTsallisEntropy(vector<vector<int> >, int, long double, long double, double);
  bool SBS(int, set<int>&, set<int>&);
  bool Replace(int, set<int>&, set<int>&);

  
  /* Auxiliary methods */
  void PrintSet(set<int>);
  void PrintMatrix(vector<vector<int> >);
	void PrintMatrix(float**);
  void salva_topX(vector<int> , double);
  void Write_top_X_preditores(int);
  void WriteOutput(set<int>);
  
private:
  vector<vector<int> > M;	// Time-series matrix of size genes x samples
  vector<double> score;	// score[i] Stores the currently max. score for a subset os size i
  float** M_SDB;	// STRING DB
  int** Chip_Set;	// Chip_Set matrix of size samples x 3
  NumericMatrix top_X_preditors;
  
  int ngenes;				// Number of genes, i.e., rows of M
  int nsamples;			// Number of samples, i.e., columns of M
  unsigned int max_k;					// The size of the final set os features
  int delta;
  float w_cod;
  float w_cod_rara;
  float w_cod_zero;
  float w_im;
  float w_im_rara;
  float w_im_zero;
  float w_tsallis;
  float w_sdb;
  std::ofstream f_resultado_com_score;
  std::string output;
  
};

#endif
