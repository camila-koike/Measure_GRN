/*
IFFS - Improved Forward Floating Selection

Carlos Higa
January, 13, 2016.
*/


#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <Rcpp.h>
#include "IFFS.hpp"




using namespace Rcpp;
using namespace std;


IFFS::IFFS(SEXP InputExpr, SEXP InputSDB, SEXP InputCS, SEXP output, SEXP resultado_com_score, int max_k, int delta,
           float w_cod, float w_cod_rara, float w_cod_zero, float w_im, float w_im_rara, float w_im_zero, float w_tsallis, float w_sdb)
{
  
  this->max_k = max_k;
  this->delta = delta;
  this->w_cod = w_cod;
  this->w_cod_rara = w_cod_rara;
  this->w_cod_zero = w_cod_zero;
  this->w_im = w_im;
  this->w_im_rara = w_im_rara;
  this->w_im_zero = w_im_zero;
  this->w_tsallis = w_tsallis;
  this->w_sdb = w_sdb;
  this->top_X_preditors = NumericMatrix(10,max_k+1);
  this->output = Rcpp::as<std::string>(output);

  std::string fname = Rcpp::as<std::string>(resultado_com_score);
  f_resultado_com_score.open(fname.c_str(),std::ios::out);
  if (!f_resultado_com_score.is_open()) {
    Rcpp::Rcout << "Could not open " << fname.c_str() << endl;
    exit(1);
  }

  
  
  this->ReadInput(InputExpr);
  this->ReadInput_SDB(InputSDB);
  this->ReadInput_ChipSet(InputCS);
  this->Run();
  
}


void IFFS::ReadInput(SEXP InputExpr)
{
	ifstream fin;
	string row;
	vector<int> v;
	string::iterator it;

  fin.open(CHAR(STRING_ELT(InputExpr,0)));
	if (!fin.is_open()) {
		cout << "Could not open file " << (CHAR(STRING_ELT(InputExpr,0))) << "." << endl;
		exit(1);
	}

	/* Reads the time series and stores it in a matrix M nxm, where n = number of genes and m = number of samples. */
	while (getline(fin, row)) {
		it = std::remove(row.begin(), row.end(), ' ');
		row.erase(it, row.end());
		v.clear();
		for (it = row.begin(); it != row.end(); ++it) {
			v.push_back(static_cast<int>(*it) - '0');
		}
		this->M.push_back(v);
	}

	

	// Number of genes and samples, respectively.
	this->ngenes = this->M.size();
	this->nsamples = this->M[0].size();


	fin.close();


}

void IFFS::ReadInput_SDB(SEXP InputSDB)
{
  ifstream fin;
  string row;
  vector<double> v;
  string::iterator it;
  
  M_SDB = new float*[this->ngenes];
  for (int i = 0; i < this->ngenes; i++)
    M_SDB[i] = new float[this->ngenes];
  
  for (int i = 0; i < this->ngenes; i++)
    for (int j = 0; j < this->ngenes; j++)
      M_SDB[i][j] = 0;
  
  fin.open((CHAR(STRING_ELT(InputSDB,0))));
  if (!fin.is_open()) {
    cout << "Could not open file " << CHAR(STRING_ELT(InputSDB,0)) << "." << endl;
    exit(1);
  }
  int linha, coluna;
  float score;
  
  
  while (!fin.eof())
  {
    fin >>linha;
    fin >>coluna;
    fin >>score;
    M_SDB[linha][coluna] = score;
    M_SDB[coluna][linha] = score;
  }
  

  fin.close();  

  
}

void IFFS::ReadInput_ChipSet(SEXP InputSDB)
{
  ifstream fin;
  Chip_Set = new int*[this->nsamples];
  for (int i = 0; i < this->nsamples; i++)
    Chip_Set[i] = new int[2];
  
  for (int i = 0; i < this->nsamples; i++)
    for (int j = 0; j < 2; j++)
      Chip_Set[i][j] = -1;
  
  fin.open((CHAR(STRING_ELT(InputSDB,0))));
  if (!fin.is_open()) {
    cout << "Could not open file " << CHAR(STRING_ELT(InputSDB,0)) << "." << endl;
    exit(1);
  }
  int c1, c2;
 
 for (int i = 0; i < this->nsamples; i++)
 {
    fin >>c1;
    fin >>c2;
    Chip_Set[i][0] = c1;
    Chip_Set[i][1] = c2;
    //cout << Chip_Set[i][0] << " " << Chip_Set[i][1] << endl;
  }
  
  fin.close();  

}

void IFFS::Run()
{
	int target;

	for (target = 0; target < this->ngenes; ++target) {
		this->RUN_IFFS (target);
	}
	f_resultado_com_score.close();
}

void IFFS::RUN_IFFS(int target)
{
	set<int> All_Candidatos, Sub_Candidatos;
	int i;

	// All_Candidatos is the set of all features, including the target
	for (i = 0; i < this->ngenes; ++i) {
		All_Candidatos.insert(i);
	}

	// Sub_Candidatos is the current set of selected features. Starts empty.
	Sub_Candidatos.clear();

	// score[i] stores the current max score for |Sub_Candidatos| = i. Position 0 is not used.
	this->score.clear();
	this->score.assign(this->max_k+1, 0.0);
	
	// IFFS algorithm
	while (true) {
		SFS(target, All_Candidatos, Sub_Candidatos);

		if (Sub_Candidatos.size() == this->max_k) {
			break;
		}

		while (SBS(target, All_Candidatos, Sub_Candidatos)) {

		}

		while (Replace(target, All_Candidatos, Sub_Candidatos)) {
			while (SBS(target, All_Candidatos, Sub_Candidatos)) {

			}
		}

	}


	//cout << "Final result for target " << target << " ";
	//PrintSet(Sub_Candidatos);
	WriteOutput(Sub_Candidatos);
	//cout <<  this->score[Sub_Candidatos.size()] << endl ;


}

/* Sequential Forward Selection */
void IFFS::SFS(int target, set<int>& All_Candidatos, set<int>& Sub_Candidatos)
{
	set<int>::iterator it;
	double cur_score, max_score;
	int best_feat;

	max_score = -1.0;
	for (it = All_Candidatos.begin(); it != All_Candidatos.end(); ++it) {
		Sub_Candidatos.insert(*it);
		cur_score = this->J(Sub_Candidatos, target);
		if (cur_score > max_score) {
			max_score = cur_score;
			best_feat = *it;
		}
		Sub_Candidatos.erase(*it);
	}

	All_Candidatos.erase(best_feat);
	Sub_Candidatos.insert(best_feat);
	this->score[Sub_Candidatos.size()] = max_score;
	//cout << "SFS: insert " << best_feat << " Score: " << max_score << endl;
	//PrintSet(All_Candidatos);
}

double IFFS::J(set<int> S, int target)
{
  double ret = 0;
  if(this->w_cod != 0)
    ret += (w_cod*COD(S,target));
  if(this->w_cod_rara != 0)
    ret += (w_cod_rara*COD_pen_rara(S,target));
  if(this->w_cod_zero != 0)
    ret += (w_cod_zero*COD_pen_zero(S,target));
  if(this->w_im != 0)
    ret += (w_im*IM(S,target));
  if(this->w_im_rara != 0)
    ret += (w_im_rara*IM_pen_rara(S,target));
  if(this->w_im_zero != 0)
    ret += (w_im_zero*IM_pen_zero(S,target));
  if(this->w_tsallis != 0)
    ret += (w_tsallis*IM_Tsallis(S,target));
  if(this->w_sdb != 0)
    ret += (w_sdb*SDB(S,target));
  
  return ret;
  
}

/* Compute the temporal CoD */
double IFFS::COD(set<int> Sub_Candidatos, int target)
{
	int num_comb, i, t, bin, tgt_one, tgt_zero;
	vector<vector<int> > table;
	vector<int> v;
	set<int>::iterator it;
	double epsilon_Z, epsilon_target;
	double samples = 0.0;

	/* Initializing the conditional probability table. */
	num_comb = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
	v.assign(2, 0);
	for (i = 0; i < num_comb; ++i) {
		table.push_back(v);
	}
	
	/* nsamples-1 because we are considering time series data. */
	for (t = 0; t < this->nsamples-1; ++t) {
		i = 0;
		bin = 1;
		for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
			i += this->M[*it][t] * bin;
			bin *= 2;
		}
		if(this->Chip_Set[t][1] == -1) //estático
		{
		  samples += 1.0;
		  if(this->M[target][t] == 1)
		    table[i][1]++;
		  else
		    table[i][0]++;
		}
		else
		{
		  if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
		  {
		    if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
		    {
		      samples += 1.0;
		      if (this->M[target][t+1] == 1) {
		        table[i][1]++;
		      } else {
		        table[i][0]++;
		      }
		    }
		  }
		}
	}

	// Compute epsilon_Z
	epsilon_Z = 0.0;
	tgt_one = 0;
	for (i = 0; i < num_comb; ++i) {
		epsilon_Z += (table[i][0] > table[i][1])? static_cast<double>(table[i][1]) : static_cast<double>(table[i][0]);

		/* To compute the epsilon_target */
		tgt_one += table[i][1];
	}

	epsilon_Z /= static_cast<double>(samples);

	// Compute epsilon_target
	tgt_zero = samples - tgt_one;
	epsilon_target = (tgt_one > tgt_zero)? static_cast<double>(tgt_zero) : static_cast<double>(tgt_one);
	epsilon_target /= static_cast<double>(samples);


	if (epsilon_target == 0.0 && epsilon_Z == 0.0) {
		return 1.0;
	} else if ((epsilon_target == 0.0 && epsilon_target < epsilon_Z) || (epsilon_target > 0.0 && epsilon_target < epsilon_Z)) {
		return 0.0;
	} else {
		return (epsilon_target - epsilon_Z)/epsilon_target;
	}
}

/* Compute the temporal CoD */
double IFFS::COD_pen_rara(set<int> Sub_Candidatos, int target)
{
  int num_comb, i, t, bin, tgt_one, tgt_zero;
  vector<vector<int> > table;
  vector<int> v;
  set<int>::iterator it;
  double epsilon_Z, epsilon_target;
  double samples = 0.0;
  
  /* Initializing the conditional probability table. */
  num_comb = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
  v.assign(2, 0);
  for (i = 0; i < num_comb; ++i) {
    table.push_back(v);
  }

  /* nsamples-1 because we are considering time series data. */
  for (t = 0; t < this->nsamples-1; ++t) {
    i = 0;
    bin = 1;
    for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
      i += this->M[*it][t] * bin;
      bin *= 2;
    }
    if(this->Chip_Set[t][1] == -1) //estático
    {
      samples += 1.0;
      if(this->M[target][t] == 1)
        table[i][1]++;
      else
        table[i][0]++;
    }
    else
    {
      if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
      {
        if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
        {
          samples += 1.0;
          if (this->M[target][t+1] == 1) {
            table[i][1]++;
          } else {
            table[i][0]++;
          }
        }
      }
    }
  }


  //Rcout << "!!!!!!!!!!!!COM PENALIZAÇÃO RARAMENTE OBSERVADA!!!!!!!"  << endl;
  double sum_Ey_linha;
  double alpha, betha;
  int M, N,s;
  alpha = 1;
  s = samples;
  M = num_comb;
  N = 0;
  betha = 0.8;

  tgt_one = 0;
  tgt_zero = 0;
  sum_Ey_linha  = 0;
  
  
  N = 0;
  
  for(int p = 0; p < num_comb; p++)
  {
    tgt_one += table[i][1];
    tgt_zero += table[i][1];
    int sum_iesima_comb = (table[p][0] + table[p][1]);
    if(sum_iesima_comb > 1)
    {
      N++;
      if(table[p][1] > table[p][0])
      {
        sum_Ey_linha += ((double)table[p][0]/(double)sum_iesima_comb)*(double)(sum_iesima_comb/samples);
      }
      else
      {
        sum_Ey_linha += ((double)table[p][1]/(double)sum_iesima_comb)*(double)(sum_iesima_comb/samples);
      }
    }
  }
  
  // Compute epsilon_target
  epsilon_target = (tgt_one > tgt_zero)? static_cast<double>(tgt_zero) : static_cast<double>(tgt_one);
  epsilon_target /= static_cast<double>(samples);
  
  // Compute epsilon_Z
  epsilon_Z = 0.0;
  epsilon_Z = (double)((((M-N)/s)*(1-betha))+(double)sum_Ey_linha);

  if (epsilon_target == 0.0 && epsilon_Z == 0.0) {
    return 1.0;
  } else if ((epsilon_target == 0.0 && epsilon_target < epsilon_Z) || (epsilon_target > 0.0 && epsilon_target < epsilon_Z)) {
    return 0.0;
  } else {
    return (epsilon_target - epsilon_Z)/epsilon_target;
  }

  
}

/* Compute the temporal CoD */
double IFFS::COD_pen_zero(set<int> Sub_Candidatos, int target)
{
  int num_comb, i, t, bin, tgt_one, tgt_zero;
  vector<vector<int> > table;
  vector<int> v;
  set<int>::iterator it;
  double epsilon_Z, epsilon_target;
  double samples = 0.0;
  
  /* Initializing the conditional probability table. */
  num_comb = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
  v.assign(2, 0);
  for (i = 0; i < num_comb; ++i) {
    table.push_back(v);
  }
  
  /* nsamples-1 because we are considering time series data. */
  for (t = 0; t < this->nsamples-1; ++t) {
    i = 0;
    bin = 1;
    for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
      i += this->M[*it][t] * bin;
      bin *= 2;
    }
    if(this->Chip_Set[t][1] == -1) //estático
    {
      samples += 1.0;
      if(this->M[target][t] == 1)
        table[i][1]++;
      else
        table[i][0]++;
    }
    else
    {
      if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
      {
        if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
        {
          samples += 1.0;
          if (this->M[target][t+1] == 1) {
            table[i][1]++;
          } else {
            table[i][0]++;
          }
        }
      }
    }
  }
  
  double sum_Ey_linha;
  double alpha, betha;
  int M, N,s;
  alpha = 1;
  s = samples;
  M = num_comb;
  N = 0;
  sum_Ey_linha = 0.0;
  tgt_one = 0;
  tgt_zero = 0;
  for(int i = 0; i < num_comb; i++)
  {
    tgt_one += table[i][1];
    tgt_zero += table[i][0];
    int sum_iesima_comb = (table[i][0] + table[i][1]);
    if(sum_iesima_comb > 0)
    {
      N++;

      if(table[i][1] > table[i][0])
      {
        sum_Ey_linha += ((double)table[i][0]/(double)sum_iesima_comb)*(double)(sum_iesima_comb+alpha);
      }
      else
      {
        sum_Ey_linha += ((double)table[i][1]/(double)sum_iesima_comb)*(double)(sum_iesima_comb+alpha);
        
      }
    }
  }
  // Compute epsilon_target
  epsilon_target = (tgt_one > tgt_zero)? static_cast<double>(tgt_zero) : static_cast<double>(tgt_one);
  epsilon_target /= static_cast<double>(samples);
  // Compute epsilon_Z
  epsilon_Z = 0.0;
  epsilon_Z = (double)((1.0/(alpha*M+s))*((alpha*(M-N)*epsilon_target)+sum_Ey_linha));
  
  if (epsilon_target == 0.0 && epsilon_Z == 0.0) {
    return 1.0;
  } else if ((epsilon_target == 0.0 && epsilon_target < epsilon_Z) || (epsilon_target > 0.0 && epsilon_target < epsilon_Z)) {
    return 0.0;
  } else {
    return (epsilon_target - epsilon_Z)/epsilon_target;
  }
}


/* Mutual information based on Tsallis entropy */
double IFFS::IM(set<int> Sub_Candidatos, int target)
{
  int num_comb, i, t, bin;
  vector<vector<int> > table;
  vector<int> v;
  set<int>::iterator it;
  long double q, tsallis, x;
  double samples = 0.0;
  

  /* Initializing the conditional probability table. */
  num_comb = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
  v.assign(2, 0);
  for (i = 0; i < num_comb; ++i) {
    table.push_back(v);
  }
  
  /* nsamples-1 because we are considering time series data. */
  for (t = 0; t < this->nsamples-1; ++t) {
    i = 0;
    bin = 1;
    for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
      i += this->M[*it][t] * bin;
      bin *= 2;
    }
    
    if(this->Chip_Set[t][1] == -1) //estático
    {
      samples += 1.0;
      if(this->M[target][t] == 1)
        table[i][1]++;
      else
        table[i][0]++;
    }
    else
    {
      if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
      {
        if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
        {
          samples += 1.0;
          if (this->M[target][t+1] == 1) {
            table[i][1]++;
          } else {
            table[i][0]++;
          }
        }
      }
    }
  }
  double Hx, Hy, Hxy, py0_xi, py1_xi, Hy_x, Hy_xi;
  double py0, py1, px;
  py0 = py1 = py0_xi = py1_xi =  0.0;
  Hx = Hy = Hxy = Hy_x = Hy_xi = 0.0;
  for(int i = 0; i < num_comb; i++)
  {
    py1 += static_cast<double>(table[i][1]);
    py0 += static_cast<double>(table[i][0]);
    
    double px = (table[i][0]+table[i][1])/samples;
    double pxy1 = table[i][1]/samples;
    double pxy0 = table[i][0]/samples;
    int tm = (table[i][0]+table[i][1]);
    Hy_xi = 0;
    
    
    if(tm > 0)
    {
      py1_xi = (double)table[i][1]/tm;
      py0_xi = (double)table[i][0]/tm;
    }
    else //combinanção não observada (não existe essa combinação)
    {
      py1_xi = 0.0;
      py0_xi = 0.0;
    }
    
    if(px != 0 && !ISNAN(px))
      Hx += (double)(px * log2l(px));
    
    //entropia conjunta
    if(pxy1 != 0 && !ISNAN(pxy1))
    {
      Hxy += (double)(pxy1 * log2l(pxy1));
    }
    if(pxy0 != 0  && !ISNAN(pxy0))
    {
      Hxy += (double)(pxy0 * log2l(pxy0));
    }
    
    //entropia condicional
    if(py1_xi != 0 && !ISNAN(py1_xi))
    {
      Hy_xi += (double)(py1_xi * log2(py1_xi));
    }
    if(py0_xi != 0  && !ISNAN(py0_xi))
    {
      Hy_xi += (double)(py0_xi * log2(py0_xi));
    }
    
    Hy_xi *= -1.0;
    Hy_x += (double)(px*Hy_xi);
  }
  
  
  
  Hx *= -1.0;
  
  Hxy *= -1.0;
  
  py1 =(double)( py1 / samples);
  py0 = (double)(py0 / samples);
  
  if(py1 != 0 && !ISNAN(py1))
    Hy = (double)(py1 * log2l(py1));
  if(py0 != 0 && !ISNAN(py0))
    Hy += (double)(py0 * log2l(py0));
  if(Hy != 0)
    Hy *= -1.0;
  
  double im2 = (double)(Hy - Hy_x);
  double im = (double)(Hx + Hy - Hxy);
 // Rcout << "IM = " << im << " IM2 = " << im2 << endl;
  
  if (im < 0.0) {
    im = 0;
  }
  
  return im;
}


/* Mutual information based on Tsallis entropy */
double IFFS::IM_pen_rara(set<int> Sub_Candidatos, int target)
{
  int num_comb, i, t, bin;
  vector<vector<int> > table;
  vector<int> v;
  set<int>::iterator it;
  long double q, tsallis, x;
  double samples = 0.0;
  
  
  /* Initializing the conditional probability table. */
  num_comb = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
  v.assign(2, 0);
  for (i = 0; i < num_comb; ++i) {
    table.push_back(v);
  }
  
  /* nsamples-1 because we are considering time series data. */
  for (t = 0; t < this->nsamples-1; ++t) {
    i = 0;
    bin = 1;
    for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
      i += this->M[*it][t] * bin;
      bin *= 2;
    }
    
    if(this->Chip_Set[t][1] == -1) //estático
    {
      samples += 1.0;
      if(this->M[target][t] == 1)
        table[i][1]++;
      else
        table[i][0]++;
    }
    else
    {
      if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
      {
        if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
        {
          samples += 1.0;
          if (this->M[target][t+1] == 1) {
            table[i][1]++;
          } else {
            table[i][0]++;
          }
        }
      }
    }
  }
  double Hy, py0_xi, py1_xi, Hy_x, Hy_xi, Hy_x_pen_raro;
  double py0, py1, px;
  py0 = py1 = py0_xi = py1_xi =  0.0;
  Hy = Hy_x = Hy_xi = Hy_x_pen_raro = 0.0;
  
  int M, N, s;
  double alpha, betha;
  M = num_comb;
  N = 0;
  s = samples;
  alpha = 1;
  betha = 0.8;
  double HY_X = 0;
  
  for(int i = 0; i < num_comb; i++)
  {
    py1 += static_cast<double>(table[i][1]);
    py0 += static_cast<double>(table[i][0]);
    
    int sum_iesima_comb = (table[i][0] + table[i][1]);
    if(sum_iesima_comb > 1)
    {
      N++;
      Hy_xi = 0;
      py1_xi = (double)table[i][1]/(double)sum_iesima_comb;
      py0_xi = (double)table[i][0]/(double)sum_iesima_comb;
      //entropia condicional
      if(py1_xi != 0 && !ISNAN(py1_xi))
      {
        Hy_xi += (double)(py1_xi * log2(py1_xi));
      }
      if(py0_xi != 0  && !ISNAN(py0_xi))
      {
        Hy_xi += (double)(py0_xi * log2(py0_xi));
      }
      px = sum_iesima_comb/samples;
      Hy_xi *= -1.0;
      Hy_x += (double)(px*Hy_xi);
    }
  }
  
  int c = 2;
  double H_delta = (double)-1*((betha * log2(betha))+(double)(((1-betha)/(c-1)) * log2((1-betha)/(c-1))));
  
  Hy_x_pen_raro = ((((double)(M-N)/s) * H_delta)  + (double)Hy_x);
  
  for(int i = 0; i < num_comb; i++)
  {
    px = (table[i][0]+table[i][1])/samples;

    HY_X += (double)(px*Hy_x_pen_raro);
  }
  
  py1 =(double)( py1 / samples);
  py0 = (double)(py0 / samples);
  
  if(py1 != 0 && !ISNAN(py1))
    Hy = (double)(py1 * log2l(py1));
  if(py0 != 0 && !ISNAN(py0))
    Hy += (double)(py0 * log2l(py0));
  if(Hy != 0)
    Hy *= -1.0;
  

  double im3 = (double)(Hy - HY_X);

  if (im3 < 0.0) {
    im3 = 0;
  }

  return im3;
}


/* Mutual information based on Tsallis entropy */
double IFFS::IM_pen_zero(set<int> Sub_Candidatos, int target)
{
  int num_comb, i, t, bin;
  vector<vector<int> > table;
  vector<int> v;
  set<int>::iterator it;
  long double q, tsallis, x;
  double samples = 0.0;
  
  
  /* Initializing the conditional probability table. */
  num_comb = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
  v.assign(2, 0);
  for (i = 0; i < num_comb; ++i) {
    table.push_back(v);
  }
  
  /* nsamples-1 because we are considering time series data. */
  for (t = 0; t < this->nsamples-1; ++t) {
    i = 0;
    bin = 1;
    for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
      i += this->M[*it][t] * bin;
      bin *= 2;
    }
    
    if(this->Chip_Set[t][1] == -1) //estático
    {
      samples += 1.0;
      if(this->M[target][t] == 1)
        table[i][1]++;
      else
        table[i][0]++;
    }
    else
    {
      if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
      {
        if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
        {
          samples += 1.0;
          if (this->M[target][t+1] == 1) {
            table[i][1]++;
          } else {
            table[i][0]++;
          }
        }
      }
    }
  }
  double Hy, py0_xi, py1_xi, Hy_xi, Hy_x;
  double py0, py1, px;
  px = py0 = py1 = py0_xi = py1_xi =  0.0;
  Hy = Hy_xi = Hy_x = 0.0;
  
  int M, N, s;
  double alpha, betha;
  M = num_comb;
  N = 0;
  s = samples;
  alpha = 1;
  betha = 0.8;

  for(int i = 0; i < num_comb; i++)
  {
    int sum_iesima_comb = (table[i][0] + table[i][1]);
    py1 += static_cast<double>(table[i][1]);
    py0 += static_cast<double>(table[i][0]);
    px = (table[i][0]+table[i][1])/samples;
    if(sum_iesima_comb > 0)
    {
      N++;
      Hy_xi = 0;
      py1_xi = (double)table[i][1]/(double)sum_iesima_comb;
      py0_xi = (double)table[i][0]/(double)sum_iesima_comb;

      if(py1_xi != 0 && !ISNAN(py1_xi))
      {
        Hy_xi += (double)(py1_xi * log2(py1_xi));
      }
      if(py0_xi != 0  && !ISNAN(py0_xi))
      {
        Hy_xi += (double)(py0_xi * log2(py0_xi));
      }
      
      Hy_xi *= -1.0;

      Hy_x += (double)(sum_iesima_comb+alpha)* px * (Hy_xi);   
    }
  }
 
  
  py1 =(double)( py1 / samples);
  py0 = (double)(py0 / samples);
  
  if(py1 != 0 && !ISNAN(py1))
    Hy = (double)(py1 * log2l(py1));
  if(py0 != 0 && !ISNAN(py0))
    Hy += (double)(py0 * log2l(py0));
  if(Hy != 0)
    Hy *= -1.0;

  double Hy_x_pen_zero = (double)(1/((alpha*M)+s))*((alpha*(M-N)*Hy) + (double)Hy_x);
  

  double im3 = (double)(Hy - Hy_x_pen_zero);
  
  if (im3 < 0.0) {
    im3 = 0;
  }
  return im3;
}

/* Mutual information based on Tsallis entropy */
double IFFS::IM_Tsallis(set<int> Sub_Candidatos, int target)
{
	int m, i, t, bin;
	vector<vector<int> > table;
	vector<int> v;
	set<int>::iterator it;
	long double q, tsallis, x;
	double samples = 0.0;
	
	q = 2.5;

	/* Initializing the conditional probability table. */
	m = static_cast<int>(pow(2.0, static_cast<double>(Sub_Candidatos.size())));
	v.assign(2, 0);
	for (i = 0; i < m; ++i) {
		table.push_back(v);
	}
	
	/* nsamples-1 because we are considering time series data. */
	for (t = 0; t < this->nsamples-1; ++t) {
		i = 0;
		bin = 1;
		for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
			i += this->M[*it][t] * bin;
			bin *= 2;
		}

		if(this->Chip_Set[t][1] == -1) //estático
		{
		  samples += 1.0;
		  if(this->M[target][t] == 1)
		    table[i][1]++;
		  else
		    table[i][0]++;
		}
		else
		{
		  if (this->Chip_Set[t][0] == this->Chip_Set[t+1][0]) 
		  {
		    if (this->Chip_Set[t][1] < this->Chip_Set[t+1][1]) 
		    {
		      samples += 1.0;
		      if (this->M[target][t+1] == 1) {
		        table[i][1]++;
		      } else {
		        table[i][0]++;
		      }
		    }
		  }
		}
	}

	tsallis = TsallisEntropy(table, m, q);
	long double cond;
	cond = ConditionalTsallisEntropy(table, m, q, tsallis, samples);
	x = tsallis - cond;
	// Round-off error treatment
	if (x < 0.0) {
		
		if (x >= -0.1) {
			x = 0.0;
		} else {
			cout << "Error " << tsallis << " - " << cond << " = " << x << endl;	
			exit(1);
		}
		
	}
	return x;
}

double IFFS::TsallisEntropy(vector<vector<int> >table, int m, long double q)
{
	int i;
	long double pr_0, pr_1;

	pr_0 = pr_1 = 0.0;
	for (i = 0; i < m; ++i) {
		pr_0 += static_cast<double>(table[i][0]);
		pr_1 += static_cast<double>(table[i][1]);
	}

	pr_0 = pr_0 / (pr_0 + pr_1);
	pr_1 = 1.0 - pr_0;

	return (pow(pr_1, q) + pow(pr_0, q) - 1.0) / (1.0 - q);
}

double IFFS::ConditionalTsallisEntropy(vector<vector<int> >table, int p, long double q, long double tsallis, double samples)
{
	long double m, n, d, alpha, rg, pr0_g, pr1_g, sum;
	int i;

	alpha = 1.0;                                 // Penalty weight   
	m = static_cast<long double>(p);                  // Number of possible instances
	d = static_cast<long double>(samples); // Number of samples

	n = sum = 0.0;
	for (i = 0; i < p; ++i) {

		// Observed instance
		if (table[i][0] + table[i][1] != 0) {
			n += 1.0;
			rg = static_cast<double>(table[i][0] + table[i][1]);

			pr0_g = static_cast<double>(table[i][0]) / rg;
			pr1_g = static_cast<double>(table[i][1]) / rg;

			sum += (rg + alpha)/(alpha * m + d) * (1.0 - pow(pr0_g, q) - pow(pr1_g, q));
		}
	}

	sum *= 1.0 / (q - 1.0);

	return (alpha * (m - n)) / (alpha * m + d) * tsallis + sum;
}

float IFFS::SDB(set<int> Sub_Candidatos, int target)
{
  float sdb = 0.0;
  int qnt = 0;
  set<int>::iterator it;
  
  for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
    sdb += (this->M_SDB[*it][target]);
    qnt++;
  }
  
  if(qnt > 0)
  {
    return (float)sdb/qnt;
  }
  return 0;
}

/* Sequential Backward Selection */
bool IFFS::SBS(int target, set<int>& All_Candidatos, set<int>& Sub_Candidatos)
{
	double best_score, cur_score;
	set<int> Z;
	set<int>::iterator it;
	int least_feat;

	if (Sub_Candidatos.size() <= 2) {
		return false;
	}

	/* Find the least significant feature. */
	best_score = this->score[Sub_Candidatos.size()-1];
	least_feat = -1;
	Z = Sub_Candidatos;
	for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
		Z.erase(*it);
		cur_score = this->J(Z, target);
		if (cur_score > best_score ) {
			best_score = cur_score;
			least_feat = *it;
		}
		Z.insert(*it);
	}

	/* If removing the least significant feature improves the currently score, then remove it. */
	if (least_feat != -1) {
		Sub_Candidatos.erase(least_feat);
		All_Candidatos.insert(least_feat);
		//cout << "SBS: remove " << least_feat << " improving from " << this->score[Sub_Candidatos.size()] << " to " << best_score << endl;
		//PrintSet(Sub_Candidatos);
		this->score[Sub_Candidatos.size()] = best_score;
		return true;
	}

	return false;
}

/* Replace step of IFFS algorithm */
bool IFFS::Replace(int target, set<int>& All_Candidatos, set<int>& Sub_Candidatos)
{
	double best_score, cur_score, max_score;
	set<int> Z;
	set<int>::iterator it, g;
	int best_feat, least_feat;

	if (Sub_Candidatos.size() < 2) {
		return false;
	}

	Z = Sub_Candidatos;
	max_score = -1.0;
	for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
		Z.erase(*it);
		for (g = All_Candidatos.begin(); g != All_Candidatos.end(); ++g) {
			Z.insert(*g);
			cur_score = this->J(Z, target);
			if (cur_score > max_score) {
				max_score = cur_score;
				best_feat = *g;
				least_feat = *it;
			}
			Z.erase(*g);
		}
		Z.insert(*it);
	}

	best_score = this->score[Sub_Candidatos.size()];
	if (max_score > best_score) {
		Sub_Candidatos.erase(least_feat);
		All_Candidatos.insert(least_feat);
		Sub_Candidatos.insert(best_feat);
		All_Candidatos.erase(best_feat);
		this->score[Sub_Candidatos.size()] = max_score;
		//cout << "Replace: " << least_feat << " by " << best_feat << " improves from " << best_score << " to " << max_score << endl;
		//PrintSet(Sub_Candidatos);
		return true;
	}

	return false;
}

void IFFS::WriteOutput(set<int> Sub_Candidatos)
{
  set<int>::iterator it;
  for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
    f_resultado_com_score << *it << "\t";
  }
  f_resultado_com_score << this->score[Sub_Candidatos.size()] << endl;
  
}

void IFFS::salva_topX(std::vector<int> Sub_Candidatos, double crit_func_eval)
{
  for(int y = 0; y < top_X_preditors.nrow(); y++)
  {
    if(crit_func_eval > top_X_preditors(y,max_k))
    {
      for(int h = top_X_preditors.nrow()-1; h > y; h--)
      {
        for(int r = 0; r <= max_k; r++)
        {
          top_X_preditors(h,r) = top_X_preditors(h-1,r);
        }
        
      }
      for(int r = 0; r < Sub_Candidatos.size(); r++)
      {
        top_X_preditors(y,r) = Sub_Candidatos[r];
      }
      top_X_preditors(y,max_k) = crit_func_eval;
      break;
    }
  }
};

void IFFS::Write_top_X_preditores(int target)
{
  std::ofstream fout;
  std::stringstream sstm;
  sstm << output << target << ".txt";
  std::string fname = sstm.str();
  fout.open(fname.c_str(),std::ios::out);
  
  if (!fout.is_open()) {
    Rcpp::Rcout << "Could not open " << fname.c_str() << endl;
    exit(1);
  }
  
  for(int h = 0; h < top_X_preditors.nrow(); h++)
  {
    for(int r = 0; r <= max_k; r++)
    {
      fout << top_X_preditors(h,r);
      if(r != max_k)
        fout << " ";
    }
    fout << endl;
  }
  fout.close();
};

/* ======================================================== */
/* Auxiliary methods */
/* ======================================================== */

void IFFS::PrintSet(set<int> Sub_Candidatos)
{
	set<int>::iterator it;
	cout << "{ ";
	for (it = Sub_Candidatos.begin(); it != Sub_Candidatos.end(); ++it) {
		cout << *it << " ";
	}
	cout << "}" << endl;
}

void IFFS::PrintMatrix(vector<vector<int> > T)
{
	vector<vector<int> >::iterator iti;
	vector<int>::iterator itj;
	for (iti = T.begin(); iti != T.end(); ++iti) {
		for (itj = (*iti).begin(); itj != (*iti).end(); ++itj) {
			cout << *itj << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void IFFS::PrintMatrix(float** T)
{
  int iti,itj;
  for (iti = 0; iti < this->ngenes; ++iti) {
    for (itj = 0; itj < this->ngenes; ++itj) {
      cout << "[" << iti << "]" << "[" <<  itj << "]" << " = "<< T[iti][itj] << "\n";
    }
    cout << endl;
  }
  cout << endl;
}


// [[Rcpp::export]]
void RunIffs(SEXP InputExpr, SEXP InputSDB, SEXP InputCS, SEXP output, SEXP resultado_com_score, int max_k, int delta,
            float w_cod, float w_cod_rara, float w_cod_zero, float w_im, float w_im_rara, float w_im_zero, float w_tsallis, float w_sdb)
{
  //Cria o objeto IFFS e a matrix de preditores
  
  IFFS *f;
  
  
  f = new IFFS(InputExpr, InputSDB,  InputCS, output, resultado_com_score,  max_k,  delta,
               w_cod,  w_cod_rara,  w_cod_zero,  w_im,  w_im_rara,  w_im_zero,  w_tsallis,  w_sdb);
  
  delete f;

  
}

