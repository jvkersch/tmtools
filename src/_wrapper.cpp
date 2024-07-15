#include <iostream>
#include <string>
#include <vector>

#include "_tmalign.h"

void _tmalign_wrapper(double **xa, double **ya, const char *seqx,
                      const char *seqy, const int xlen, const int ylen,
                      double t0[3], double u0[3][3], double &TM1, double &TM2, double &rmsd0,
                      std::string &seqM, std::string &seqxA, std::string &seqyA)

{
  // user alignment
  std::vector<std::string> sequence; // get value from alignment file

  // sequences and coordinates
  char *secx = new char[xlen + 1];
  char *secy = new char[ylen + 1];

  make_sec(xa, xlen, secx);
  make_sec(ya, ylen, secy);

  // parameters (currently hardcoded)
  int i_opt = 0;         // 1 for -i, 3 for -I
  int a_opt = 0;         // flag for -a, do not normalized by average length
  bool u_opt = false;    // flag for -u, normalized by user specified length
  bool d_opt = false;    // flag for -d, user specified d0
  bool fast_opt = false; // flags for -fast, fTM-align algorithm
  double Lnorm_ass = 0, d0_scale = 0; // options that can be set by user

  int mol_type = 0;

  // output variables
  double TM3, TM4, TM5; // for a_opt, u_opt, d_opt
  double d0_0, TM_0;
  double d0A, d0B, d0u, d0a;
  double d0_out = 5.0;
  int L_ali; // Aligned length in standard_TMscore
  double Liden = 0;
  double TM_ali; // TMscore in standard_TMscore
  double rmsd_ali; // TMscore in standard_TMscore
  int n_ali = 0;
  int n_ali8 = 0;

  TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0, TM1, TM2, TM3, TM4, TM5,
               d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
               rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8, xlen, ylen,
               sequence, Lnorm_ass, d0_scale, i_opt, a_opt, u_opt, d_opt,
               fast_opt, mol_type);

  delete[] secx;
  delete[] secy;
}
