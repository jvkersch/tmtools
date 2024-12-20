#ifndef _TMALIGN_H_
#define _TMALIGN_H_

void make_sec(double **x, int len, char *sec);

int TMalign_main(double **xa, double **ya, const char *seqx, const char *seqy,
                 const char *secx, const char *secy, double t0[3],
                 double u0[3][3], double &TM1, double &TM2, double &TM3,
                 double &TM4, double &TM5, double &d0_0, double &TM_0,
                 double &d0A, double &d0B, double &d0u, double &d0a,
                 double &d0_out, std::string &seqM, std::string &seqxA,
                 std::string &seqyA, double &rmsd0, int &L_ali, double &Liden,
                 double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
                 const int xlen, const int ylen,
                 const std::vector<std::string> sequence,
                 const double Lnorm_ass, const double d0_scale, const int i_opt,
                 const int a_opt, const bool u_opt, const bool d_opt,
                 const bool fast_opt, const int mol_type,
                 const double TMcut = -1);

void print_version();

#endif // _TMALIGN_H_
