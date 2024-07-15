#ifndef _TMALIGN_WRAPPER_H_
#define _TMALIGN_WRAPPER_H_

void _tmalign_wrapper(double **xa, double **ya, const char *seqx,
                      const char *seqy, const int xlen, const int ylen,
                      double t0[3], double u0[3][3], double &TM1, double &TM2, double &rmsd0,
                      std::string &seqM, std::string &seqxA, std::string &seqyA);

#endif // _TMALIGN_WRAPPER_H_