#ifndef ATTACKLDA_LDAALPHA_H_
#define ATTACKLDA_LDAALPHA_H_

#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "lda.h"
#include "utils.h"

#define NEWTON_THRESH 1e-5
#define MAX_ALPHA_ITER 1000

namespace attack_lda {
double AlphaLikelihood(double a, double ss, int DD, int KK);
double DAlphaLikelihood(double a, double ss, int DD, int KK);
double D2AlphaLikelihood(double a, int DD, int KK);
double OptimalAlpha(double ss, int DD, int KK);
} // namespace attack_lda
#endif
