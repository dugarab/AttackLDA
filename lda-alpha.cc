// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "lda-alpha.h"

namespace attack_lda {
/*
 * objective function and its derivatives
 *
 */
double AlphaLikelihood(double a, double ss, int DD, int KK)
{ return(DD * (lgamma(KK * a) - KK * lgamma(a)) + (a - 1) * ss); }

double DAlphaLikelihood(double a, double ss, int DD, int KK)
{ return(DD * (KK * Digamma(KK * a) - KK * Digamma(a)) + ss); }

double D2AlphaLikelihood(double a, int DD, int KK)
{ return(DD * (KK * KK * Trigamma(KK * a) - KK * Trigamma(a))); }


/*
 * newtons method
 *
 */

double OptimalAlpha(double ss, int DD, int KK)
{
    double a, log_a, init_a = 100;
    double f, df, d2f;
    int iter = 0;

    log_a = log(init_a);
    do
    {
        iter++;
        a = exp(log_a);
        if (isnan(a))
        {
            init_a = init_a * 10;
            printf("warning : alpha is nan; new init = %5.5f\n", init_a);
            a = init_a;
            log_a = log(a);
        }
        f = AlphaLikelihood(a, ss, DD, KK);
        df = DAlphaLikelihood(a, ss, DD, KK);
        d2f = D2AlphaLikelihood(a, DD, KK);
        log_a = log_a - df/(d2f * a + df);
        printf("alpha maximization : %5.5f   %5.5f\n", f, df);
    }
    while ((fabs(df) > NEWTON_THRESH) && (iter < MAX_ALPHA_ITER));
    return(exp(log_a));
}
} // namespace attack_lda