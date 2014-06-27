#ifndef ATTACKLDA_LDAESTIMATE_H_
#define ATTACKLDA_LDAESTIMATE_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <algorithm>

#include "lda.h"
#include "lda-data.h"
#include "lda-alpha.h"
#include "utils.h"
#include "cokus.h"
#include "info_attack.h"

namespace attack_lda {

class InitialInfo {
 public:
  LdaSuffstats* ss;
  LdaModel* model;
 public:
  void InitializeLdaModel(int num_attackwords, char* start, char* directory, Corpus* corpus, ParamAttack *pm_att);
  void InitializeSuffstatsCorpus(Corpus* c);
  void InitializeSuffstatsRandom();
  LdaSuffstats* NewLdaSuffstats(LdaModel* model);
  LdaModel* NewLdaModel(int num_topics, Corpus* corpus, int num_attackwords, double alpha, double beta);
  LdaModel* LoadLdaModel(char* model_root, Corpus* corpus, int num_attackwords);
};
} // namespace attack_lda;
#endif


