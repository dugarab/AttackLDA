#ifndef ATTACKLDA_UTILTOPICS_H_
#define ATTACKLDA_UTILTOPICS_H_

#include "lda.h"
#include "utils.h"
#include "gibbs-sampling.h"
#include <cstdio>
#include <queue>
#include <vector>

using std::vector;
using std::priority_queue;

namespace attack_lda {

typedef struct {
  int attack_word;
  int attack_doc;
  double delta_w; 
} AttackPair;

class AttackPairsLessthanByGain {
 public:
  bool operator() (const  AttackPair& lhs, const AttackPair& rhs) const {
    return (fabs(lhs.delta_w) > fabs(rhs.delta_w));
  }
};

class AttackPairsMorethanByGain {
 public:
  bool operator() (const AttackPair& lhs, const AttackPair& rhs) const {
    return (fabs(lhs.delta_w) < fabs(rhs.delta_w));
  }
};

/** 
  * @brief pop elements in qp to control the size <= lower_size
  * @param lower_size 
  * @param qp
  * @return 
  * @type 
   */
void ResizeQueueAttackPairs(int lower_size, priority_queue<AttackPair, vector<AttackPair>, AttackPairsLessthanByGain>* qp);

/** 
  * @brief insert q1 into q2 and delete q1
  * @param q1
  * @param q2
  * @return 
  * @type 
  */
void InsertQueueAttackPairs(priority_queue<AttackPair, vector<AttackPair>, AttackPairsLessthanByGain>* q1, 
							 priority_queue<AttackPair, vector<AttackPair>, AttackPairsLessthanByGain>* q2);

/** 
  * @brief copy the minimal root heap to maximal root heap
  * @param min_q
  * @param max_q
  * @return 
  * @type 
  */
void CopyToMaxQueueFromMinQueue(priority_queue<AttackPair, vector<AttackPair>, AttackPairsLessthanByGain>* min_q, 
								 priority_queue<AttackPair, vector<AttackPair>, AttackPairsMorethanByGain>* max_q);

int CompareWordWeightPair (const void * a, const void * b);

class LdaModel {
 public:
  double alpha_;
  double beta_;
  double** eta_;
  double* sumeta_;
  double** init_eta_; // initial value of eta
  double** phi_;

  WordWeightPair** sorted_topics_;
  double** mislead_probw_;
  double** log_origin_topics_;
  double** gamma_; // variational parameters for topic proportion
  int num_documents_;
  int num_topics_;
  int num_terms_;

  GibbsSampling gibbs_sampling_;

  void FindBestTopicPermutation(int* permute);
  void PermuteMisleadTopics(int* permute);
  double L2DistanceMisleadEstimateTopics(FILE* fout, ParamAttack* pm_att);
  void SortWordProbEachTopic(double** log_beta);
  double DifferenceSpecificWordOriginMisleadTopic(int w);

  void InitializeSuffstatsZero(LdaSuffstats* ss);

  void CalcDigammaSumEta(std::vector<double>& digamma_sumeta);

  double ExpectationStepOneDoc(Document* doc, int d, std::vector<double>& digamma_sumeta, ParamAttack* pm_att, LdaSuffstats* ss);

  void RunVariationalInference(char* directory, Corpus* corpus, ParamAttack *pm_att, LdaSuffstats* ss );

  double ComputeVariationalProb();

  void LdaMle(LdaSuffstats* ss);

  double LdaInference(Document*, int d, ParamAttack* pm_att, std::vector<double>& digamma_sumeta);

  double ComputeLikelihood(Document*, int d, std::vector<double>& digamma_sumeta);

  void FreeLdaModel();
  void SaveLdaModel(char*);
};

} // namespace attack_lda
#endif
