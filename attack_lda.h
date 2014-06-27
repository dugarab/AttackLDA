#ifndef ATTACKLDA_ATTACKLDA_H_
#define ATTACKLDA_ATTACKLDA_H_

#include "info_attack.h"
#include "util_topics.h"
#include "lda-estimate.h"

#include <vector>
#include <queue>

namespace attack_lda {

/**
 * @brief 
 * @param 
 * @param 
 *
 */
void ComputeMisleadTopics( Corpus* corpus, ParamAttack* pm_att , FILE* fout,  LdaModel* model );
void AttackAllPairs( Corpus* corpus, LdaModel* model , ParamAttack* pm_att, 
							std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q );
void AttackPairsInQueue( Corpus* corpus, LdaModel* model , ParamAttack* pm_att, 
							std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q );
void LossChangeAttackPair( Corpus* corpus , LdaModel* model , ParamAttack* pm_att, 
	 					std::vector<double>& digamma_sumeta , std::vector<double>& phinew ,
						AttackPair* tmp_attack_p );
void ModifyCorpus( AttackPair attpair , Corpus* corpus );
void EstimateGradientStepLength( Corpus* corpus , LdaModel* model , std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q ,
							ParamAttack* paramattack );
void RandomAttack( Corpus* corpus , ParamAttack* pm_att, 
									std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q );
void AttackLda(int argc, char* argv[]);

} // namespace attack_lda
#endif // ATTACKLDA_ATTACKLDA_H_
