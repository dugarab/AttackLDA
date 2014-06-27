#ifndef ATTACKLDA_INFOATTACK_H_
#define ATTACKLDA_INFOATTACK_H_
#include <cstdio>
#include <cstdlib>
#include "lda.h"
#include "util_topics.h"

namespace attack_lda {
void ReadAttackSettings(char* filename , ParamAttack *paramattack);
void ReadEffortSettings(char* filename , ParamAttack *paramattack);
void PrintAttackSettings( FILE* fout , ParamAttack *paramattack);
void SaveWordAssignment(FILE* f, Document* doc, double** phi, LdaModel* model);
void SaveGamma(char* filename, double** gamma, int num_docs, int num_topics);
void PrintTopWordsInTopic( FILE* fout , double** log_beta , int num_topwords , LdaModel* model , Corpus* corpus , char* topicname );
void PrintAttackMeasure( FILE* fout , Corpus* corpus , LdaModel* model , ParamAttack* pm_att );
void PrintDocumentWordStr( FILE* fout_selectdocs , int doc_id , Document doc , char** dict );
} // namespace attack_lda
#endif