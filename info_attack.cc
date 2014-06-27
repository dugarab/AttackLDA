#include "info_attack.h"

namespace attack_lda {
/*
 * read settings.
 *
 */
void ReadAttackSettings(char* filename , ParamAttack *paramattack )
{
    FILE* fileptr;
    char alpha_action[100];
    fileptr = fopen(filename, "r");
	paramattack->lag = 5;
    fscanf(fileptr, "var max iter %d\n", &paramattack->e_step_limit_iterations);
    fscanf(fileptr, "var convergence %f\n", &paramattack->e_step_converge_threshold);
    fscanf(fileptr, "em max iter %d\n", &paramattack->em_limit_iterations);
    fscanf(fileptr, "em convergence %f\n", &paramattack->em_converge_threshold);
    fscanf(fileptr, "alpha %s\n", alpha_action);
	if( strcmp( alpha_action , "fixed" ) == 0 )
		paramattack->if_estimate_alpha = 0;
	else 
		paramattack->if_estimate_alpha = 1;
	fscanf(fileptr, "total attack words %d\n" , &paramattack->num_attack_words);
	fscanf(fileptr , "max attack iteration %d\n" , &paramattack->limit_attack_iterations);
	fscanf(fileptr , "limit words per doc %d\n" , &paramattack->limit_perdoc_attackwords );
	
	char tmp_word[300];
	char tmp_topictheme[300];
	int tmp_sortid;
	while( fscanf(fileptr , "inject word %s\n" , tmp_word ) == 1 )
	{
		fscanf(fileptr , "attack topic %s\n" , tmp_topictheme );
		fscanf(fileptr , "word sortid attack %d\n" , &tmp_sortid );
		paramattack->vector_attack_wordstr.push_back( tmp_word );
		paramattack->vector_attack_topictheme.push_back( tmp_topictheme );
		paramattack->vector_attack_word_sortid.push_back( tmp_sortid );
	}
	fclose(fileptr);
}
/*
* read the total effort and the estimated effort 
* ONLY for random attacker
*/
void ReadEffortSettings(char* filename , ParamAttack *paramattack)
{
	FILE* fileptr;
    char word[100];
    fileptr = fopen(filename, "r");
	int n_att = paramattack->vector_attack_wordstr.size();
	paramattack->vector_expected_attack_effort.resize( n_att );
	double sum_effort = 0.0;
	double tmp_effort;

	fscanf( fileptr , "%lf\n" , &sum_effort );
	for(int i = 0; i < n_att; i++)
	{
		fscanf( fileptr , "%s %lf\n" , word , &tmp_effort );
		for(int j = 0; j < n_att; j++)
		{
			printf( "word %s and real word %s\n" , word , paramattack->vector_attack_wordstr[j].c_str() );
		
			if( strcmp( paramattack->vector_attack_wordstr[j].c_str() , word ) == 0 )
			{
				printf( "found pair %d for word %s\n" , j , word );
				paramattack->vector_expected_attack_effort[j] = tmp_effort;
			}
		}
	}
	double real_sum_effort = 0.0;
	for(int i = 0; i < n_att; i++)
		real_sum_effort += fabs(paramattack->vector_expected_attack_effort[i]);
	for(int i = 0; i < n_att; i++)
		paramattack->vector_expected_attack_effort[i] *= sum_effort / real_sum_effort;

	fclose(fileptr);	
}

/*
 * print the setting parameters to the outinfo_file
 *
 */
void PrintAttackSettings( FILE *fout , ParamAttack *paramattack )
{
	fprintf(fout, "var max iter %d\n", paramattack->e_step_limit_iterations);
    fprintf(fout, "var convergence %f\n", paramattack->e_step_converge_threshold);
    fprintf(fout, "em max iter %d\n", paramattack->em_limit_iterations);
    fprintf(fout, "em convergence %f\n", paramattack->em_converge_threshold);
	fprintf(fout , "ALPHA = %0.5lf" , paramattack->initial_alpha_value );
	fprintf(fout , "BETA = %0.5lf" , paramattack->initial_beta_value );

	fprintf(fout, "total attack words %d\n" , paramattack->num_attack_words);
	fprintf(fout , "max attack iteration %d\n" , paramattack->limit_attack_iterations);
	fprintf(fout , "limit words per doc %d\n" , paramattack->limit_perdoc_attackwords );
	for(int i = 0; i < paramattack->vector_attack_wordid.size(); i++)
	{
		fprintf(fout , "inject wordid %d in topicid %d\n" 
					, paramattack->vector_attack_wordid[i] , paramattack->vector_attack_topicid[i] );
	}
}

/*
 * writes the word assignments line for a document to a file
 *
 */

void SaveWordAssignment(FILE* f, Document* doc, double** phi, LdaModel* model)
{
    int n;

    fprintf(f, "%03d", doc->length);
    for (n = 0; n < doc->length; n++)
    {
        fprintf(f, " %04d:%02d",
                (*doc->words)[n], IndexOfMax(phi[n], model->num_topics_));
    }
    fprintf(f, "\n");
    fflush(f);
}

/*
 * saves the gamma parameters of the current dataset
 *
 */
void SaveGamma(char* filename, double** gamma, int num_docs, int num_topics)
{
    FILE* fileptr;
    int d, k;
    fileptr = fopen(filename, "w");

    for (d = 0; d < num_docs; d++)
    {
	fprintf(fileptr, "%5.10f", gamma[d][0]);
	for (k = 1; k < num_topics; k++)
	{
	    fprintf(fileptr, " %5.10f", gamma[d][k]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}


/*
* print top words of topics log_beta 
* @param model  sorted_topics to sort the log_beta,
* @param topicname  the name of the topics: "original" or "mislead" 
*/
void PrintTopWordsInTopic( FILE* fout , double** log_beta , int num_topwords , LdaModel* model , Corpus* corpus , char* topicname )
{
	//printf( "start print toptopic log_beta = %d , model = %d" ,log_beta , model );
	int k = 0 , i = 0;
	model->SortWordProbEachTopic(log_beta);

	fprintf(fout , "%s\n" , topicname);
	for(k = 0; k < model->num_topics_; k++)
	{
		double norm = 0.0;
		for(i = 0; i < model->num_terms_; i++)
			norm += model->sorted_topics_[k][i].proportion;
	 
		for(i = 0; i < num_topwords; i++)
		{
			fprintf(fout , "%s:%0.5lf " , corpus->dict[model->sorted_topics_[k][i].word_id] , model->sorted_topics_[k][i].proportion/norm );
		}
		fprintf(fout , " TotalWeight=%0.5lf\n" , norm);
	}
	fprintf(fout , "\n");
}

/*
* print two measures : epsilon insensitive l2 distance and the proportion of target word and Rank in target topic
* @param param_attack : the attack setting
*/
void PrintAttackMeasure( FILE* fout , Corpus* corpus , LdaModel* model , ParamAttack* pm_att )
{
	fprintf(fout , "epsilon_insen_l2_dis( phi* | mean(eta) ) = %0.10lf\n" , model->L2DistanceMisleadEstimateTopics(fout, pm_att) );
	double sum_effort = 0.0;
	for(int d = 0; d < corpus->num_docs; d++)
		sum_effort += corpus->docs[d].attack_total;
	fprintf( fout , "sum_effort = %0.5lf\n" , sum_effort );

	// output the proportion of attack word in the topic
	for(int i = 0; i < pm_att->vector_attack_wordid.size(); i++)
	{
		double norm = model->sumeta_[pm_att->vector_attack_topicid[i]];
		double total_attackword = 0.0;
		for(int k = 0; k < model->num_topics_; k++) 
			total_attackword += model->eta_[k][pm_att->vector_attack_wordid[i]];

		fprintf(fout , "AttackLog (word = %s) TotalProb = %0.5lf Prob_in(attack_topic %d) = %0.5lf Rank = %d Proportion = %0.5lf\n" 
								, corpus->dict[pm_att->vector_attack_wordid[i]] 
								, total_attackword
								, pm_att->vector_attack_topicid[i]
								, model->eta_[pm_att->vector_attack_topicid[i]][pm_att->vector_attack_wordid[i]]
								, WordRankInTopic(model->num_terms_ , pm_att->vector_attack_topicid[i] , pm_att->vector_attack_wordid[i] , model->eta_)
								, model->eta_[pm_att->vector_attack_topicid[i]][pm_att->vector_attack_wordid[i]] / norm
			);
	}
}

/*
* print the attack documents in words (not in word index)
*/
void PrintDocumentWordStr( FILE* fout_selectdocs , int doc_id , Document doc , char** dict )
{
	fprintf( fout_selectdocs , "DocId=%d" , doc_id );
	for(int i = 0; i < doc.length; i++)
		fprintf( fout_selectdocs , " %s:%0.3lf" , dict[(*doc.words)[i]] , (*doc.counts)[i] );
	fprintf( fout_selectdocs , "\n" );
}

} // namespace attack_lda