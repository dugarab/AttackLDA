#include "attack_lda.h"

namespace attack_lda { 

  ParamAttack paramattack;
  std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain> pq;
  std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsMorethanByGain> max_q;
  const int kLimitSizeQueue = 50000;
  const double kEpsilon = 1e-10;

  /*
  * get the mislead topics from original topics, attack words, from topic, to topic
  *
  */
  void ComputeMisleadTopics(Corpus* corpus,  ParamAttack* pm_att, FILE* fout, LdaModel* model) {
    model->sorted_topics_ = (WordWeightPair**) malloc(sizeof(WordWeightPair*)*(model->num_topics_));
    // sort the learned topics
    for (int k = 0; k < model->num_topics_; k++) {
      model->sorted_topics_[k] = (WordWeightPair*) malloc(sizeof(WordWeightPair) * (model->num_terms_));
    }
    model->SortWordProbEachTopic(model->eta_);

    // init mislead topics as original topics
    model->log_origin_topics_ = (double**) malloc(sizeof(double*)*(model->num_topics_));
    model->mislead_probw_ = (double**) malloc(sizeof(double*)*(model->num_topics_));
    for (int k = 0; k < model->num_topics_; k++) {
      model->log_origin_topics_[k] = (double*) malloc(sizeof(double) * (model->num_terms_));
      model->mislead_probw_[k] = (double*) malloc(sizeof(double) * (model->num_terms_));
      for(int v = 0; v < model->num_terms_; v++) {
        model->log_origin_topics_[k][v] = model->eta_[k][v];
        model->mislead_probw_[k][v] = model->eta_[k][v];
      }
    }

    //printf("original proportion for word %s : %0.5lf\n", corpus->dict[sorted_topics_[1][0].word_id], mislead_probw_[0][sorted_topics_[1][0].word_id]);
    // make the attack words at the desired rank in attack topic
    for(int i = 0; i < pm_att->vector_attack_wordid.size(); i++) {
      model->mislead_probw_[pm_att->vector_attack_topicid[i]][pm_att->vector_attack_wordid[i]] 
      = model->sorted_topics_[pm_att->vector_attack_topicid[i]][pm_att->vector_attack_word_sortid[i]].proportion;
    }
    for(int k = 0; k < model->num_topics_; k++) {
      NormalizeArray(model->mislead_probw_[k], model->num_terms_);
    }

    for(int i = 0; i < pm_att->vector_attack_wordid.size(); i++) {
      printf("mislead proportion for word %s in topic %d : %0.5lf\n", corpus->dict[pm_att->vector_attack_wordid[i]] 
             ,pm_att->vector_attack_topicid[i] 
             , model->mislead_probw_[pm_att->vector_attack_topicid[i]][pm_att->vector_attack_wordid[i]]);
    }

    PrintTopWordsInTopic(fout, model->eta_, 15, model, corpus, "original");
    PrintTopWordsInTopic(fout, model->mislead_probw_, 15, model, corpus, "mislead");
  }

void LossChangeAttackPair(Corpus* corpus, LdaModel* model, ParamAttack* pm_att, std::vector<double>& digamma_sumeta, std::vector<double>& phinew,
                          AttackPair* tmp_attack_p) {
                            // record the sum_prod value as gain in AttackPair
                            int d = tmp_attack_p->attack_doc;
                            int w = tmp_attack_p->attack_word;
                            // calculate the topic assignment probability of inserted new word (d,w)
                            for(int k = 0; k < model->num_topics_; k++) {
                              phinew[k] = exp(Digamma(model->gamma_[d][k])); // + std::max(Digamma(eta_[k][w]),20.0) - digamma_sumeta[k]);
                            }
                            // normalize the probability
                            NormalizeArray(phinew, model->num_topics_);

                            double delta_logprob = 0;
                            for(int k = 0; k < model->num_topics_; k++) {
                              //printf("T%d ori-mis=%0.5lf sumeta[k]=%0.5lf eta=%0.5lf phi=%0.5lf\n", k, eta_[k][w]/sumeta_[k] - mislead_probw_[k][w],
                              //	sumeta[k], eta_[k][w], phinew[k]);
                              delta_logprob -= (model->eta_[k][w]/model->sumeta_[k] - model->mislead_probw_[k][w]) 
                              *(model->sumeta_[k] - model->eta_[k][w])/(model->sumeta_[k]*model->sumeta_[k])
                              * phinew[k];
                            }
                            // if the attack does not use out of budget per document
                            if(corpus->docs[d].attack_total + fabs(pm_att->step_length_insert * delta_logprob) <= pm_att->limit_perdoc_attackwords) {
                              if(delta_logprob > 0)
                              tmp_attack_p->delta_w = pm_att->step_length_insert * delta_logprob;
                              else 
                              tmp_attack_p->delta_w = pm_att->step_length_delete * delta_logprob;
                            }
                            else
                            tmp_attack_p->delta_w = 0.0; // if the total number of attack words added exceeds the Limit, the gain is set as zero.
                          }

void AttackAllPairs(Corpus* corpus, LdaModel* model, 
                    ParamAttack* pm_att, std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q) {
                      std::vector<double> digamma_sumeta;
                      digamma_sumeta.resize(model->num_topics_);
                      std::vector<double> phinew;
                      phinew.resize(model->num_topics_);

                      model->CalcDigammaSumEta(digamma_sumeta);

                      std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* delete_qp 
                      = new std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>();
                      std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* insert_qp 
                      = new std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>();

                      AttackPair tmp_attack_p;
                      for(int w = 0; w < model->num_terms_; w++) {
                        if(model->DifferenceSpecificWordOriginMisleadTopic(w) > 1e-4) {
                          for(int d = 0; d < model->num_documents_; d++) {
                            tmp_attack_p.attack_doc = d;
                            tmp_attack_p.attack_word = w;
                            // calculat the delta of w
                            LossChangeAttackPair(corpus, model, pm_att, digamma_sumeta, phinew, &tmp_attack_p);

                            if(tmp_attack_p.delta_w > kEpsilon) {
                              // add current insert pair into priority queue (minimal heap)
                              insert_qp->push(tmp_attack_p); 
                            }
                            else if(tmp_attack_p.delta_w < -kEpsilon) {
                              // add current delete pair into priority queue (minimal heap)
                              delete_qp->push(tmp_attack_p);
                            }			
                          }
                        }
                        ResizeQueueAttackPairs(kLimitSizeQueue/2, delete_qp);
                        ResizeQueueAttackPairs(kLimitSizeQueue/2, insert_qp);
                      }

                      InsertQueueAttackPairs(insert_qp, att_q);
                      InsertQueueAttackPairs(delete_qp, att_q);
                      ResizeQueueAttackPairs(0, delete_qp);
                      ResizeQueueAttackPairs(0, insert_qp);
                    }

void RandomAttack(Corpus* corpus, ParamAttack* pm_att, 
                  std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q) {
                    // clear queue of attack pairs
                    ResizeQueueAttackPairs(0, att_q);
                    AttackPair tmp_attack_p;
                    // only select words that appear in attack pairs
                    std::vector<int> vec_docid;
                    for(int i = 0; i < pm_att->vector_attack_wordid.size(); i++)
                    {
                      vec_docid.clear();
                      int w =	pm_att->vector_attack_wordid[i];
                      int k = pm_att->vector_attack_topicid[i];
                      tmp_attack_p.attack_word = w;

                      for(int d = 0; d < corpus->num_docs; d++)
                      if(corpus->docs[d].attack_total + 1.0 <= pm_att->limit_perdoc_attackwords)
                      {
                        // want to increase the proportion of word w
                        if(pm_att->vector_expected_attack_effort[i] > 0)
                        {
                          vec_docid.push_back(d);
                        }
                        // want to decrease the proportion of word w
                        else 
                        {
                          bool exist_deleteword = false;
                          for(int n = 0; n < corpus->docs[d].length; n++)
                          if((*corpus->docs[d].words)[n] == w && (*corpus->docs[d].counts)[n] >= 1.0)
                          {
                            exist_deleteword = true;
                            break;
                          }
                          if(exist_deleteword)
                          vec_docid.push_back(d);
                        }
                      }

                      for(int j = 1; j < vec_docid.size(); j++)
                      {
                        int p = rand()%j;
                        std::swap(vec_docid[j], vec_docid[p]);
                      }
                      printf("value1 = %0.10lf   value2 = %0.10lf\n", (double)vec_docid.size(), fabs(pm_att->vector_expected_attack_effort[i])/(pm_att->limit_attack_iterations-1));
                      for(int j = 0; j < std::min((double)vec_docid.size(), fabs(pm_att->vector_expected_attack_effort[i])/(pm_att->limit_attack_iterations-1)); j++)
                      {
                        tmp_attack_p.attack_doc = vec_docid[j];
                        if(pm_att->vector_expected_attack_effort[i] > 0)
                        {
                          tmp_attack_p.delta_w = 1.0;
                        }
                        // want to decrease the proportion of word w
                        else 
                        {
                          tmp_attack_p.delta_w = -1.0;
                        }
                        att_q->push(tmp_attack_p);
                      }
                      ResizeQueueAttackPairs(kLimitSizeQueue, att_q);	
                    }
}


void AttackPairsInQueue(Corpus* corpus, LdaModel* model, 
                        ParamAttack* pm_att, std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q)
{
  std::vector<double> digamma_sumeta;
  digamma_sumeta.resize(model->num_topics_);
  std::vector<double> phinew;
  phinew.resize(model->num_topics_);

  model->CalcDigammaSumEta(digamma_sumeta);

  std::vector<AttackPair> tmp_vec_attpair;
  tmp_vec_attpair.resize(att_q->size());
  int id = 0;
  // pop all elements to the temp vector of attack pairs
  while(!att_q->empty())
  {
    tmp_vec_attpair[id++] = att_q->top();
    att_q->pop();
  }

  for(int i = 0; i < tmp_vec_attpair.size();i++)
  {
    // calculat the new delta of w
    LossChangeAttackPair(corpus, model, pm_att, digamma_sumeta, phinew,  &tmp_vec_attpair[i]);
    // add current attack pair into priority queue (minimal heap)
    att_q->push(tmp_vec_attpair[i]);
  }
  tmp_vec_attpair.clear();
}


void ModifyCorpus(AttackPair attpair, Corpus* corpus)
{
  // real delta is the actuall number of added/removed words from W
  double real_delta = 0;
  int d = attpair.attack_doc;
  int w = attpair.attack_word;
  if(attpair.delta_w > 0)
  {
    real_delta = attpair.delta_w;
    // flag : if the insert word exist in doc.words
    bool flag = false;
    for(int i = 0; i < corpus->docs[d].words->size(); i++)
    if((*corpus->docs[d].words)[i] == attpair.attack_word)
    {
      (*corpus->docs[d].counts)[i] += real_delta;
      flag = true;
    }
    if(!flag)
    {
      corpus->docs[d].words->push_back(w);
      corpus->docs[d].counts->push_back(attpair.delta_w);
    }
  }
  else if(attpair.delta_w < 0)
  {
    for(int i = 0; i < corpus->docs[d].words->size(); i++)
    if((*corpus->docs[d].words)[i] == w)
    while(real_delta > attpair.delta_w && (*(corpus->docs[d].counts))[i] > 0)
    {
      double tmp_delta = std::min((*(corpus->docs[d].counts))[i], real_delta - attpair.delta_w);
      (*(corpus->docs[d].counts))[i] -= tmp_delta;
      real_delta -= tmp_delta;
    }
  }

  corpus->docs[d].length = corpus->docs[attpair.attack_doc].words->size();
  corpus->docs[d].total += real_delta;
  corpus->docs[d].attack_total += fabs(real_delta);
}

void EstimateGradientStepLength(Corpus* corpus, LdaModel* model, std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* att_q 
                                , ParamAttack* paramattack)
{
  CopyToMaxQueueFromMinQueue(att_q, &max_q);
  double effort_periter_delete = 0.0;
  double effort_periter_insert = 0.0;
  // calculate the delta of delete/insert per iteration
  while(!max_q.empty())
  {
    if(max_q.top().delta_w > 0) 
    effort_periter_insert += fabs(max_q.top().delta_w);
    else
    effort_periter_delete += fabs(max_q.top().delta_w);
    max_q.pop();
  }

  double sum_effort_delete = 0.0;
  double sum_effort_insert = 0.0;
  for(int i = 0; i < paramattack->vector_attack_wordid.size(); i++)
  {
    int k = paramattack->vector_attack_topicid[i];
    int w = paramattack->vector_attack_wordid[i];
    if(model->eta_[k][w] < model->mislead_probw_[k][w]*model->sumeta_[k])
    sum_effort_insert += fabs(model->eta_[k][w] - model->mislead_probw_[k][w]*model->sumeta_[k]);
    else
    sum_effort_delete += fabs(model->eta_[k][w] - model->mislead_probw_[k][w]*model->sumeta_[k]);
  }

  // set step length as the total effort needed divide the some of effort of the first iteration
  paramattack->step_length_delete = 16.0 * sum_effort_delete / (effort_periter_delete * paramattack->limit_attack_iterations);
  paramattack->step_length_insert = 8.0 * sum_effort_insert / (effort_periter_insert * paramattack->limit_attack_iterations);
}

void AttackLda(int argc, char* argv[]) {
   // (est / inf) alpha k settings data (random / seed/ model) (directory / out) (#words for attacking) (#iterations)

  Corpus* corpus;
  double t1 = atof(argv[11]);
  seedMT(t1);
  paramattack.step_length_delete = 1.0;
  paramattack.step_length_insert = 1.0;

  paramattack.initial_alpha_value = atof(argv[2]);
  paramattack.initial_beta_value = atof(argv[8]);
  paramattack.num_topics = atoi(argv[3]);
  ReadAttackSettings(argv[4], &paramattack);
  corpus = ReadCorpus(argv[5], argv[9]);
  MakeDirectory(argv[7]);
  char* outinfo_filename = argv[10];
  char output_selecteddocs_fn[200];
  strcpy(output_selecteddocs_fn, outinfo_filename);
  strcat(output_selecteddocs_fn, "_selectedDocs");
  bool is_variational_attack = true;
  FILE* fout = fopen(outinfo_filename, "ab");

  if(strcmp(argv[12], "random") == 0) {
    printf("Random Strategy\n");
    is_variational_attack = false;
    ReadEffortSettings(argv[13], &paramattack);
    RandomAttack(corpus, &paramattack, &pq);
    for(int i = 0; i < paramattack.vector_expected_attack_effort.size(); i++)
    fprintf(fout, "effort_word %s : %0.10lf\n", paramattack.vector_attack_wordstr[i].c_str(), paramattack.vector_expected_attack_effort[i]);
  }

  LdaModel* model = NULL;
  LdaSuffstats* ss = NULL;
  InitialInfo init_info;
  init_info.InitializeLdaModel(paramattack.num_attack_words, argv[6], argv[7], corpus, &paramattack);
  model = init_info.model;
  ss = init_info.ss;

  model->RunVariationalInference(argv[7], corpus, &paramattack, ss);
 
  printf("finish EM\n");
  paramattack.vector_attack_wordid.resize(paramattack.vector_attack_topictheme.size());
  paramattack.vector_attack_topicid.resize(paramattack.vector_attack_topictheme.size());
  for(int i = 0; i < paramattack.vector_attack_topictheme.size(); i++) {
    paramattack.vector_attack_wordid[i] = WordIdDictionary(paramattack.vector_attack_wordstr[i], corpus->dict, model->num_terms_);
    paramattack.vector_attack_topicid[i] = FindAttackTopicId(paramattack.vector_attack_topictheme[i], corpus->dict, model->num_terms_, model->num_topics_, model->eta_);
  }
  printf("finish ATTACKword\n");

  PrintAttackSettings(fout, &paramattack);
  ComputeMisleadTopics(corpus, &paramattack, fout, model);
  PrintAttackMeasure(fout, corpus, model, &paramattack);
  if(is_variational_attack) {
    AttackAllPairs(corpus, model , &paramattack, &pq);
    EstimateGradientStepLength(corpus, model, &pq, &paramattack);
    fprintf(fout, "VariAttack step_length_delete = %0.5lf\n", paramattack.step_length_delete);
    fprintf(fout, "VariAttack step_length_insert = %0.5lf\n", paramattack.step_length_insert);
  }

  // automatically set the step_length
  fclose(fout);

  int permute[model->num_topics_]; // permute[i] the corresponding topic for i-th mislead topic;
  
  for(int iter = 0; iter < paramattack.limit_attack_iterations; iter++) {
    FILE* fout = fopen(outinfo_filename, "ab");
    FILE* fout_selectdocs = fopen(output_selecteddocs_fn, "ab");

    if(is_variational_attack) {
		AttackPairsInQueue(corpus, model, &paramattack, &pq);
	}
	else {
      RandomAttack(corpus, &paramattack, &pq);
      fprintf(fout, "queue size = %d\n", pq.size());	
    }

    CopyToMaxQueueFromMinQueue(&pq, &max_q);
    while(!max_q.empty()) {
      ModifyCorpus(max_q.top(), corpus);
      fprintf(fout_selectdocs, "(word=%s,doc=%d) delta=%0.10lf : ", 
              corpus->dict[max_q.top().attack_word], max_q.top().attack_doc, max_q.top().delta_w);
      PrintTopicProportions(fout_selectdocs, max_q.top().attack_doc,  
                            model->gamma_[max_q.top().attack_doc], 
                            model->num_topics_);
      max_q.pop();
    }

    // allocate a new model will cause segmentation fault after RunVariationalInference(), I do not know why:(
    // instead, we just initialize the topics as below
    for (int k = 0; k < model->num_topics_; k++) {
      for (int w = 0; w < model->num_terms_; w++) {
        model->eta_[k][w] = model->init_eta_[k][w];  // make the initialization the same       
      }
    }				

    CopyToMaxQueueFromMinQueue(&pq, &max_q);
    while(!max_q.empty()) {
      fprintf(fout_selectdocs, "(after attack) docs topic proportion (%d,%s): ", 
              max_q.top().attack_doc, corpus->dict[max_q.top().attack_word]);
      PrintDocumentWordStr(fout_selectdocs, max_q.top().attack_doc, 
                           corpus->docs[max_q.top().attack_doc],
                           corpus->dict);
      max_q.pop();
    }

//	printf("I am here\n");
	model->gibbs_sampling_.RunGibbsSampling(corpus, model->num_topics_, 1000, model->alpha_, model->beta_, model->eta_);
//	printf("I am here1\n");
	
    // find best permutation of the learned topic after attack
    model->FindBestTopicPermutation(permute);
	PrintTopicPermutation(fout, permute, model->num_topics_);
	model->PermuteMisleadTopics(permute);
	PrintTopWordsInTopic(fout, model->eta_, 15, model, corpus, "permuted");
    PrintAttackMeasure(fout, corpus, model, &paramattack);

	// allocate a new model will cause segmentation fault after RunVariationalInference(), I do not know why:(
    // instead, we just initialize the topics as below
    for (int k = 0; k < model->num_topics_; k++) {
      for (int w = 0; w < model->num_terms_; w++) {
        model->eta_[k][w] = model->init_eta_[k][w];  // make the initialization the same       
      }
    }				
	
	model->RunVariationalInference(argv[7], corpus, &paramattack, ss);

    fclose(fout);
    fclose(fout_selectdocs);
    }
}
} // namespace attack_lda


