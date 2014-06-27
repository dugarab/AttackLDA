#include "util_topics.h"

namespace attack_lda {
  void ResizeQueueAttackPairs(int lower_size, std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* qp) {
    while(qp->size() > lower_size) {
      qp->pop();
    }
  }

  void InsertQueueAttackPairs(std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* q1, 
                              std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* q2) {
                                while(!q1->empty()) {
                                  q2->push(q1->top());
                                  q1->pop();
                                }
                              }

  void CopyToMaxQueueFromMinQueue(std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsLessthanByGain>* min_q , 
                                  std::priority_queue<AttackPair, std::vector<AttackPair>, AttackPairsMorethanByGain>* max_q) {
                                    std::vector<AttackPair> tmp_vec;
                                    tmp_vec.resize(min_q->size());
                                    for(int i = 0; i < tmp_vec.size(); i++) {
                                      tmp_vec[i] = min_q->top();
                                      min_q->pop();
                                    }
                                    // clear the maximal heap
                                    while(!max_q->empty()) max_q->pop();

                                    for(int i = 0; i < tmp_vec.size(); i++) {
                                      max_q->push(tmp_vec[i]);
                                      min_q->push(tmp_vec[i]);
                                    }
                                    tmp_vec.clear();
                                  }

  int CompareWordWeightPair(const void *a, const void *b) {
    if(((WordWeightPair*)a)->proportion > ((WordWeightPair*)b)->proportion) {
      return -1;
    }
    else if (((WordWeightPair*)a)->proportion < ((WordWeightPair*)b)->proportion) {
      return 1;
    }
    else {
      return 0;
    }
  }

  void LdaModel::FindBestTopicPermutation(int* permute) {
    int used[num_topics_];
    int k = 0 , i = 0;
    for(k = 0; k < num_topics_; k++) {
      used[k] = 0;
    }

    for(k = 0; k < num_topics_; k++) {
      double min_dis = 1e9;
      int select_id = 0;
      for(i = 0; i < num_topics_; i++) {
        if(used[i] == 0) {
          if(NumOverlapTopWordTopics(log_origin_topics_[k] , eta_[i] , num_terms_ , 10) < min_dis) {
            min_dis = NumOverlapTopWordTopics(log_origin_topics_[k] , eta_[i] , num_terms_ , 10);
            select_id = i;
          }
        }
      }

      used[select_id] = 1;
      permute[k] = select_id;
    }
  }

void LdaModel::PermuteMisleadTopics(int* permute) {
  double tmp[num_topics_];
  int k = 0 , v = 0 , d = 0;
  for(v = 0; v < num_terms_; v++) {
    for(k = 0; k < num_topics_; k++) {
      tmp[k] = eta_[permute[k]][v];
    }
    for(k = 0; k < num_topics_; k++) {
      eta_[k][v] = tmp[k];
    }
  }

  for(d = 0; d < num_documents_; d++) {
    for(k = 0; k < num_topics_; k++) {
      tmp[k] = gamma_[d][permute[k]];
    }
    for(k = 0; k < num_topics_; k++) {
      gamma_[d][k] = tmp[k];
    }
  }
}

double LdaModel::L2DistanceMisleadEstimateTopics(FILE* fout, ParamAttack* pm_att) {
  int k = 0 , v = 0;
  double res = 0.0;

  // use lgamma() not LogGamma(), LogGamma(x) has problem for small positive x
  for(int i = 0; i < pm_att->vector_attack_wordid.size(); i++) {
    int k = pm_att->vector_attack_topicid[i];
    int w = pm_att->vector_attack_wordid[i];
    double sum_eta = 0;
    for(v = 0; v < num_terms_; v++) {
      sum_eta += eta_[k][v];
    }

    //for(v = 0; v < num_terms_; v++)
    //{
    res += EpsilonInsensitiveSquare(eta_[k][w]/sum_eta - mislead_probw_[k][w]);
    if(EpsilonInsensitiveSquare(eta_[k][w]/sum_eta - mislead_probw_[k][w]) > 1e-8) {
      fprintf(fout , "T = %d , w = %d  eta/sum = %0.10lf  misprob = %0.10lf\n" , k , w , eta_[k][w]/sum_eta , mislead_probw_[k][w]);
    }
    //}
  }
  return res;
}

void LdaModel::SortWordProbEachTopic(double** log_beta) {
  int k = 0 , v = 0;
  for (k = 0; k < num_topics_; k++) {
    for(v = 0; v < num_terms_; v++) {
      sorted_topics_[k][v].proportion = log_beta[k][v];
      sorted_topics_[k][v].word_id = v;
    }
    qsort(sorted_topics_[k] , num_terms_ , sizeof(WordWeightPair), CompareWordWeightPair);
  }
}

double LdaModel::DifferenceSpecificWordOriginMisleadTopic(int w) {
  double diff = 0.0;
  for(int k = 0; k < num_topics_; k++) {
    diff += fabs(eta_[k][w]/sumeta_[k] - mislead_probw_[k][w]);
  }
  return diff;
}

/*
* run_em
*
*/
void LdaModel::RunVariationalInference(char* directory, Corpus* corpus , ParamAttack *pm_att , LdaSuffstats* ss ) {
  int i = 0 , d = 0;
  double likelihood, likelihood_old = 0, converged = 1;
  //char filename[100];
  //sprintf(filename, "%s/likelihood.dat", directory);

  //FILE* likelihood_file = fopen(filename, "w");
  std::vector<double> digamma_sumeta;
  digamma_sumeta.resize(num_topics_);

  while (((converged < 0) || (converged > pm_att->em_converge_threshold) || (i <= 2)) && (i <= pm_att->em_limit_iterations)) {
    i++; printf("**** em iteration %d ****\n", i);
    likelihood = 0;
    InitializeSuffstatsZero(ss);
    // e-step
    CalcDigammaSumEta( digamma_sumeta );
    double vari_prob_topics = ComputeVariationalProb();

    printf( "vari prob for topics = %0.5lf\n" , vari_prob_topics );

    for (d = 0; d < corpus->num_docs; d++) {
      // if ((d % 1000) == 0) printf("document %d\n",d);
      likelihood += ExpectationStepOneDoc(&(corpus->docs[d]),
                                          d,
                                          digamma_sumeta,
                                          pm_att,
                                          ss);
    }
    likelihood += vari_prob_topics;
    // m-step
    LdaMle(ss);
    printf("finish-m-step\n")
    ;
    // check for convergence
    converged = (likelihood_old - likelihood) / (likelihood_old);
    if (converged < 0) pm_att->e_step_limit_iterations = pm_att->e_step_limit_iterations * 2;
    likelihood_old = likelihood;

    // output model and likelihood

    /*fprintf(likelihood_file, "%10.10f\t%5.5e\n", likelihood, converged);
    fflush(likelihood_file);
    if ((i % LAG) == 0)
    {
    sprintf(filename,"%s/%03d",directory, i);
    SaveLdaModel(model, filename);
    sprintf(filename,"%s/%03d.gamma",directory, i);
    SaveGamma(filename, var_gamma, corpus->num_docs, num_topics_);
  }
    */
  }
}

/*
* various intializations for the sufficient statistics
*/
void LdaModel::InitializeSuffstatsZero(LdaSuffstats* ss) {
  int k, w;
  for (k = 0; k < num_topics_; k++) {
    ss->class_total[k] = 0;
    for (w = 0; w < num_terms_; w++) {
      ss->class_word[k][w] = 0;
    }
  }
  ss->num_docs = 0;
  ss->alpha_suffstats = 0;
}

void LdaModel::CalcDigammaSumEta(std::vector<double>& digamma_sumeta) {
  for(int k = 0; k < num_topics_; k++) {
    digamma_sumeta[k] = 0;
    for(int v = 0; v < num_terms_; v++)
    digamma_sumeta[k] += eta_[k][v];
    digamma_sumeta[k] = Digamma(digamma_sumeta[k]);
  }
}

/*
* compute the E_{q(topics)}[log p(topics|beta)] - E_{q(topics)}[log q(topics|eta)]
*/
double LdaModel::ComputeVariationalProb() {
  int k = 0 , v = 0;
  double res_prob = 0.0;
  std::vector<double> digamma_sumeta;
  digamma_sumeta.resize(num_topics_);
  CalcDigammaSumEta(digamma_sumeta );

  // use lgamma() not LogGamma(), LogGamma(x) has problem for small positive x
  for(k = 0; k < num_topics_; k++) {
    res_prob += lgamma(beta_ * num_terms_) - num_terms_ * lgamma(beta_);
    double sum_eta = 0;
    //		printf( "res_prob = %0.5lf\n" , res_prob );
    //		printf( "digamma_sumeta = %0.5lf\n" , digamma_sumeta[k] );

    double tmp = 0;
    for(v = 0; v < num_terms_; v++) {
      sum_eta += eta_[k][v];
      tmp += lgamma(eta_[k][v]) + (beta_ - eta_[k][v]) * (Digamma(eta_[k][v]) - digamma_sumeta[k]);
      //			printf( "%0.5lf  %0.5lf  %0.5lf\n" , lgamma(eta_[k][v]) , beta_ - eta_[k][v] , Digamma(eta_[k][v]) - digamma_sumeta[k] );
    }
    //		printf( "tmp = %0.5lf\n" , tmp );
    //		printf( "sum_eta = %0.5lf\n" , sum_eta);
    res_prob -= lgamma(sum_eta);
    res_prob += tmp;

    //	printf( "res_prob_full = %0.5lf\n" , res_prob );
  }

  return (res_prob);
}

/*
* compute MLE lda model from sufficient statistics
*
*/
void LdaModel::LdaMle(LdaSuffstats* ss) {
  for (int k = 0; k < num_topics_; k++) {
    for (int w = 0; w < num_terms_; w++) {
      eta_[k][w] = beta_ + ss->class_word[k][w];        
      // record the initialization value of eta
	  init_eta_[k][w] = eta_[k][w];
    }
  }
  for(int k = 0; k < num_topics_; k++) {
	sumeta_[k] = SumVector( eta_[k] , num_terms_ );
  }
}

/*
* perform inference on a document and update sufficient statistics
*
*/
double LdaModel::ExpectationStepOneDoc(Document* doc, int d, std::vector<double>& digamma_sumeta, ParamAttack* pm_att, LdaSuffstats* ss) {
  double likelihood;
  int n, k;

  // posterior inference
  likelihood = LdaInference(doc, d , pm_att, digamma_sumeta );

  // update sufficient statistics
  double gamma_sum = 0;
  for(k = 0; k < num_topics_; k++) {
    gamma_sum += gamma_[d][k];
    ss->alpha_suffstats += Digamma(gamma_[d][k]);
  }
  ss->alpha_suffstats -= num_topics_ * Digamma(gamma_sum);

  for(n = 0; n < doc->length; n++) {
    for(k = 0; k < num_topics_; k++) {
      ss->class_word[k][(*(doc->words))[n]] += (*doc->counts)[n]*phi_[n][k];
      ss->class_total[k] += (*(doc->counts))[n]*phi_[n][k];
    }
  }
  ss->num_docs = ss->num_docs + 1;
  return(likelihood);
}

/*
* variational inference
*
*/
double LdaModel::LdaInference(Document* doc, int d , ParamAttack* pm_att, std::vector<double>& digamma_sumeta) {
  double converged = 1;
  double phisum = 0, likelihood = 0;
  double likelihood_old = 0, oldphi[num_topics_];
  int k, n , var_iter;
  double digamma_gam[num_topics_];

  // compute posterior dirichlet

  for (k = 0; k < num_topics_; k++) {
    gamma_[d][k] = alpha_ + (doc->total/((double) num_topics_));
    digamma_gam[k] = Digamma(gamma_[d][k]);
    for (n = 0; n < doc->length; n++)
    phi_[n][k] = 1.0/num_topics_;
  }
  var_iter = 0;

  while ((converged > pm_att->e_step_converge_threshold) &&
         ((var_iter < pm_att->e_step_limit_iterations) || (pm_att->e_step_limit_iterations == -1))) {
           var_iter++;
           for (n = 0; n < doc->length; n++) {
             phisum = 0;
             for (k = 0; k < num_topics_; k++) {
               oldphi[k] = phi_[n][k];
               phi_[n][k] =
               digamma_gam[k] +
               Digamma(eta_[k][(*(doc->words))[n]]) -
               digamma_sumeta[k];

               if (k > 0)
               phisum = LogSum(phisum, phi_[n][k]);
               else
               phisum = phi_[n][k]; // note, phi_ is in log space
             }

             for (k = 0; k < num_topics_; k++) {
               phi_[n][k] = exp(phi_[n][k] - phisum);
               gamma_[d][k] =
               gamma_[d][k] + (*doc->counts)[n]*(phi_[n][k] - oldphi[k]);
               // !!! a lot of extra digamma's here because of how we're computing it
               // !!! but its more automatically updated too.
               digamma_gam[k] = Digamma(gamma_[d][k]);
             }
           }

           likelihood = ComputeLikelihood(doc, d, digamma_sumeta);
           //assert(!isnan(likelihood));
           converged = (likelihood_old - likelihood) / likelihood_old;
           likelihood_old = likelihood;
           // printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);
         }
return(likelihood);
}

/*
* compute likelihood bound
*
*/
double LdaModel::ComputeLikelihood(Document* doc, int d, std::vector<double>& digamma_sumeta ) {
  double likelihood = 0, digsum = 0, var_gamma_sum = 0, dig[num_topics_];

  for (int k = 0; k < num_topics_; k++) {
    dig[k] = Digamma(gamma_[d][k]);
    var_gamma_sum += gamma_[d][k];
  }
  digsum = Digamma(var_gamma_sum);

  likelihood = lgamma(alpha_ * num_topics_) - num_topics_ * lgamma(alpha_)
  - (lgamma(var_gamma_sum));

  for (int k = 0; k < num_topics_; k++) {
    likelihood += (alpha_ - 1)*(dig[k] - digsum) + lgamma(gamma_[d][k])
    - (gamma_[d][k] - 1)*(dig[k] - digsum);
    for (int n = 0; n < doc->length; n++) {
      if (phi_[n][k] > 0) {
        likelihood += (*doc->counts)[n]*
        (phi_[n][k]*((dig[k] - digsum) - log(phi_[n][k])
                     +  Digamma(eta_[k][(*(doc->words))[n]]) - digamma_sumeta[k] ));
      }
    }
  }
  return(likelihood);
}

/*
* deallocate new lda model
*
*/
void LdaModel::FreeLdaModel() {
  for (int i = 0; i < num_topics_; i++) {
    free(eta_[i]);
  }
  free(eta_);
}

/*
* save an lda model
*
*/
void LdaModel::SaveLdaModel(char* model_root) {
  char filename[100];
  FILE* fileptr;

  sprintf(filename, "%s.beta", model_root);
  fileptr = fopen(filename, "w");
  for (int i = 0; i < num_topics_; i++) {
    for (int j = 0; j < num_terms_; j++) {
      fprintf(fileptr, " %5.10f", eta_[i][j]);
    }
    fprintf(fileptr, "\n");
  }
  fclose(fileptr);

  sprintf(filename, "%s.other", model_root);
  fileptr = fopen(filename, "w");
  fprintf(fileptr, "num_topics %d\n", num_topics_);
  fprintf(fileptr, "num_terms %d\n", num_terms_);
  fprintf(fileptr, "alpha %5.10f\n", alpha_);
  fclose(fileptr);
}
} // namespace attack_lda
