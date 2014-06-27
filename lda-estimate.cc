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

#include "lda-estimate.h"

namespace attack_lda {

  void InitialInfo::InitializeLdaModel(int num_attackwords, char* start, char* directory, Corpus* corpus, ParamAttack *pm_att) {
    // initialize model

    char filename[100];
    if (strcmp(start, "seeded")==0)
    {
	  model = NewLdaModel(pm_att->num_topics, corpus, num_attackwords, pm_att->initial_alpha_value, pm_att->initial_beta_value);
      ss = NewLdaSuffstats(model);
      InitializeSuffstatsCorpus(corpus);
      model->LdaMle(ss);
	}
    else if (strcmp(start, "random")==0)
    {
      printf("start here");
	  model = NewLdaModel(pm_att->num_topics, corpus, num_attackwords, pm_att->initial_alpha_value, pm_att->initial_beta_value);
      ss = NewLdaSuffstats(model);
      InitializeSuffstatsRandom();
      model->LdaMle(ss);
      printf("end here");
	}
    else
    {
      model = LoadLdaModel(start, corpus, num_attackwords);
      ss = NewLdaSuffstats(model);
    }

    sprintf(filename,"%s/000",directory);
    model->SaveLdaModel(filename);
  }

/*
* allocate sufficient statistics
*
*/
LdaSuffstats* InitialInfo::NewLdaSuffstats(LdaModel* model) {
  int num_topics = model->num_topics_;
  int num_terms = model->num_terms_;
  int i,j;

  LdaSuffstats* ss = new LdaSuffstats();
  ss->class_total = (double*) malloc(sizeof(double)*num_topics);
  ss->class_word = (double**) malloc(sizeof(double*)*num_topics);
  for (i = 0; i < num_topics; i++) {
    ss->class_total[i] = 0;
    ss->class_word[i] = (double*) malloc(sizeof(double)*num_terms);
    for (j = 0; j < num_terms; j++) {
      ss->class_word[i][j] = 0;
    }
  }
  return(ss);
}

  void InitialInfo::InitializeSuffstatsCorpus(Corpus* c)
  {
    int num_topics = model->num_topics_;
    int i, k, d, n;
    Document* doc;

    for (k = 0; k < num_topics; k++)
    {
      for (i = 0; i < NUM_INIT; i++)
      {
        d = floor(myrand() * c->num_docs);
        printf("initialized with document %d\n", d);
        doc = &(c->docs[d]);
        for (n = 0; n < doc->length; n++)
        {
          ss->class_word[k][(*(doc->words))[n]] += (*(doc->counts))[n];
        }
      }
      for (n = 0; n < model->num_terms_; n++)
      {
        ss->class_word[k][n] += 1.0;
        ss->class_total[k] = ss->class_total[k] + ss->class_word[k][n];
      }
    }
  }


void InitialInfo::InitializeSuffstatsRandom() {
  int num_topics = model->num_topics_;
  int num_terms = model->num_terms_;
  for (int k = 0; k < num_topics; k++) {
    for (int n = 0; n < num_terms; n++) {
      ss->class_word[k][n] += 1.0/num_terms + myrand();
      ss->class_total[k] += ss->class_word[k][n];
    }
  }
}

/*
* allocate new lda model
*
*/
LdaModel* InitialInfo::NewLdaModel(int num_topics, Corpus* corpus, int num_attackwords, double alpha, double beta) {
  int i,j;
  LdaModel* model;

  model = new LdaModel();
  model->num_topics_ = num_topics;
  model->num_terms_ = corpus->num_terms;
  model->num_documents_ = corpus->num_docs;

  model->alpha_ = alpha;
  model->beta_ = beta;
  model->eta_ = (double**) malloc(sizeof(double*)*num_topics);
  model->sumeta_ = (double*) malloc(sizeof(double)*num_topics);
  model->init_eta_ = (double**) malloc(sizeof(double*)*num_topics);
 
  model->gamma_ = (double**) malloc(sizeof(double*)*(corpus->num_docs));
  for (int d = 0; d < corpus->num_docs; d++) {
    model->gamma_[d] = (double*) malloc(sizeof(double) * num_topics);
  }

  int max_length = MaxDocumentLength(corpus);
  model->phi_ = (double**) malloc(sizeof(double*)*(max_length + num_attackwords));
  for (int n = 0; n < (max_length+num_attackwords); n++) {
    model->phi_[n] = (double*) malloc(sizeof(double) * num_topics);    
  }

  for (i = 0; i < num_topics; i++) {
    model->eta_[i] = (double*) malloc(sizeof(double)*model->num_terms_);
    model->init_eta_[i] = (double*) malloc(sizeof(double)*model->num_terms_);
    for (j = 0; j < model->num_terms_; j++) {
      model->eta_[i][j] = 0;
      model->init_eta_[i][j] = 0;
    }
  }

  model->sorted_topics_ = (WordWeightPair**) malloc(sizeof(WordWeightPair*)*(model->num_topics_));
  // sort the learned topics
  for (int k = 0; k < model->num_topics_; k++) {
    model->sorted_topics_[k] = (WordWeightPair*) malloc(sizeof(WordWeightPair) * (model->num_terms_));
  }

  return(model);
}

LdaModel* InitialInfo::LoadLdaModel(char* model_root, Corpus* corpus, int num_attackwords) {
  char filename[100];
  FILE* fileptr;
  int i, j, num_terms, num_topics;
  float x, alpha;

  sprintf(filename, "%s.other", model_root);
  printf("loading %s\n", filename);
  fileptr = fopen(filename, "r");
  fscanf(fileptr, "num_topics %d\n", &num_topics);
  fscanf(fileptr, "num_terms %d\n", &num_terms);
  fscanf(fileptr, "alpha %f\n", &alpha);
  fclose(fileptr);

  LdaModel* model = NewLdaModel(num_topics, corpus, num_attackwords, 0, 0);
  
  sprintf(filename, "%s.beta", model_root);
  printf("loading %s\n", filename);
  fileptr = fopen(filename, "r");
  for (i = 0; i < num_topics; i++) {
    for (j = 0; j < num_terms; j++) {
      fscanf(fileptr, "%f", &x);
      model->eta_[i][j] = x;
    }
  }
  fclose(fileptr);
  return(model);
}

} // namespace attack_lda
