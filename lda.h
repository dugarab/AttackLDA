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

#ifndef ATTACKLDA_LDA_H_
#define ATTACKLDA_LDA_H_
#include <vector>
#include <string>

namespace attack_lda {
typedef struct
{
    std::vector<int>* words;
    std::vector<double>* counts;
    int length;
    double total;
	double attack_total; // number of attack words added.
} Document;

typedef struct
{
    Document* docs;
	char** dict;
    int num_terms;
    int num_docs;
} Corpus;

typedef struct
{
	double proportion;
	int word_id;
}WordWeightPair;

typedef struct
{
    double** class_word;
    double* class_total;
    double alpha_suffstats;
    int num_docs;
} LdaSuffstats;

typedef struct
{
	int lag;

float e_step_converge_threshold;
int e_step_limit_iterations;
float em_converge_threshold;
int em_limit_iterations;
int if_estimate_alpha;
double initial_alpha_value;
double initial_beta_value;
int num_topics;

	int num_attack_words;
	int limit_attack_iterations;
	std::vector<int> vector_attack_wordid;
	std::vector<int> vector_attack_topicid;
	int limit_perdoc_attackwords;
	std::vector<std::string> vector_attack_wordstr;
	std::vector<std::string> vector_attack_topictheme;
	std::vector<int> vector_attack_word_sortid;

	// the expected attack effort for each attack pair
	// It is ONLY used in random attack.
	std::vector<double> vector_expected_attack_effort;
double step_length_delete;
double step_length_insert;
} ParamAttack;

} // namespace attack_lda
#endif
