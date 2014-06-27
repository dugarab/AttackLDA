#ifndef ATTACKLDA_UTILS_H_
#define ATTACKLDA_UTILS_H_

#include <cstdio>
#include <cmath>
#include <float.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)
#define NUM_INIT 1

double LogSum(double log_a, double log_b);
double Trigamma(double x);
double Digamma(double x);
double LogGamma(double x);
void MakeDirectory(char* name);
int IndexOfMax(double* x, int n);
void NormalizeArray(double* x , int len);
void NormalizeArray(std::vector<double>& x , int len);
double NormalizedL2Distance( double* a , double* b , int len );
double NumOverlapTopWordTopics(double* a , double* b , int len , int num_topwords );
double ProbOnRankInTopic( double* a , int len , int num_topwords );
void PrintTopicProportions( FILE* fout ,int doc_id ,  double* gamma , int len );
void PrintTopicPermutation( FILE* fout , int* permute , int len );
double CalcTopicL2Distance( double** a , double** b , int num_topics , int num_terms );
double SumVector( double* x, int len);
int WordIdDictionary( std::string word_str , char** dict , int len_dict );
int FindAttackTopicId( std::string topic_theme , char** dict , int num_terms , int num_topics , double** eta );
double EpsilonInsensitiveSquare(double x);
int WordRankInTopic( int num_terms , int topic_id , int word_id , double** eta );
#endif
