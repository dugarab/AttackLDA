#include "utils.h"

/*
 * given log(a) and log(b), return log(a + b)
 *
 */

double LogSum(double log_a, double log_b)
{
  double v;

  if (log_a < log_b)
  {
      v = log_b+log(1 + exp(log_a-log_b));
  }
  else
  {
      v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}

 /**
   * Proc to calculate the value of the trigamma, the second
   * derivative of the loggamma function. Accepts positive matrices.
   * From Abromowitz and Stegun.  Uses formulas 6.4.11 and 6.4.12 with
   * recurrence formula 6.4.6.  Each requires workspace at least 5
   * times the size of X.
   *
   **/

double Trigamma(double x)
{
    double p;
    int i;

    x=x+6;
    p=1/(x*x);
    p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
    for (i=0; i<6 ;i++)
    {
        x=x-1;
        p=1/(x*x)+p;
    }
    return(p);
}


/*
 * taylor approximation of first derivative of the log gamma function
 *
 */

double Digamma(double x)
{
	if( x < 0.1 ) x = 0.1;

    double p;
    x=x+6;
    p=1/(x*x);
    p=(((0.004166666666667*p-0.003968253986254)*p+
	0.008333333333333)*p-0.083333333333333)*p;
    p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
    return p;
}

// LogGamma(x) gets extremely large value for small positive x
// use lgamma(x) instead
double LogGamma(double x)
{
	double z=1/(x*x);

    x=x+6;
    z=(((-0.000595238095238*z+0.000793650793651)
	*z-0.002777777777778)*z+0.083333333333333)/x;
    z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
	log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
    return z;
}



/*
 * make directory
 *
 */

void MakeDirectory(char* name)
{
    mkdir(name, S_IRUSR|S_IWUSR|S_IXUSR);
}


/*
 * argmax
 *
 */

int IndexOfMax(double* x, int n)
{
    int i;
    double max = x[0];
    int argmax = 0;
    for (i = 1; i < n; i++)
    {
        if (x[i] > max)
        {
            max = x[i];
            argmax = i;
        }
    }
    return(argmax);
}

void NormalizeArray(std::vector<double>& x , int len)
{
	double sum = 0;
	int i = 0;
	for(i = 0; i < len; i++)
		sum += x[i];
	for(i = 0; i < len; i++)
		x[i] /= sum;
}

void NormalizeArray(double* x , int len)
{
	double sum = 0;
	int i = 0;
	for(i = 0; i < len; i++)
		sum += x[i];
	for(i = 0; i < len; i++)
		x[i] /= sum;
}

double NormalizedL2Distance( double* a , double* b , int len )
{
	double res = 0.0;
	double suma = 0.0 , sumb = 0.0;
	int i = 0;
	for(i = 0; i < len; i++) suma += a[i];
	for(i = 0; i < len; i++) sumb += b[i];

	// l2 norm of difference of normalized vector a,b 
	for(i = 0; i < len; i++)
		res += (a[i]/suma - b[i]/sumb) * (a[i]/suma - b[i]/sumb);
	return res;
}

void PrintTopicProportions( FILE* fout , int doc_id , double* gamma , int len )
{
	fprintf( fout , "(docs %d) topic proportion : " , doc_id  );
	int k = 0;
	for(k = 0; k < len; k++)
		fprintf( fout , "T%d:%0.5lf " , k , gamma[k] );
	fprintf( fout , "\n");
}

double CalcTopicL2Distance( double** a , double** b , int num_topics , int num_terms )
{
	double res_dis = 0.0;
	int k = 0;
	for(k = 0; k < num_topics; k++)
		res_dis += NormalizedL2Distance( a[k] , b[k] , num_terms );
	return res_dis;
}

double ProbOnRankInTopic( double* a , int len , int num_topwords )
{
	double st = 0 , en = 1e9 , mid;
	int i , num_higher = 0;
	while( en - st > 1e-5 )
	{
		mid = (st + en) / 2;
		num_higher = 0;
		for(i = 0; i < len; i++) 
			if( a[i] > mid )	
				num_higher += 1;

		if( num_higher > num_topwords )
			st = mid;
		else
			en = mid;
	}
	return (en);
}

double NumOverlapTopWordTopics(double* a , double* b , int len , int num_topwords )
{
	double thres_a = ProbOnRankInTopic( a , len , num_topwords );
	double thres_b = ProbOnRankInTopic( b , len , num_topwords );

	double res = num_topwords;
	int i = 0;

	// calculate num_topwords - overlap in topwords 
	for(i = 0; i < len; i++)
		if( a[i]>=thres_a && b[i]>=thres_b )
			res -= 1.0;
	return (res);
}

double SumVector( double* x, int len)
{
	double sum = 0.0;
	int i = 0;
	for(i = 0; i < len; i++)
		sum += x[i];
	return (sum);
}

int WordIdDictionary( std::string word_str , char** dict , int len_dict )
{
	int i = 0;
	for(i = 0; i < len_dict; i++)
		if( strcmp( word_str.c_str() , dict[i] ) == 0 )
			return (i);
}

/*
select the topic id which contains word topic_theme at the highest rank.
*/
int FindAttackTopicId( std::string topic_theme , char** dict , int num_terms , int num_topics , double** eta )
{
	int select_topic = 0 , highest_rank = num_terms+1;
	int word_id = WordIdDictionary( topic_theme , dict , num_terms );
	int k = 0;
	for(k = 0; k < num_topics; k++)
	{
		int tmp_rank = WordRankInTopic( num_terms , k , word_id , eta );
		if( tmp_rank < highest_rank )
		{
			highest_rank = tmp_rank;
			select_topic = k;
		}
	}
	return (select_topic);
}

int WordRankInTopic( int num_terms , int topic_id , int word_id , double** eta )
{
	int rank = 0 , v = 0;
	for(v = 0; v < num_terms; v++)
		if( eta[topic_id][v] > eta[topic_id][word_id] )
			rank += 1;
	return (rank);
}

void PrintTopicPermutation( FILE* fout , int* permute , int len )
{
	fprintf(fout , "permute:\n");
	int k = 0;
	for(k = 0; k < len; k++)
		fprintf( fout , "%d " , permute[k] );
	fprintf(fout , "\n");
}

const double kEpsilon = 0.001;
double EpsilonInsensitiveSquare(double x) 
{ 
	if( x > -kEpsilon  && x < kEpsilon )
		return 0.0;
	else
	{	
		//x = fabs(x) - kEpsilon;
		return (fabs(x)-kEpsilon)*(fabs(x)-kEpsilon); 
	}
}