#include <vector>
#include "lda-data.h"
#include "cokus.h"
#include "utils.h"

namespace attack_lda {

class GibbsSampling {
 public:
	int num_topics_;
	int num_terms_;
	int num_docs_;
	double alpha_, beta_;
	int n_iters_;
	std::vector<std::vector<int> > nw_;
	std::vector<std::vector<int> > nd_;
	std::vector<int> nw_sum_;
	std::vector<int> nd_sum_;
	std::vector<std::vector<int> > z_;
	std::vector<double> p_;

	std::vector<std::vector<double> > theta_;
	std::vector<std::vector<double> > phi_;
	
	int InitEstimation(Corpus* corpus, int num_topics, int n_iters, double alpha, double beta);
	void Estimate(Corpus* corpus);
	int Sampling(int w, int d, int n);
	void ComputeTheta();
	void ComputePhi(double** phi);
	void RunGibbsSampling(Corpus* corpus, int num_topics, int n_iters, double alpha, double beta, double** phi);
};
} // namespace attack_lda