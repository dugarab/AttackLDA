#include "gibbs-sampling.h"

namespace attack_lda {

int GibbsSampling::InitEstimation(Corpus* corpus, int num_topics, int n_iters, double alpha, double beta) {	
    num_docs_ = corpus->num_docs;
	num_terms_ = corpus->num_terms;
	num_topics_ = num_topics;
	n_iters_ = n_iters;
	alpha_ = alpha;
	beta_ = beta;
	
	// + allocate memory and assign values for variables
	nw_.resize(num_topics_, std::vector<int>());
	phi_.resize(num_topics_, std::vector<double>());
    for(int k = 0; k < num_topics_; k++) {
		nw_[k].resize(num_terms_);
		phi_[k].resize(num_terms_);
        for (int v = 0; v < num_terms_; v++) {
    	    nw_[k][v] = 0;
			phi_[k][v] = 0.0;
        }
    }
	
	nd_.resize(num_docs_, std::vector<int>());
	theta_.resize(num_docs_, std::vector<double>());
    for(int d = 0; d < num_docs_; d++) {
		nd_[d].resize(num_topics_);
		theta_[d].resize(num_topics_);
        for (int k = 0; k < num_topics_; k++) {
    	    nd_[d][k] = 0;
			theta_[d][k] = 0.0;
        }
    }

	nw_sum_.resize(num_topics_);
    for (int k = 0; k < num_topics_; k++) {
		nw_sum_[k] = 0;
    }
    
	nd_sum_.resize(num_docs_);
    for (int d = 0; d < num_docs_; d++) {
		nd_sum_[d] = 0;
    }

	z_.resize(num_docs_);
    for(int d = 0; d < num_docs_; d++) {
		int len_doc = corpus->docs[d].total + 1;
		printf("%d %d\n", d, len_doc);
		z_[d].clear();
		// initialize for z_
		for (int n = 0; n < (*corpus->docs[d].counts).size(); n++) {
			for(int j = 0; j < (*corpus->docs[d].counts)[n]; j++) {
				int topic = (int)((double)myrand() * num_topics_);
				z_[d].push_back(topic);
				// number of instances of word i assigned to topic j
				nw_[topic][(*corpus->docs[d].words)[n]] += 1;
    			// number of words in document i assigned to topic j
    			nd_[d][topic] += 1;
    			// total number of words assigned to topic j
    			nw_sum_[topic] += 1;
			} 
		}
        // total number of words in document i
        nd_sum_[d] = len_doc;      
	}
		printf("I am here 8	\n");
	p_.resize(num_topics_);
    return 0;
}

void GibbsSampling::Estimate(Corpus* corpus) {
    printf("Sampling %d iterations!\n", n_iters_);

    for(int liter = 1; liter <= n_iters_; liter++) {
		printf("Iteration %d ...\n", liter);
	
		// for all z_i
		for(int d = 0; d < num_docs_; d++) {
			int idx = 0;
			for(int n = 0; n < corpus->docs[d].length; n++) {
				for(int j = 0; j < (*corpus->docs[d].counts)[n]; j++) {
		
			// (z_i = z_[d][n])
			// sample from p_(z_i|z_-i, w)
					z_[d][idx] = Sampling((*corpus->docs[d].words)[n], d, idx);
					idx += 1;
				}
			}
		}
	}
    
    printf("Gibbs sampling completed!\n");
}

int GibbsSampling::Sampling(int w, int d, int n) {
    // remove z_i from the count variables
    int topic = z_[d][n];
    nw_[topic][w] -= 1;
    nd_[d][topic] -= 1;
    nw_sum_[topic] -= 1;
    nd_sum_[d] -= 1;

    double Vbeta = num_terms_ * beta_;
    double Kalpha = num_topics_ * alpha_;    
    // do multinomial sampling via cumulative method
    for(int k = 0; k < num_topics_; k++) {
		p_[k] = (nw_[k][w] + beta_) / (nw_sum_[k] + Vbeta) *
			    (nd_[d][k] + alpha_) / (nd_sum_[d] + Kalpha);
    }
    // cumulate multinomial parameters
    for(int k = 1; k < num_topics_; k++) {
		p_[k] += p_[k - 1];
    }
    // scaled sample because of unnormalized p_[]
    double u = ((double)myrand()) * p_[num_topics_ - 1];
	// sweep the topics to find the sample
    for (topic = 0; topic < num_topics_; topic++) {
		if (p_[topic] > u) {
			break;
		}
    }
	// add newly estimated z_i to count variables
    nw_[topic][w] += 1;
    nd_[d][topic] += 1;
    nw_sum_[topic] += 1;
    nd_sum_[d] += 1;    
    
    return topic;
}

void GibbsSampling::ComputeTheta() {
    for (int d = 0; d < num_docs_; d++) {
		for (int k = 0; k < num_topics_; k++) {
		    theta_[d][k] = (nd_[d][k] + alpha_) / (nd_sum_[d] + num_topics_ * alpha_);
		}
    }
}

void GibbsSampling::ComputePhi(double** phi) {
    for (int k = 0; k < num_topics_; k++) {
		for (int w = 0; w < num_terms_; w++) {
			phi[k][w] = nw_[k][w] + beta_;
		}
    }
}

void GibbsSampling::RunGibbsSampling(Corpus* corpus, int num_topics, int n_iters, double alpha, double beta, double** phi) {
	printf("I am here 2\n");
	InitEstimation(corpus, num_topics, n_iters, alpha, beta);
	printf("I am here 3\n");
	Estimate(corpus);
	ComputeTheta();
	ComputePhi(phi);
	for (int k = 0; k < num_topics_; k++) {
		for (int w = 0; w < 100; w++) {
			printf("%0.3lf", phi[k][w]);
		}
		printf("\n");
    }
}

} // namespace attack_lda