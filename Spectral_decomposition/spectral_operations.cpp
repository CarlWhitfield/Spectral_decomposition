#include"spectral_operations.h"

//specialisation for DIFFUSION case
template<> double edge_weight_calc<DIFFUSION_NETWORK>(const double & rad, const double & length, const double & scale)
{
	return (scale*rad*rad/length);
}

void sort_and_return_by_evalues(Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	std::vector<int> indices(evalues.size());
	std::iota(indices.begin(), indices.end(), 0);

	// sort indexes based on comparing values in evalues (descending)
	std::sort(indices.begin(), indices.end(),  [&evalues](int i1, int i2) {return evalues[i1] > evalues[i2];});
	Eigen::VectorXd eval_copy = evalues;
	std::vector<Eigen::VectorXd> evec_copy = evectors;
	for(size_t i = 0; i < indices.size(); i++)
	{
		evalues[i] = eval_copy[indices[i]];
		evectors[i] = evec_copy[indices[i]];
	}
}