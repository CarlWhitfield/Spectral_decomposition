#include<network_3D.h>
#include<list_template.h>
#include<Eigen\Dense>
#include<Eigen\Sparse>
#include<SymEigsSolver.h>
#include <SymEigsShiftSolver.h>
#include<GenEigsSolver.h>
#include<MatOp\SparseGenMatProd.h>
#include<MatOp\SparseSymMatProd.h>
#include<MatOp\SparseSymShiftSolve.h>
#include<algorithm>
#include<numeric>
#include<chrono>
#include"r_matrices.h"

typedef network::Network<network::Node,network::Edge<network::Node>> Tree;
//choices for network TYPE
#define RESISTANCE_NETWORK 0
#define DIFFUSION_NETWORK 1
#define VISCOSITY 1.93E-07  //cmH20s
#define DIFFUSIVITY 1

template<class SpectraSolver, class SpectraOp> int do_solve(SpectraOp * op, const int & Neigs, const int & eig_param,
							     const int & mat_dim, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	int ncv = eig_param;
	bool repeat = true;
	while(repeat && eig_param <= mat_dim)
	{
		SpectraSolver solver(op, Neigs, ncv);
		solver.init();
		solver.compute();
		if(solver.info() == Spectra::SUCCESSFUL)
		{
			repeat = false;
			get_evals_and_evecs<SpectraSolver>(&solver, mat_dim, evalues, evectors);
			return 0;
		}
		else
		{
			ncv = std::min(2*ncv, mat_dim);
			std::cout << "Computation failed, repeating...\n";
		}
	}
	return 1;
}

template<class SpectraSolver, class SpectraOp> int do_shift_solve(SpectraOp * op, const int & Neigs, const int & eig_param,
					const int & mat_dim, const double &s, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	int ncv = eig_param;
	bool repeat = true;
	while(repeat && eig_param <= mat_dim)
	{
		SpectraSolver solver(op, Neigs, ncv, s);
		solver.init();
		solver.compute();
		if(solver.info() == Spectra::SUCCESSFUL)
		{
			repeat = false;
			get_evals_and_evecs<SpectraSolver>(&solver, mat_dim, evalues, evectors);
			return 0;
		}
		else
		{
			ncv = std::min(2*ncv, mat_dim);
			std::cout << "Computation failed, repeating...\n";
		}
	}
	return 1;
}

template<class SpectraSolver> void get_evals_and_evecs(SpectraSolver * solver, const int & mat_dim, Eigen::VectorXd & evalues, 
													                                   std::vector<Eigen::VectorXd> & evectors)
{
	evalues = solver->eigenvalues();
	evectors.resize(evalues.size());
	for(int k = 0; k < evalues.size(); k++)
	{
		evectors[k] = solver->eigenvectors().block(0, k, mat_dim, 1);
	}
}

template<class SpectraMatrixClass, class SpectraShiftSolveClass> class SymmSolver    //class when SpectraMatrix class is a r matric
{
private:
	SpectraMatrixClass matrix, fs_mat;
	SpectraShiftSolveClass *fs_solve;
	void sort_and_return(Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
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
public:
	SymmSolver(const SpectraMatrixClass & mat)
	{
		this->matrix = mat;
		Eigen::VectorXd conductances = Eigen::VectorXd::Ones(this->matrix.return_rvec().size()).array() 
			                         / this->matrix.return_rvec().array();
		this->fs_mat =  SpectraMatrixClass(this->matrix.return_tree_pointer(), -conductances);
		this->fs_solve = new SpectraShiftSolveClass(this->matrix.return_tree_pointer(), -conductances);
	}

	int calc_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
	{
		int eig_param = std::min(2*Neigs, this->matrix.rows());
		int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SpectraMatrixClass>, SpectraMatrixClass>
			          (&this->matrix, Neigs, eig_param, int(this->matrix.rows()), evalues, evectors);
		this->sort_and_return(evalues, evectors);
		return flag;
	}

	int flip_and_shift_and_calc_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, 
		                                            std::vector<Eigen::VectorXd> & evectors, const double & s)
	{
		int eig_param = std::min(5*Neigs, this->matrix.rows());
		this->fs_mat.set_shift(s);
		int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SpectraMatrixClass>, SpectraMatrixClass>
			               (&this->fs_mat, Neigs, eig_param, int(this->matrix.rows()), evalues, evectors);
		//return evalues of matrix
		evalues += s*Eigen::VectorXd::Ones(Neigs);
		evalues = -evalues;
		this->sort_and_return(evalues, evectors);
		return flag;
	}

	int flip_shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, 
		                                   std::vector<Eigen::VectorXd> & evectors, const double & s)
	{
		int eig_param = std::min(5*Neigs, this->matrix.rows());
		int flag = do_shift_solve<Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_ALGE, SpectraShiftSolveClass>, 
			SpectraShiftSolveClass>(this->fs_solve, Neigs, eig_param, int(this->matrix.rows()), s, evalues, evectors);
		//flip back
		evalues = -evalues;
		this->sort_and_return(evalues, evectors);
		return flag;
	}
};

//specialisation for Sparse case
template<> SymmSolver<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>>::SymmSolver
	(const Eigen::SparseMatrix<double> & mat)
{
	this->matrix = mat;
	this->fs_mat = -mat;
	this->fs_solve = new Spectra::SparseSymShiftSolve<double>(this->fs_mat);
}

template<> int SymmSolver<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>>::calc_largest_evalues
	                          (const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	int eig_param = std::min(2*Neigs, int(this->matrix.rows()));
	Spectra::SparseSymMatProd<double> prod_op(this->matrix);
	int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseSymMatProd<double>>, 
		          Spectra::SparseSymMatProd<double>>(&prod_op, Neigs, eig_param, int(this->matrix.rows()), evalues, evectors);
	this->sort_and_return(evalues, evectors);
	return flag;
}

template<> int SymmSolver<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>>::flip_and_shift_and_calc_largest_evalues
	                  (const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors, const double & s)
{
	int eig_param = std::min(5*Neigs, int(this->matrix.rows()));
	Eigen::SparseMatrix<double> id(this->matrix.rows(), this->matrix.cols());
	id.setIdentity();
	Eigen::SparseMatrix<double> shifted_mat = this->fs_mat - s*id;
	Spectra::SparseSymMatProd<double> prod_op(shifted_mat);
	int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseSymMatProd<double>>, 
		         Spectra::SparseSymMatProd<double>>(&prod_op, Neigs, eig_param, int(this->matrix.rows()), evalues, evectors);
	//return evalues of matrix
	evalues += s*Eigen::VectorXd::Ones(Neigs);
	evalues = -evalues;
	this->sort_and_return(evalues, evectors);
	return flag;
}

template<> int SymmSolver<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>>::flip_shift_solve_largest_evalues
	                (const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors, const double & s)
{
	int eig_param = std::min(5*Neigs, int(this->matrix.rows()));
	bool repeat = true;
	this->fs_solve->set_shift(s);
	int flag = do_shift_solve<Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseSymShiftSolve<double>>, 
		           Spectra::SparseSymShiftSolve<double>>(this->fs_solve, Neigs, eig_param, int(this->matrix.rows()), s, evalues, evectors);
	//flip back
	evalues = -evalues;
	this->sort_and_return(evalues, evectors);
	return flag;
}

template<class SpectraMatrixClass, class SpectraShiftSolveClass> void compute_partial_spectrum(SpectraMatrixClass & op, const size_t & mat_dim, const int & Nsmallest,
	                              const int & Nlargest, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	Eigen::VectorXd largest_evals;
	std::vector<Eigen::VectorXd> largest_evecs;
	int Neigs = Nlargest;
	evalues.resize(Nsmallest + Nlargest);
	if(Nlargest == 0)
	{
		Neigs = 1;
	}
	SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass> solver(op);
	solver.calc_largest_evalues(Neigs, largest_evals, largest_evecs);
	Eigen::VectorXd smallest_evals;
	std::vector<Eigen::VectorXd> smallest_evecs;
	if(Nsmallest > 0)
	{
		//get largest eval
		double lambda0 = largest_evals[0];
		solver.flip_and_shift_and_calc_largest_evalues(Nsmallest, smallest_evals, smallest_evecs, -lambda0); 
	}

	//concatenate
	if(Nsmallest > 0 && Nlargest > 0)
	{
		evalues.head(Nlargest) = largest_evals;
		evalues.tail(Nsmallest) = smallest_evals;
		evectors = largest_evecs;
		evectors.insert(evectors.end(), smallest_evecs.begin(), smallest_evecs.end());
	}
	else
	{
		if(Nsmallest > 0)
		{
			evalues = smallest_evals;
			evectors = smallest_evecs;
		}
		if(Nlargest > 0)
		{
			evalues = largest_evals;
			evectors = largest_evecs;
		}
	}
}

template<class SpectraMatrixClass, class SpectraShiftSolveClass> void compute_full_spectrum(SpectraMatrixClass & op, const size_t & mat_dim, 
													Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	auto start = std::chrono::system_clock::now();
	if(mat_dim <= 1000)
	{
		compute_partial_spectrum<SpectraMatrixClass, SpectraShiftSolveClass>(op, mat_dim, size_t(0.51*mat_dim), mat_dim - size_t(0.51*mat_dim), evalues, evectors);
	}
	else
	{
		//get largest first, then get next largst by shifting
		int step = 10;
		SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass> solver(op);
		evalues = Eigen::VectorXd::Ones(step+1);
		while(std::fabs(evalues[evalues.size()-1] - evalues[evalues.size()-2]) < 1E-06*std::fabs(evalues[evalues.size()-1]))
		{
			solver.calc_largest_evalues(step+1,evalues,evectors);
			if(std::fabs(evalues[evalues.size()-1] - evalues[evalues.size()-2]) < 1E-06*std::fabs(evalues[evalues.size()-1]))
			{
				step += 1;
			}
		}
		
		while(evalues.size() < int(mat_dim))
		{
			std::cout << evalues.size() << "\n";
			//start again between last two
		    double lambda0 = 0.5*(evalues[evalues.size()-1] + evalues[evalues.size()-2]);
			//remove last
			evalues.conservativeResize(evalues.size()-1);
			evectors.pop_back();
			Eigen::VectorXd h_evals = Eigen::VectorXd::Ones(step+1);
			std::vector<Eigen::VectorXd> h_evecs;
			if(step > int(mat_dim - evalues.size()) - 1) 
			{
				step = int(mat_dim - evalues.size()) - 1;
			}
			while(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1])
				  && step < int(mat_dim - evalues.size()))
			{
				if(solver.flip_shift_solve_largest_evalues(step+1, h_evals, h_evecs, -lambda0))
				{
					std::cout << "Error computing spectrum, aborting...\n";
					abort_on_failure();
				}
				if(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1]))
				{
					step += 1;
				}
			}
			Eigen::VectorXd evnew(evalues.size() + h_evals.size());
			evnew << evalues, h_evals;
			evalues = evnew;
			evectors.insert(evectors.end(), h_evecs.begin(), h_evecs.end());
		}
	}
	auto end = std::chrono::system_clock::now();
	std::cout << "Spectral calc took: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "s.\n";
}

template<int TYPE> double edge_weight_calc(const double & rad, const double & length, const double & scale)
{
	return (scale*rad*rad*rad*rad/length);
}

template<> double edge_weight_calc<DIFFUSION_NETWORK>(const double & rad, const double & length, const double & scale)
{
	return (scale*rad*rad/length);
}

template<int TYPE> class SpectralNetwork: public Tree
{
private:
	double Ceff;

	Eigen::VectorXd truncated_laplacian_evalues, maury_evalues, adjacency_evalues, truncated_laplacian_reffs, maury_ceffs;
	std::vector<Eigen::VectorXd> truncated_laplacian_evectors, maury_evectors, adjacency_evectors;

	Eigen::VectorXd edge_weight, flux, pressure, Nt, laplacian_flux, laplacian_pressure, maury_flux, maury_pressure;

	Eigen::SparseMatrix<double> Incidence, Degree, Adjacency, Laplacian, Truncated_Laplacian;
	Rmat<RMAT1> Maury_matrix;
	void build_asymmetric_network(const size_t & Ngens, const double & rad0, const double & length0,
			                          const double & scale_factor, const double & asymm_factor);
	void fill_matrices();
	void solve_effective_conductance_problem();
	void fill_fluxes_from_term_vals(Eigen::VectorXd & flux, const Eigen::VectorXd & term_flux);
	void calc_laplacian_reffs();
	void calc_maury_ceffs();
	void print_dominant_vectors_vtk(const std::string & filename, const size_t & Nvecs, const Eigen::VectorXd & evalues,
									const Eigen::VectorXd & contributions, const std::vector<Eigen::VectorXd> & evectors);
public:
	SpectralNetwork():Tree(){};
	SpectralNetwork(const std::vector<network::Node*> & n, const std::vector<network::Edge<network::Node>*> & e):Tree(n,e)
	{
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}
	SpectralNetwork(const std::map<int, network::Node*> & n, const std::map<int, network::Edge<network::Node>*> & e):Tree(n,e)
	{
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}
	SpectralNetwork(const std::string & node_fname, const std::string & edge_fname, const std::string & term_fname):Tree(node_fname, edge_fname, term_fname, 1.0)
	{
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}
	SpectralNetwork(const size_t & Ngens, const double & rad0, const double & length0,
			          const double & scale_factor, const double & asymm_factor)
	{
		this->build_asymmetric_network(Ngens, rad0, length0, scale_factor, asymm_factor);
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}

	inline const Eigen::VectorXd& get_edge_weights(){ return this->edge_weight; }

	void compute_truncated_laplacian_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		size_t mat_dim = this->Truncated_Laplacian.rows();
		compute_partial_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>>
			                                 (this->Truncated_Laplacian, mat_dim, Nsmallest, Nlargest,
								                     this->truncated_laplacian_evalues, this->truncated_laplacian_evectors);
		this->calc_laplacian_reffs();
	}
	void compute_maury_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		size_t mat_dim = this->Maury_matrix.rows();
		compute_partial_spectrum<Rmat<RMAT1>,Rmatinv<RMAT1>>(this->Maury_matrix, mat_dim, Nsmallest, Nlargest,
			                                  this->maury_evalues, this->maury_evectors);
		//do ceff calc
		this->calc_maury_ceffs();
	}

	void compute_full_truncated_laplacian_spectrum()
	{
		size_t mat_dim = this->Truncated_Laplacian.rows();
		compute_full_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>>
			         (this->Truncated_Laplacian, mat_dim,  this->truncated_laplacian_evalues, this->truncated_laplacian_evectors);
		this->calc_laplacian_reffs();
	}
	void compute_full_maury_spectrum()
	{
		size_t mat_dim = this->Maury_matrix.rows();
		compute_full_spectrum<Rmat<RMAT1>,Rmatinv<RMAT1>>(this->Maury_matrix, mat_dim, this->maury_evalues, this->maury_evectors);
		//do ceff calc
		this->calc_maury_ceffs();
	}
	void compute_adjacency_specturm(const size_t & Nsmallest, const size_t & Nlargest)
	{
		compute_partial_spectrum<Eigen::SparseMatrix<double>>(this->Adjacency, this->count_nodes(),
		            Nsmallest, Nlargest, this->adjacency_evalues, this->adjacencey_evectors);
	}
	void compute_full_adjacency_specturm()
	{
		compute_full_spectrum<Eigen::SparseMatrix<double>>(this->Adjacency, this->count_nodes(),
		                            this->adjacency_evalues, this->adjacencey_evectors);
	}

	void print_dominant_laplacian_evectors_vtk(const std::string & filename, const size_t & Nvecs)
	{
		//copy truncated eigenvectors to full network
		std::vector<Eigen::VectorXd> node_evecs;
		node_evecs.resize(this->truncated_laplacian_evectors.size());
		for(size_t n = 0; n < this->truncated_laplacian_evectors.size(); n++)
		{
			node_evecs[n] = Eigen::VectorXd::Zero(this->count_nodes());
			//all internal nodes remain the same
			node_evecs[n].head(this->truncated_laplacian_evectors[n].size()) = this->truncated_laplacian_evectors[n];
			//all term nodes have same value
			node_evecs[n].tail(this->count_nodes()-this->truncated_laplacian_evectors[n].size())
				= this->truncated_laplacian_evectors[n][this->truncated_laplacian_evectors[n].size()-1]
			    * Eigen::VectorXd::Ones(this->count_nodes()-this->truncated_laplacian_evectors[n].size());
		}
		//print
		this->print_dominant_vectors_vtk(filename, Nvecs, this->truncated_laplacian_evalues, this->truncated_laplacian_reffs, node_evecs);
	}
	void print_dominant_maury_evectors_vtk(const std::string & filename, const size_t & Nvecs)
	{
		//convert maury eigenvectors into edge values for visualisation purposes
		std::vector<Eigen::VectorXd> edge_evecs;
		edge_evecs.resize(this->maury_evectors.size());
		for(size_t n = 0; n < this->maury_evectors.size(); n++)
		{
			this->fill_fluxes_from_term_vals(edge_evecs[n], this->maury_evectors[n]);
			edge_evecs[n] = edge_evecs[n].array() / this->Nt.array();
		}
		//print
		this->print_dominant_vectors_vtk(filename, Nvecs, this->maury_evalues, this->maury_ceffs, edge_evecs);
	}
	void print_unit_pressure_drop_solutions_vtk(const std::string & filename);
	void print_laplacian_modes_csv(const std::string & filename);
	void print_maury_modes_csv(const std::string & filename);
};

template<int TYPE> void SpectralNetwork<TYPE>::build_asymmetric_network(const size_t & Ngens, const double & rad0, const double & length0,
			                          const double & scale_factor, const double & asymm_factor)
{
	std::vector<network::Node*> nodes;
	std::vector<network::Edge<network::Node>*> edges;
	std::vector<double> theta_vals;
	nodes.resize(size_t(pow(2,Ngens)));
	theta_vals.resize(nodes.size());
	edges.resize(nodes.size()-1);

	nodes[0] = new network::Node(0,0,length0);
	size_t node_count = 1;
	for(size_t j = 0; j < Ngens; j++)
	{
		for(size_t k = 0; k < pow(2,j); k++)
		{
			network::Position dir;
			size_t node_in;
			double radius, length;
			if(j==0)
			{
				node_in = 0;
				radius = rad0;
				length = length0;
				dir = network::Position(0,0,-1);
			}
			else
			{
				node_in = size_t(pow(2,j-1) + int(k/2));
				radius = pow(0.5*(1.0 + pow(-1.0,k)*asymm_factor),1.0/scale_factor) * edges[node_in - 1]->get_geom()->get_inner_radius();
				length = pow(0.5*(1.0 + pow(-1.0,k)*asymm_factor),1.0/scale_factor) * edges[node_in - 1]->get_geom()->get_length();

				network::Position pos0 = nodes[node_in]->get_pos();
				double r0 = pos0.magnitude();
				double theta0 = 0;
				if(pos0.x[0] != 0)
				{
					theta0 = atan2(pos0.x[1],pos0.x[0]);
				}
				double dtheta = M_PI*(-1.0 + 2*(k+0.5)/pow(2,j)) - theta0;
                double r1;
				if (length < r0*sqrt(1.0 - cos(dtheta)*cos(dtheta)))
                {
					r1 = r0*cos(dtheta);
				}
                else
				{
                    r1 = r0*cos(dtheta) + sqrt(length*length - r0*r0*(1-cos(dtheta)*cos(dtheta)));
				}

				dir = network::Position(r1*cos(theta0 + dtheta) - r0*cos(theta0), r1*sin(theta0 + dtheta) - r0*sin(theta0), 0);
				dir.normalise();
			}
			nodes[node_count] = new network::Node(nodes[node_in]->get_pos() + dir*length);
			edges[node_count-1] = new network::Edge<network::Node>(nodes[node_in], nodes[node_count], 1, radius);
			node_count++;
		}
	}
	this->Tree::Network<network::Node,network::Edge<network::Node>>(nodes,edges);
}

template<int TYPE> void SpectralNetwork<TYPE>::fill_matrices()
{
	std::cout << "Filling matrices...\n";
	auto start = std::chrono::system_clock::now();
	//fill vector of edge weights
	this->edge_weight = Eigen::VectorXd::Zero(this->count_edges());
	std::cout << this->get_edge(0)->get_geom()->get_inner_radius() << ' ' << this->get_edge(0)->get_geom()->get_length() << '\n';
	for(size_t j = 0; j < this->count_edges(); j++)
	{
		double scale =  this->get_edge(j)->branch_count();
		if(TYPE == RESISTANCE_NETWORK)
		{
			scale *= M_PI*1E-06 / (8*VISCOSITY);   //conductance in L / cmH20 s, airway geom in mm
		}
		if(TYPE == DIFFUSION_NETWORK)
		{
			scale *= M_PI * DIFFUSIVITY;
		}
		this->edge_weight[j] = edge_weight_calc<TYPE>(this->get_edge(j)->get_geom()->get_inner_radius(),
				                                      this->get_edge(j)->get_geom()->get_length(), scale);
	}
	//initialise various matrices
	this->Degree = Eigen::SparseMatrix<double>(this->count_nodes(), this->count_nodes());
	this->Adjacency = Eigen::SparseMatrix<double>(this->count_nodes(), this->count_nodes());
	this->Incidence = Eigen::SparseMatrix<double>(this->count_nodes(), this->count_edges());
	std::vector<Eigen::Triplet<double>> degree_fill, adjacency_fill, incidence_fill;
	degree_fill.reserve(this->count_nodes());
	adjacency_fill.reserve(3*this->count_nodes());
	incidence_fill.reserve(3*this->count_nodes());

	//fill matrices - loop over all nodes
	for(size_t k = 0; k < this->count_nodes(); k++)
	{
		double deg = 0;
		//loop over edges in
		for(size_t ji = 0; ji < this->count_edges_in(k); ji++)
		{
			size_t j = this->get_edge_in_index(k,ji);
			incidence_fill.push_back(Eigen::Triplet<double>(int(k), int(j), -1.0)); //incidence matrix

			//adjacency contribution
			size_t ki = this->get_node_in_index(j);
			adjacency_fill.push_back(Eigen::Triplet<double>(int(k), int(ki), this->edge_weight[j]));

			deg += this->edge_weight[j];    //count degree
		}
		//loop over edges out
		for(size_t jo = 0; jo < this->count_edges_out(k); jo++)
		{
			size_t j = this->get_edge_out_index(k,jo);
			incidence_fill.push_back(Eigen::Triplet<double>(int(k), int(j), 1.0));  //incidence matrix

			//adjacency contribution
			size_t ko = this->get_node_out_index(j);
			adjacency_fill.push_back(Eigen::Triplet<double>(int(k), int(ko), this->edge_weight[j]));

			deg += this->edge_weight[j];    //count degree
		}
		degree_fill.push_back(Eigen::Triplet<double>(int(k), int(k), deg));
	}
	//build matrices
	this->Adjacency.setFromTriplets(adjacency_fill.begin(),adjacency_fill.end());
	this->Degree.setFromTriplets(degree_fill.begin(),degree_fill.end());
	this->Incidence.setFromTriplets(incidence_fill.begin(),incidence_fill.end());
	this->Laplacian = this->Degree - this->Adjacency;

	//create tuncated laplacian where all terminal nodes are identified into a single node
	//              |            |   |
	//              |      A     | B |
	//    TL   =    |            |   |   where A = Laplacian.topLeftCorner(TLsize,Tlsize), B=Laplacian.topRightCorner(TLsize,Nterm)*Ones(Nterm)  (Tlsize x 1)
	//              |____________|___|         C = Ones(Nterm).transpose()*Laplacian.bottomLeftCorner(Nterm,TLsize)    (1 x TLsize)
	//              |      C     | D |         D = Ones(Nterm).transpose()*Laplacian.bottomRightCorner(Nterm,Nterm)*Ones(Nterm)  (1 x 1)

	size_t TLsize = this->count_nodes() - this->count_term_nodes();
	this->Truncated_Laplacian = Eigen::SparseMatrix<double>(TLsize+1, TLsize+1);
	std::vector<Eigen::Triplet<double>> TL_fill;
	TL_fill.reserve(4*this->count_nodes());
	//fill A and C
	for (int k=0; k<TLsize; ++k)
	{
		double row_sum = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(this->Laplacian,k); it; ++it)
		{
			if(size_t(it.row()) < TLsize)  //fill A
			{
				TL_fill.push_back(Eigen::Triplet<double>(int(it.row()), int(it.col()), it.value()));
			}
			else  //sum over rows
			{
				row_sum += it.value();
			}
		}
		//fill C
		TL_fill.push_back(Eigen::Triplet<double>(int(TLsize), k, row_sum));
	}
	Eigen::VectorXd row_sum = Eigen::VectorXd::Zero(TLsize);
	double diag_sum = 0;
	//top right is sum of rows, bottom right is sum of rest
	for (size_t k=TLsize; k<this->count_nodes(); k++) //loop over columns
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(this->Laplacian,int(k)); it; ++it)  //loop over rows
		{
			if(it.row() < int(TLsize))
			{
				row_sum[it.row()] += it.value();
			}
			else
			{
				diag_sum += it.value();
			}
		}
	}
	//fill B
	for(size_t k=0; k<TLsize; k++) //loop over rows
	{
		TL_fill.push_back(Eigen::Triplet<double>(int(k), int(TLsize), row_sum[k]));
	}
	//fill D
	TL_fill.push_back(Eigen::Triplet<double>(int(TLsize), int(TLsize), diag_sum));
	this->Truncated_Laplacian.setFromTriplets(TL_fill.begin(),TL_fill.end());
	this->Degree.makeCompressed();
	this->Adjacency.makeCompressed();
	this->Laplacian.makeCompressed();
	this->Truncated_Laplacian.makeCompressed();

	//std::cout << this->Truncated_Laplacian.toDense() << "\n\n";

	//this->Laplacian = this->Incidence * Eigen::DiagonalMatrix(this->edge_weight) *  this->Incidence.transpose();
	this->Maury_matrix = Rmat<RMAT1>(this, this->edge_weight);
	this->fill_fluxes_from_term_vals(this->Nt, Eigen::VectorXd::Ones(this->count_term_nodes()));
	auto end = std::chrono::system_clock::now();
	std::cout << "Matrix fill took " << (std::chrono::duration<double>(end-start)).count() << '\n';
}

template<int TYPE> void SpectralNetwork<TYPE>::solve_effective_conductance_problem()
{
	//define matrix
	std::cout << "Solving pressure drop case...\n";
	auto start = std::chrono::system_clock::now();
	Eigen::SparseMatrix<double> A((this->count_nodes()+this->count_term_nodes()+1), (this->count_nodes()+this->count_term_nodes()+1));
	//cout << this->count_nodes() << '\n';
	//cout << this->count_nodes()+this->count_term_nodes()+1 << '\n';
	std::vector<Eigen::Triplet<double>> A_fill;
	A_fill.reserve(4*this->count_nodes() + 3*this->count_term_nodes());
	//first Nnodes*Nnodes block is laplacian
	for(size_t i = 0; i < this->count_nodes(); i++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(this->Laplacian,i); it; ++it)
		{
			A_fill.push_back(Eigen::Triplet<double>(int(it.row()), int(it.col()), it.value()));
		}
	}
	//enforce pressure boundary conditions
	for(size_t i = this->get_first_term_index(); i <this->count_nodes(); i++)
	{
		A_fill.push_back(Eigen::Triplet<double>(int(this->count_nodes() + i - this->get_first_term_index()), int(i), 1.0));
		A_fill.push_back(Eigen::Triplet<double>(int(i), int(this->count_nodes() + i - this->get_first_term_index()), 1.0));
	}
	int Msize = ((int) (this->count_nodes() + this->count_term_nodes()));
	A_fill.push_back(Eigen::Triplet<double>(0, Msize, -1.0));
	A_fill.push_back(Eigen::Triplet<double>(Msize, 0, 1.0));
	A.setFromTriplets(A_fill.begin(),A_fill.end());
	A.makeCompressed();

	//std::cout << A.toDense() << "\n\n";
	Eigen::VectorXd Bvec = Eigen::VectorXd::Zero(this->count_nodes() + this->count_term_nodes() + 1);
	Bvec[this->count_nodes() + this->count_term_nodes()] = 1.0;

	//cout << "Populating matrix.\n";
	Eigen::SparseLU<Eigen::SparseMatrix<double>> LUsolve;
	LUsolve.compute(A);
	//cout << LUsolve.lastErrorMessage() << '\n';
	//cout << "Solving.\n";
	Eigen::VectorXd Xsol = LUsolve.solve(Bvec);
	//cout << "Solved.\n";

	this->pressure = Xsol.head(this->count_nodes());
	this->fill_fluxes_from_term_vals(this->flux, Xsol.segment(this->count_nodes(), this->count_term_nodes()));
	this->Ceff = Xsol[this->count_nodes() + this->count_term_nodes()];
	auto end = std::chrono::system_clock::now();
	std::cout << "Pressure drop solve took " << (std::chrono::duration<double>(end-start)).count() << '\n';
}

template<int TYPE> void SpectralNetwork<TYPE>::fill_fluxes_from_term_vals(Eigen::VectorXd & flux, const Eigen::VectorXd & term_flux)
{
	flux = Eigen::VectorXd::Zero(this->count_edges());
	//check input is valid
	if(term_flux.size() != this->count_term_nodes())
	{
		std::cout << "Error, flux vector is not correct size, cannot compute fluxes.\n";
	}
	else
	{
		for(size_t k = this->get_first_term_index(); k < this->count_nodes(); k++)
		{
			size_t kt = k - this->get_first_term_index();
			size_t j = this->get_edge_in_index(k,0);
			flux[j] = term_flux[kt];
		}

		//fill in rest of tree by looping over horsfield orders
		for(size_t ho = 1; ho < this->count_horsfield_orders(); ho++)  //loop over horsfield orders
		{
			for(size_t ji = 0; ji < this->count_edges_in_horsfield_order(ho); ji++)    //loop over edges in each order
			{
				size_t j = this->get_edge_index_from_horsfield_order(ho,ji);     //index of this edge
				size_t k_out = this->get_node_out_index(j);
				for (size_t jo = 0; jo < this->count_edges_out(k_out); jo++)
				{
					size_t eo_index = this->get_edge_out_index(this->get_node_out_index(j),jo);
					flux[j] += flux[eo_index];
				}
			}
		}
	}

}

template<int TYPE> void SpectralNetwork<TYPE>::calc_laplacian_reffs()
{
	this->truncated_laplacian_reffs = Eigen::VectorXd::Zero(this->truncated_laplacian_evalues.size());
	size_t TLsize = this->count_nodes() - this->count_term_nodes();
	for(size_t n = 0; n < size_t(this->truncated_laplacian_evalues.size()); n++)
	{
		this->truncated_laplacian_reffs[n] = (this->truncated_laplacian_evectors[n][0] - this->truncated_laplacian_evectors[n][TLsize])*
			                                 (this->truncated_laplacian_evectors[n][0] - this->truncated_laplacian_evectors[n][TLsize])/
								              this->truncated_laplacian_evalues[n];
	}
	double Q = 1.0 / this->truncated_laplacian_reffs.sum();
	this->laplacian_pressure = Eigen::VectorXd::Zero(this->count_nodes());
	for(size_t n = 0; n < size_t(this->truncated_laplacian_evalues.size()); n++)
	{
		Eigen::VectorXd evect_full = Eigen::VectorXd::Zero(this->count_nodes());
		evect_full.head(TLsize+1) = this->truncated_laplacian_evectors[n];
		evect_full.tail(this->count_nodes()-TLsize) =  this->truncated_laplacian_evectors[n][TLsize]*Eigen::VectorXd::Ones(this->count_nodes()-TLsize);
		this->laplacian_pressure = this->laplacian_pressure + Q*(this->truncated_laplacian_evectors[n][0]
		                         - this->truncated_laplacian_evectors[n][TLsize])*evect_full/this->truncated_laplacian_evalues[n];

	}
	this->laplacian_pressure = this->laplacian_pressure - this->laplacian_pressure[this->count_nodes()-1]*Eigen::VectorXd::Ones(this->count_nodes());   //offset
	this->laplacian_flux = Eigen::VectorXd::Zero(this->count_edges());
	for(size_t j = 0; j < this->count_edges(); j++)
	{
		size_t k_in = this->get_node_in_index(j);
		size_t k_out = this->get_node_out_index(j);
		this->laplacian_flux[j] = (this->laplacian_pressure[k_in] - this->laplacian_pressure[k_out])*this->edge_weight[j];
	}
}

template<int TYPE> void SpectralNetwork<TYPE>::calc_maury_ceffs()
{
	this->maury_ceffs = Eigen::VectorXd::Zero(this->maury_evalues.size());
	Eigen::VectorXd term_fluxes = Eigen::VectorXd::Zero(this->count_term_nodes());
	for(size_t n = 0; n < size_t(this->maury_evalues.size()); n++)
	{
		double vec_sum = this->maury_evectors[n].sum();
		this->maury_ceffs[n] = (vec_sum *vec_sum)/this->maury_evalues[n];
		term_fluxes = term_fluxes + (vec_sum/this->maury_evalues[n])*this->maury_evectors[n];
	}
	this->fill_fluxes_from_term_vals(this->maury_flux, term_fluxes);
	this->maury_pressure = Eigen::VectorXd::Zero(this->count_nodes());
	this->maury_pressure.tail(this->count_term_nodes()) = Eigen::VectorXd::Zero(this->count_term_nodes());
	//fill in rest of tree by looping over horsfield orders
	for(size_t ho = 1; ho < this->count_horsfield_orders(); ho++)  //loop over horsfield orders
	{
		for(size_t ji = 0; ji < this->count_edges_in_horsfield_order(ho); ji++)    //loop over edges in each order
		{
			size_t j = this->get_edge_index_from_horsfield_order(ho,ji);     //index of this edge
			size_t k_out = this->get_node_out_index(j);
			size_t j_out = this->get_edge_out_index(this->get_node_out_index(j),0);
			size_t k_out_out = this->get_node_out_index(j_out);
			this->maury_pressure[k_out] = this->maury_pressure[k_out_out] + this->maury_flux[j_out] / this->edge_weight[j_out];
		}
	}
	size_t j_out = this->get_edge_out_index(0,0);
	this->maury_pressure[0] = this->maury_pressure[this->get_node_out_index(j_out)] + this->maury_flux[j_out] / this->edge_weight[j_out];
}

template<int TYPE> void SpectralNetwork<TYPE>::print_dominant_vectors_vtk(const std::string & filename, const size_t & Nvecs, const Eigen::VectorXd & evalues,
																		  const Eigen::VectorXd & contributions, const std::vector<Eigen::VectorXd> & evectors)
{
	size_t to_print = std::min(Nvecs, size_t(evalues.size()));
	std::vector<size_t> indices(size_t(evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(), [&contributions](size_t i1, size_t i2){ return contributions[i1] > contributions[i2]; });

	std::unordered_map<std::string, std::vector<double>> extra_vals;
	for(size_t n = 0; n < to_print; n++)
	{
		std::stringstream name;
		name << "Dominant_mode_" << n << "_eigenvalue_" << evalues[indices[n]];

		extra_vals[name.str()] = std::vector<double>(evectors[indices[n]].data(), evectors[indices[n]].data() + evectors[indices[n]].size());
	}
	this->print_vtk(filename, 1.0, extra_vals);
}

template<int TYPE> void SpectralNetwork<TYPE>::print_unit_pressure_drop_solutions_vtk(const std::string & filename)
{
	std::unordered_map<std::string, std::vector<double>> extra_vals;
	extra_vals["Flux"] = std::vector<double>(this->flux.data(), this->flux.data() + this->flux.size());
	extra_vals["Pressure"] = std::vector<double>(this->pressure.data(), this->pressure.data() + this->pressure.size());
	if(truncated_laplacian_evalues.size() > 0)
	{
		extra_vals["Laplace_flux_approx"] = std::vector<double>(this->laplacian_flux.data(), this->laplacian_flux.data() + this->laplacian_flux.size());
		extra_vals["Laplace_pressure_approx"] = std::vector<double>(this->laplacian_pressure.data(), this->laplacian_pressure.data() + this->laplacian_pressure.size());
	}
	if(maury_evalues.size() > 0)
	{
		extra_vals["Maury_flux_approx"] = std::vector<double>(this->maury_flux.data(), this->maury_flux.data() + this->maury_flux.size());
		extra_vals["Maury_pressure_approx"] = std::vector<double>(this->maury_pressure.data(), this->maury_pressure.data() + this->maury_pressure.size());
	}

	this->print_vtk(filename, 1.0, extra_vals);
}

template<int TYPE> void SpectralNetwork<TYPE>::print_laplacian_modes_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(2);
	headers[0] = std::string("Eigenvalue");
	headers[1] = std::string("Relative_effective_resistance");
	int nevals = int(this->truncated_laplacian_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals);
	for(size_t n = 0; n < nevals; n++)
	{
		data[n].resize(2);
		data[n][0] = this->truncated_laplacian_evalues[n];
		data[n][1] = this->truncated_laplacian_reffs[n] * this->Ceff;
	}
	write_csv_file<double>(filename, headers, data);
}

template<int TYPE> void SpectralNetwork<TYPE>::print_maury_modes_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(2);
	headers[0] = std::string("Eigenvalue");
	headers[1] = std::string("Relative_effective_conductance");
	int nevals = int(this->maury_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals);
	for(size_t n = 0; n < nevals; n++)
	{
		data[n].resize(2);
		data[n][0] = this->maury_evalues[n];
		data[n][1] = this->maury_ceffs[n] / this->Ceff;
	}

	write_csv_file<double>(filename, headers, data);
}
