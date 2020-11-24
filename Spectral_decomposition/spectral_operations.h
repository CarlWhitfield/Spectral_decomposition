#ifndef SPECTRAL_OPERATIONS_H
#define SPECTRAL_OPERATIONS_H

#include<network_3D.h>
#include<list_template.h>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Spectra/SymEigsSolver.h>
#include<Spectra/SymEigsShiftSolver.h>
#include<Spectra/GenEigsSolver.h>
#include<Spectra/MatOp/SparseGenMatProd.h>
#include<Spectra/MatOp/SparseSymMatProd.h>
#include<Spectra/MatOp/SparseSymShiftSolve.h>
#include<algorithm>
#include<numeric>
#include<chrono>
#include<memory>
#include"r_matrices.h"

typedef network::Network<network::Node,network::Edge<network::Node>> Tree;
//choices for network TYPE
#define RESISTANCE_NETWORK 0
#define DIFFUSION_NETWORK 1
#define VISCOSITY 1.93E-07  //cmH20s
#define DIFFUSIVITY 1

//sort options
#define LARGEST_SORT 1
#define SMALLEST_SORT 2
#define DOMINANT_SORT 3

void sort_and_return_by_evalues(Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors);

//custom spectra solver
template<class SpectraMatrixClass, class SpectraShiftSolveClass, int RTYPE> class SymmSolver
{
protected:
	SpectraMatrixClass matrix, fs_mat;
	SpectraShiftSolveClass *fs_solve;
	SpectraShiftSolveClass *shift_solve;
public:
	SymmSolver(const Eigen::SparseMatrix<double> & mat){};
	~SymmSolver(){};
	int calc_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues,
							 std::vector<Eigen::VectorXd> & evectors){};
	int calc_smallest_evalues(const int & Neigs, Eigen::VectorXd & evalues,
							  std::vector<Eigen::VectorXd> & evectors){};
	int shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues,
									std::vector<Eigen::VectorXd> & evectors, const double & s){};
	int flip_shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues,
										 std::vector<Eigen::VectorXd> & evectors, const double & s){};
};

//retrieve eigenvalues and eigenvectors from solver
template<class SpectraSolver> void get_evals_and_evecs(SpectraSolver * solver, const int & mat_dim, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	evalues = solver->eigenvalues();
	evectors.resize(evalues.size());
	for(int k = 0; k < evalues.size(); k++)
	{
		evectors[k] = solver->eigenvectors().block(0, k, mat_dim, 1);
	}
}

//template functions for performing spectra solve operations
//straight forward eigenvalue solver
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

//shift-solve operation
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

//function to compute parts of the spectrum
template<class SpectraMatrixClass, class SpectraShiftSolveClass, int RTYPE> void compute_partial_spectrum
	                              (SpectraMatrixClass & op, const size_t & mat_dim, const int & Nsmallest,
	                              const int & Nlargest, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	Eigen::VectorXd largest_evals;
	std::vector<Eigen::VectorXd> largest_evecs;
	Eigen::VectorXd smallest_evals;
	std::vector<Eigen::VectorXd> smallest_evecs;
	auto start = std::chrono::system_clock::now();
	if(mat_dim <= 1000)
	{
		int Neigs = Nlargest;
		evalues.resize(Nsmallest + Nlargest);
		if(Nlargest == 0)
		{
			Neigs = 1;
		}	
		SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE> solver(op);
		solver.calc_largest_evalues(Neigs, largest_evals, largest_evecs);
	
		if(Nsmallest > 0)
		{
			//get largest eval
			//double lambda0 = largest_evals[0];
			solver.calc_smallest_evalues(Nsmallest, smallest_evals, smallest_evecs); 
		}
	}
	else
	{
		if(Nlargest > 0)
		{
			//get largest first, then get next largst by shifting
			int step = 10;
			SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE> solver(op);
			largest_evals = Eigen::VectorXd::Ones(step+1);
			while(std::fabs(largest_evals[largest_evals.size()-1] - largest_evals[largest_evals.size()-2]) < 1E-06*std::fabs(largest_evals[largest_evals.size()-1]))
			{
				solver.calc_largest_evalues(step+1,largest_evals,largest_evecs);
				if(std::fabs(largest_evals[largest_evals.size()-1] - largest_evals[largest_evals.size()-2]) < 1E-06*std::fabs(largest_evals[largest_evals.size()-1]))
				{
					step += 1;
					//std::cout << step << std::endl;
				}
			}
		
			while(largest_evals.size() < Nlargest)
			{
				std::cout << largest_evals.size() << "\n";
				//start again between last two
				double lambda0 = 0.5*(largest_evals[largest_evals.size()-1] + largest_evals[largest_evals.size()-2]);
				//remove last
				largest_evals.conservativeResize(largest_evals.size()-1);
				largest_evecs.pop_back();
				step = 10;
				if(step > int(Nlargest - largest_evals.size()) - 1) 
				{
					step = int(Nlargest - largest_evals.size()) - 1;
				}
				Eigen::VectorXd h_evals = Eigen::VectorXd::Ones(step+1);
				std::vector<Eigen::VectorXd> h_evecs;
				while(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1])
					  && step < int(mat_dim - largest_evals.size()))
				{
					if(solver.flip_shift_solve_largest_evalues(step+1, h_evals, h_evecs, -lambda0))
					{
						std::cout << "Error computing spectrum, aborting...\n";
						abort_on_failure();
					}
					if(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1]))
					{
						step += 11;
					}
				}
				Eigen::VectorXd evnew(largest_evals.size() + h_evals.size());
				evnew << largest_evals, h_evals;
				largest_evals = evnew;
				largest_evecs.insert(largest_evecs.end(), h_evecs.begin(), h_evecs.end());
			}
		}
		if(Nsmallest > 0)
		{
			//double lambda0 = largest_evals[0];
			//get smallest first, then get next smallest by shifting
			int step = 10;
			SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE> solver(op);
			smallest_evals = Eigen::VectorXd::Ones(step+1);
			while(std::fabs(smallest_evals[smallest_evals.size()-1] - smallest_evals[smallest_evals.size()-2]) < 1E-06*std::fabs(smallest_evals[smallest_evals.size()-1]))
			{
				solver.calc_smallest_evalues(step, smallest_evals, smallest_evecs); 
				if(std::fabs(smallest_evals[smallest_evals.size()-1] - smallest_evals[smallest_evals.size()-2]) < 1E-06*std::fabs(smallest_evals[smallest_evals.size()-1]))
				{
					step += 1;
					//std::cout << step << std::endl;
				}
			}
			while(smallest_evals.size() < Nsmallest)
			{
				std::cout << smallest_evals.size() << "\n";
				//start again between last two
				double lambda0 = 0.5*(smallest_evals[0] + smallest_evals[1]);
				//remove first
				smallest_evals = smallest_evals.tail(smallest_evals.size()-1);
				smallest_evecs.erase(smallest_evecs.begin());
				step = 10;
				if(step > int(Nsmallest - smallest_evals.size()) - 1) 
				{
					step = int(Nsmallest - smallest_evals.size()) - 1;
				}
				Eigen::VectorXd h_evals = Eigen::VectorXd::Ones(step+1);
				std::vector<Eigen::VectorXd> h_evecs;
				while(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1])
					  && step < int(mat_dim - smallest_evals.size()))
				{
					if(solver.shift_solve_largest_evalues(step+1, h_evals, h_evecs, lambda0))
					{
						std::cout << "Error computing spectrum, aborting...\n";
						abort_on_failure();
					}
					if(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1]))
					{
						step += 11;
					}
				}
				Eigen::VectorXd evnew(smallest_evals.size() + h_evals.size());
				evnew << h_evals, smallest_evals;
				smallest_evals = evnew;
				smallest_evecs.insert(smallest_evecs.begin(), h_evecs.begin(), h_evecs.end());
			}
		}
	}
	//concatenate
	if(Nsmallest > 0 && Nlargest > 0)
	{
		evalues.head(Nlargest) = largest_evals.head(Nlargest);
		evalues.tail(Nsmallest) = smallest_evals.tail(Nsmallest);
		evectors.insert(evectors.end(), largest_evecs.begin(), largest_evecs.begin() + Nlargest);
		evectors.insert(evectors.end(), smallest_evecs.begin(), smallest_evecs.begin() + Nsmallest);
	}
	else
	{
		if(Nsmallest > 0)
		{
			evalues = smallest_evals.tail(Nsmallest);
			evectors.insert(evectors.end(), smallest_evecs.begin(), smallest_evecs.begin() + Nsmallest);
		}
		if(Nlargest > 0)
		{
			evalues = largest_evals;
			evectors.insert(evectors.end(), largest_evecs.begin(), largest_evecs.begin() + Nlargest);
		}
	}
	auto end = std::chrono::system_clock::now();
	std::cout << "Spectral calc took: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << "s.\n";
}

//compute all eigenvalues above a particular cutoff
template<class SpectraMatrixClass, class SpectraShiftSolveClass, int RTYPE> void compute_spectrum_above_cutoff
	                             (SpectraMatrixClass & op, const size_t & mat_dim, const double & cutoff,
								  Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	//auto start = std::chrono::system_clock::now();
	int step = std::min(int(mat_dim), 10);
	SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE> solver(op);
	evalues = Eigen::VectorXd::Ones(step+1);
	//in order to calc next lot, need to ensure we are not going to repeat modes, check last 2 are not degenerate
	while(std::fabs(evalues[evalues.size()-1] - evalues[evalues.size()-2]) < 1E-06*std::fabs(evalues[evalues.size()-1]))
	{
		solver.calc_largest_evalues(step+1,evalues,evectors);
		if(std::fabs(evalues[evalues.size()-1] - evalues[evalues.size()-2]) < 1E-06*std::fabs(evalues[evalues.size()-1]))
		{
			step += 1;
			//std::cout << step << std::endl;
		}
	}

	while(evalues[evalues.size()-1] > cutoff && evalues.size() < int(mat_dim))   //continue calculating until we hit cutoff
	{
		std::cout << evalues.size() << "\n";
		//start again between last two
		double lambda0 = 0.5*(evalues[evalues.size()-1] + evalues[evalues.size()-2]);
		//remove last
		evalues.conservativeResize(evalues.size()-1);
		evectors.pop_back();
		Eigen::VectorXd h_evals = Eigen::VectorXd::Ones(step+1);
		std::vector<Eigen::VectorXd> h_evecs;
		step = std::min(10, int(mat_dim - evalues.size()) - 1);
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
				std::cout << step << std::endl;
			}
		}
		Eigen::VectorXd evnew(evalues.size() + h_evals.size());
		evnew << evalues, h_evals;
		evalues = evnew;
		evectors.insert(evectors.end(), h_evecs.begin(), h_evecs.end());
	}
}

//compute all eigenvalues below a particular cutoff
template<class SpectraMatrixClass, class SpectraShiftSolveClass, int RTYPE> void compute_spectrum_below_cutoff
	                             (SpectraMatrixClass & op, const size_t & mat_dim, const double & cutoff,
								  Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	//auto start = std::chrono::system_clock::now();
	int step = std::min(int(mat_dim), 10);
	SymmSolver<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE> solver(op);
	evalues = Eigen::VectorXd::Ones(step+1);
	//in order to calc next lot, need to ensure we are not going to repeat modes, check last 2 are not degenerate
	while(std::fabs(evalues[evalues.size()-1] - evalues[evalues.size()-2]) < 1E-06*std::fabs(evalues[evalues.size()-1]))
	{
		solver.calc_smallest_evalues(step+1,evalues,evectors);
		if(std::fabs(evalues[0] - evalues[1]) < 1E-06*std::fabs(evalues[0]))
		{
			step += 1;
			//std::cout << step << std::endl;
		}
	}

	while(evalues[0] < cutoff && evalues.size() < int(mat_dim))   //continue calculating until we hit cutoff
	{
		std::cout << evalues.size() << "\n";
		//start again between last two
		double lambda0 = 0.5*(evalues[0] + evalues[1]);
		//remove first
		Eigen::VectorXd evalnew = evalues.tail(evalues.size()-1);
		evalues = evalnew;
		evectors.erase(evectors.begin());
		Eigen::VectorXd h_evals = Eigen::VectorXd::Ones(step+1);
		std::vector<Eigen::VectorXd> h_evecs;
		step = std::min(10, int(mat_dim - evalues.size()) - 1);
		while(std::fabs(h_evals[0] - h_evals[1]) < 1E-06*std::fabs(h_evals[0])
				  && step < int(mat_dim - evalues.size()))
		{
			if(solver.shift_solve_largest_evalues(step+1, h_evals, h_evecs, -lambda0))
			{
				std::cout << "Error computing spectrum, aborting...\n";
				abort_on_failure();
			}
			if(std::fabs(h_evals[h_evals.size()-1] - h_evals[h_evals.size()-2]) < 1E-06*std::fabs(h_evals[h_evals.size()-1]))
			{
				step += 1;
				//std::cout << step << std::endl;
			}
		}
		Eigen::VectorXd evnew(evalues.size() + h_evals.size());
		evnew << h_evals, evalues;
		evalues = evnew;
		evectors.insert(evectors.begin(), h_evecs.begin(), h_evecs.end());
	}
}

//compute whole spectrum of operator
template<class SpectraMatrixClass, class SpectraShiftSolveClass, int RTYPE> void compute_full_spectrum(SpectraMatrixClass & op, const size_t & mat_dim, 
													Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
{
	if(mat_dim <= 1000)
	{
		compute_partial_spectrum<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE>(op, mat_dim, int(0.51*mat_dim), 
			                                               int(mat_dim - size_t(0.51*mat_dim)), evalues, evectors);
	}
	else
	{
		compute_partial_spectrum<SpectraMatrixClass, SpectraShiftSolveClass, RTYPE>(op, mat_dim, 0, int(mat_dim), 
			                                                                        evalues, evectors);
	}
}


	
//template specifiers
template<int RTYPE> class SymmSolver<Rmat<RTYPE>, Rmatinv<RTYPE>, RTYPE>    //functions defined for SpectraMatrix class is a r matrix
{
protected:
	std::shared_ptr<Rmat<RTYPE>> matrix, fs_mat;
	std::shared_ptr<Rmatinv<RTYPE>> fs_solve, shift_solve;
	void initialise(const Rmat<RTYPE> & mat)
	{
		this->matrix = std::make_shared<Rmat<RTYPE>>(mat);
		Eigen::VectorXd conductances = Eigen::VectorXd::Ones(this->matrix->return_rvec().size()).array() 
										/ this->matrix->return_rvec().array();
		this->fs_mat =  std::make_shared<Rmat<RTYPE>>(this->matrix->return_tree_pointer(), -conductances);
		this->fs_solve = std::make_shared<Rmatinv<RTYPE>>(this->matrix->return_tree_pointer(), -conductances);
		this->shift_solve = std::make_shared<Rmatinv<RTYPE>>(this->matrix->return_tree_pointer(), conductances);
	}
public:
	SymmSolver(const Rmat<RTYPE> & mat)
	{
		this->initialise(mat);
	}

	int calc_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
	{
		int eig_param = std::min(2*Neigs, int(this->matrix->rows()));
		int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Rmat<RTYPE>>, Rmat<RTYPE>>
			          (this->matrix.get(), Neigs, eig_param, int(this->matrix->rows()), evalues, evectors);
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

	int calc_smallest_evalues(const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
	{
		int eig_param = std::min(2*Neigs, int(this->matrix->rows()));
		int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Rmat<RTYPE>>, Rmat<RTYPE>>
			          (this->matrix.get(), Neigs, eig_param, int(this->matrix->rows()), evalues, evectors);
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

	int shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, 
		                                   std::vector<Eigen::VectorXd> & evectors, const double & s)
	{
		int eig_param = std::min(5*Neigs, int(this->matrix->rows()));
		int flag = do_shift_solve<Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_ALGE, Rmatinv<RTYPE>>, 
			Rmatinv<RTYPE>>(this->shift_solve.get(), Neigs, eig_param, int(this->matrix->rows()), s, evalues, evectors);
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

	int flip_shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, 
		                                   std::vector<Eigen::VectorXd> & evectors, const double & s)
	{
		int eig_param = std::min<int>(5*Neigs, this->matrix->rows());
		int flag = do_shift_solve<Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_ALGE, Rmatinv<RTYPE>>, 
			Rmatinv<RTYPE>>(this->fs_solve.get(), Neigs, eig_param, int(this->matrix->rows()), s, evalues, evectors);
		//flip back
		evalues = -evalues;
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}
};

template<int RTYPE> class SymmSolver<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>    
//functions defined for SpectraMatrix class is Eigen Sparse double
{
protected:
	Eigen::SparseMatrix<double> matrix, fs_mat;
	std::shared_ptr<Spectra::SparseSymShiftSolve<double>> fs_solve;
	std::shared_ptr<Spectra::SparseSymShiftSolve<double>> shift_solve;
	void initialise(const Eigen::SparseMatrix<double> & mat)
	{
		this->matrix = mat;
		this->fs_mat = -mat;
		this->fs_solve = std::make_shared<Spectra::SparseSymShiftSolve<double>>(this->fs_mat);
		this->shift_solve = std::make_shared<Spectra::SparseSymShiftSolve<double>>(this->matrix);
	}
public:
	SymmSolver(const Eigen::SparseMatrix<double> & mat)
	{
		this->initialise(mat);
	}

	int calc_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
	{
		int eig_param = std::min(2*Neigs, int(this->matrix.rows()));
		std::shared_ptr<Spectra::SparseSymMatProd<double>> prod_op
			              = std::make_shared<Spectra::SparseSymMatProd<double>>(this->matrix);
		int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseSymMatProd<double>>, 
					  Spectra::SparseSymMatProd<double>>(prod_op.get(), Neigs, eig_param, int(this->matrix.rows()), evalues, evectors);
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

	int calc_smallest_evalues(const int & Neigs, Eigen::VectorXd & evalues, std::vector<Eigen::VectorXd> & evectors)
	{
		int eig_param = std::min(2*Neigs, int(this->matrix.rows()));
		std::shared_ptr<Spectra::SparseSymMatProd<double>> prod_op 
			                   = std::make_shared<Spectra::SparseSymMatProd<double>>(this->matrix);
		int flag = do_solve<Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double>>, 
									Spectra::SparseSymMatProd<double>>
						(prod_op.get(), Neigs, eig_param, int(this->matrix.rows()), evalues, evectors);
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

	int shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, 
		                                   std::vector<Eigen::VectorXd> & evectors, const double & s)
	{
		int eig_param = std::min(5*Neigs, int(this->matrix.rows()));
		this->shift_solve->set_shift(s);
		int flag = do_shift_solve<Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseSymShiftSolve<double>>, 
					   Spectra::SparseSymShiftSolve<double>>(this->shift_solve.get(), Neigs, eig_param, int(this->matrix.rows()), s, evalues, evectors);
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

	int flip_shift_solve_largest_evalues(const int & Neigs, Eigen::VectorXd & evalues, 
		                                   std::vector<Eigen::VectorXd> & evectors, const double & s)
	{
		int eig_param = std::min(5*Neigs, int(this->matrix.rows()));
		this->fs_solve->set_shift(s);
		int flag = do_shift_solve<Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_ALGE, Spectra::SparseSymShiftSolve<double>>, 
					   Spectra::SparseSymShiftSolve<double>>(this->fs_solve.get(), Neigs, eig_param, int(this->matrix.rows()), s, evalues, evectors);
		//flip back
		evalues = -evalues;
		sort_and_return_by_evalues(evalues, evectors);
		return flag;
	}

};

//function to caluclate edge weights -- default poiseille resistance
template<int TYPE> double edge_weight_calc(const double & rad, const double & length, const double & scale)
{
	return (scale*rad*rad*rad*rad/length);
}

//or 1D diffusion
template<> double edge_weight_calc<DIFFUSION_NETWORK>(const double & rad, const double & length, const double & scale);

//custom network class
template<int RTYPE = RMAT1> class SpectralNetwork: public Tree
{
private:
	typedef network::Node Snode;
	typedef network::Edge<network::Node> Sedge;

	double Ceff, P_int_sqrd, Q_ext_sqrd, Lap_11;

	Eigen::VectorXd graph_laplacian_evalues, truncated_laplacian_evalues, modified_laplacian_evalues,
		            maury_evalues, adjacency_evalues, truncated_laplacian_reffs, maury_ceffs, 
					modified_laplacian_ceffs, modified_laplacian_ceffs_alt, 
					modified_laplacian_sqrd_P_int_contribution, maury_sqrd_Q_ext_contribution;
	std::vector<Eigen::VectorXd> graph_laplacian_evectors, truncated_laplacian_evectors, 
		            modified_laplacian_evectors,  maury_evectors, adjacency_evectors;
	Eigen::VectorXd edge_weight, flux, pressure, Nt, laplacian_flux, laplacian_pressure, 
		maury_flux, maury_pressure, modified_laplacian_flux, modified_laplacian_pressure, Lap_avec;

	Eigen::SparseMatrix<double> Degree, Adjacency, Laplacian, Truncated_Laplacian, Modified_Laplacian,
		                        Lap_OffDiag, Lap_Ext;
	Eigen::SparseMatrix<int> Incidence;

	double maury_scale_factor;
	Rmat<RTYPE> Maury_matrix;
	int tree_type;   //either RESISTANCE_NETWORK or DIFFUSION_NETWORK

	void build_asymmetric_network(const size_t & Ngens, const double & rad0, const double & length0,
			                          const double & scale_factor, const double & asymm_factor);
	void fill_matrices();
	void solve_effective_conductance_problem();
	void fill_fluxes_from_term_vals(Eigen::VectorXd & flux, const Eigen::VectorXd & term_flux);
	void get_term_count_edges(Eigen::VectorXd & term_count);
	void calc_truncated_laplacian_reffs();
	void calc_modified_laplacian_ceffs();
	void calc_maury_ceffs();
	void print_evectors_vtk(const std::string & filename, const size_t & Nvecs, const Eigen::VectorXd & evalues,
							const Eigen::VectorXd & contributions, const std::vector<Eigen::VectorXd> & evectors,
							const std::vector<Eigen::VectorXd> & solution_contributions, const int & sort_option);
	void print_evectors_csv(const std::string & filename, const size_t & Nvecs, const Eigen::VectorXd & evalues,
									const Eigen::VectorXd & contributions, const std::vector<Eigen::VectorXd> & evectors,
									const int & sort_option);
public:
	SpectralNetwork():Tree(){};
	
	SpectralNetwork(const std::vector<std::shared_ptr<Snode>> & n, 
		            const std::vector<std::shared_ptr<Sedge>> & e, const int & type):Tree(n,e)
	{
		this->tree_type = type;
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}

	SpectralNetwork(const std::map<long int, std::shared_ptr<Snode>> & n, 
		            const std::map<long int, std::shared_ptr<Sedge>> & e,
					const int & type):Tree(n,e)
	{
		this->tree_type = type;
		std::cout << "Filling matrices..." << std::endl;
		this->fill_matrices();
		std::cout << "Solving flow problem..." << std::endl;
		this->solve_effective_conductance_problem();
	}

	SpectralNetwork(const std::string & node_fname, const std::string & edge_fname, const std::string & term_fname, const int & type):Tree(node_fname, edge_fname, term_fname, 1.0)
	{
		this->tree_type = type;
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}

	SpectralNetwork(const size_t & Ngens, const double & rad0, const double & length0,
			        const double & scale_factor, const double & asymm_factor, const int & type)
	{
		this->tree_type = type;
		this->build_asymmetric_network(Ngens, rad0, length0, scale_factor, asymm_factor);
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}

	SpectralNetwork(const size_t & Ngens, const double & rad0, const double & length0,
			        const double & scale_factor, const double & asymm_factor, const double & dr,
					const size_t & pert_gen, const int & type)
	{
		this->tree_type = type;
		this->build_asymmetric_network(Ngens, rad0, length0, scale_factor, asymm_factor);
		double new_rad = this->get_edge(size_t(pow(2,pert_gen))-1)->get_geom()->get_inner_radius()*pow(1+dr,-1.0/4.0);
		this->get_edge(size_t(pow(2,pert_gen))-1)->update_inner_radius(new_rad);
		this->fill_matrices();
		this->solve_effective_conductance_problem();
	}

	inline const Eigen::VectorXd& get_edge_weights(){ return this->edge_weight; }

	inline void compute_partial_graph_laplacian_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		size_t mat_dim = size_t(this->Laplacian.rows());
		compute_partial_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			                    (this->Laplacian, mat_dim, Nsmallest, Nlargest,
								 this->graph_laplacian_evalues, this->graph_laplacian_evectors);
	}
	inline void compute_partial_truncated_laplacian_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		size_t mat_dim = this->Truncated_Laplacian.rows();
		compute_partial_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			                                 (this->Truncated_Laplacian, mat_dim, Nsmallest, Nlargest,
								                     this->truncated_laplacian_evalues, this->truncated_laplacian_evectors);
		this->calc_truncated_laplacian_reffs();
	}
	inline void compute_partial_modified_laplacian_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		size_t mat_dim = this->Modified_Laplacian.rows();
		compute_partial_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			                            (this->Modified_Laplacian, mat_dim, Nsmallest, Nlargest,
								         this->modified_laplacian_evalues, this->modified_laplacian_evectors);
		this->calc_modified_laplacian_ceffs();
	}
	inline void compute_partial_maury_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		size_t mat_dim = this->Maury_matrix.rows();
		compute_partial_spectrum<Rmat<RTYPE>,Rmatinv<RTYPE>,RTYPE>(this->Maury_matrix, mat_dim, Nsmallest, Nlargest,
			                                  this->maury_evalues, this->maury_evectors);
		//do ceff calc
		this->maury_evalues = this->maury_scale_factor*this->maury_evalues;
		this->calc_maury_ceffs();
	}
	inline void compute_partial_adjacency_spectrum(const int & Nsmallest, const int & Nlargest)
	{
		compute_partial_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			(this->Adjacency, this->count_nodes(), Nsmallest, Nlargest, this->adjacency_evalues, this->adjacency_evectors);
	}

	inline void compute_full_graph_laplacian_spectrum()
	{
		size_t mat_dim = this->Laplacian.rows();
		compute_full_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Laplacian, mat_dim,  this->graph_laplacian_evalues, this->graph_laplacian_evectors);
	}
	inline void compute_full_truncated_laplacian_spectrum()
	{
		size_t mat_dim = this->Truncated_Laplacian.rows();
		compute_full_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Truncated_Laplacian, mat_dim,  this->truncated_laplacian_evalues, this->truncated_laplacian_evectors);
		this->calc_truncated_laplacian_reffs();
	}
	inline void compute_full_modified_laplacian_spectrum()
	{
		size_t mat_dim = this->Modified_Laplacian.rows();
		compute_full_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Modified_Laplacian, mat_dim,  this->modified_laplacian_evalues, this->modified_laplacian_evectors);
		this->calc_modified_laplacian_ceffs();
	}
	inline void compute_full_maury_spectrum()
	{
		size_t mat_dim = this->Maury_matrix.rows();
		compute_full_spectrum<Rmat<RTYPE>,Rmatinv<RTYPE>, RTYPE>(this->Maury_matrix, mat_dim, this->maury_evalues, 
			                                                     this->maury_evectors);
		//do ceff calc
		this->maury_evalues = this->maury_scale_factor*this->maury_evalues;
		this->calc_maury_ceffs();
	}
	inline void compute_full_adjacency_spectrum()
	{
		compute_full_spectrum<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			(this->Adjacency, this->count_nodes(), this->adjacency_evalues, this->adjacency_evectors);
	}

	inline void compute_all_graph_laplacian_modes_above_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Laplacian.rows();
		compute_spectrum_above_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Laplacian, mat_dim, cutoff, this->graph_laplacian_evalues,
					  this->graph_laplacian_evectors);
	}
	inline void compute_all_truncated_laplacian_modes_above_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Truncated_Laplacian.rows();
		compute_spectrum_above_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Truncated_Laplacian, mat_dim, cutoff, this->truncated_laplacian_evalues,
					  this->truncated_laplacian_evectors);
		this->calc_truncated_laplacian_reffs();
	}
	inline void compute_all_modified_laplacian_modes_above_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Modified_Laplacian.rows();
		compute_spectrum_above_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Modified_Laplacian, mat_dim, cutoff, this->modified_laplacian_evalues,
					  this->modified_laplacian_evectors);
		this->calc_modified_laplacian_ceffs();
	}
	inline void compute_all_maury_modes_above_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Maury_matrix.rows();
		compute_spectrum_above_cutoff<Rmat<RTYPE>,Rmatinv<RTYPE>,RTYPE>(this->Maury_matrix, mat_dim, cutoff,
			                                                      this->maury_evalues, this->maury_evectors);
		//do ceff calc
		this->maury_evalues = this->maury_scale_factor*this->maury_evalues;
		this->calc_maury_ceffs();
	}
	inline void compute_all_adjacency_modes_above_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Adjacency.rows();
		compute_spectrum_above_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Adjacency, mat_dim, cutoff, this->adjacency_evalues,
					  this->adjacency_evectors);
	}

	inline void compute_all_graph_laplacian_modes_below_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Laplacian.rows();
		compute_spectrum_below_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Laplacian, mat_dim, cutoff, this->graph_laplacian_evalues,
					  this->graph_laplacian_evectors);
	}
	inline void compute_all_maury_modes_below_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Maury_matrix.rows();
		compute_spectrum_below_cutoff<Rmat<RTYPE>,Rmatinv<RTYPE>, RTYPE>(this->Maury_matrix, mat_dim, cutoff,
			                                                      this->maury_evalues, this->maury_evectors);
		//do ceff calc
		this->maury_evalues = this->maury_scale_factor*this->maury_evalues;
		this->calc_maury_ceffs();
	}
	inline void compute_all_truncated_laplacian_modes_below_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Truncated_Laplacian.rows();
		compute_spectrum_below_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Truncated_Laplacian, mat_dim, cutoff, this->truncated_laplacian_evalues,
					  this->truncated_laplacian_evectors);
		this->calc_truncated_laplacian_reffs();
	}
	inline void compute_all_modified_laplacian_modes_below_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Modified_Laplacian.rows();
		compute_spectrum_below_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Modified_Laplacian, mat_dim, cutoff, this->modified_laplacian_evalues,
					  this->modified_laplacian_evectors);
		this->calc_modified_laplacian_ceffs();
	}
	inline void compute_all_adjacency_modes_below_cutoff(const double & cutoff)
	{
		size_t mat_dim = this->Adjacency.rows();
		compute_spectrum_below_cutoff<Eigen::SparseMatrix<double>, Spectra::SparseSymShiftSolve<double>, RTYPE>
			         (this->Adjacency, mat_dim, cutoff, this->adjacency_evalues,
					  this->adjacency_evectors);
	}

	void print_graph_laplacian_evectors_vtk(const std::string & filename, const size_t & Nvecs, const int & sort_option);
	inline void print_graph_laplacian_evectors_csv(const std::string & filename, const size_t & Nvecs, const int & sort_option)
	{
		//convert maury eigenvectors into edge values for visualisation purposes
		Eigen::VectorXd blank = Eigen::VectorXd::Zero(this->graph_laplacian_evalues.size());
		this->print_evectors_csv(filename, Nvecs, this->graph_laplacian_evalues, blank, 
			                     this->truncated_laplacian_evectors, sort_option);
	}

	void print_truncated_laplacian_evectors_vtk(const std::string & filename, const size_t & Nvecs, const int & sort_option);
	inline void print_truncated_laplacian_evectors_csv(const std::string & filename, const size_t & Nvecs, const int & sort_option)
	{
		//convert maury eigenvectors into edge values for visualisation purposes
		this->print_evectors_csv(filename, Nvecs, this->truncated_laplacian_evalues, this->truncated_laplacian_reffs, 
			                     this->truncated_laplacian_evectors, sort_option);
	}

	void print_modified_laplacian_evectors_vtk(const std::string & filename, const size_t & Nvecs, const int & sort_option);
	inline void print_modified_laplacian_evectors_csv(const std::string & filename, const size_t & Nvecs, const int & sort_option)
	{
		//convert maury eigenvectors into edge values for visualisation purposes
		this->print_evectors_csv(filename, Nvecs, this->modified_laplacian_evalues, this->modified_laplacian_sqrd_P_int_contribution/this->P_int_sqrd, 
			                     this->modified_laplacian_evectors, sort_option);
	}

	void print_maury_evectors_vtk(const std::string & filename, const size_t & Nvecs, const int & sort_option);
	inline void print_maury_evectors_csv(const std::string & filename, const size_t & Nvecs, const int & sort_option)
	{
		//convert maury eigenvectors into edge values for visualisation purposes
		this->print_evectors_csv(filename, Nvecs, this->maury_evalues, this->maury_sqrd_Q_ext_contribution/this->Q_ext_sqrd, 
			                     this->maury_evectors, sort_option);
	}

	void print_adjacency_evectors_vtk(const std::string & filename, const size_t & Nvecs, const int & sort_option);
	inline void print_adjacency_evectors_csv(const std::string & filename, const size_t & Nvecs, const int & sort_option)
	{
		//convert maury eigenvectors into edge values for visualisation purposes
		Eigen::VectorXd blank = Eigen::VectorXd::Zero(this->adjacency_evalues.size());
		this->print_evectors_csv(filename, Nvecs, this->adjacency_evalues, blank, 
			                     this->adjacency_evectors, sort_option);
	}

	void print_unit_pressure_drop_solutions_vtk(const std::string & filename);

	void print_graph_laplacian_modes_summary_csv(const std::string & filename);
	void print_truncated_laplacian_modes_summary_csv(const std::string & filename);
	void print_modified_laplacian_modes_summary_csv(const std::string & filename);
	void print_maury_modes_summary_csv(const std::string & filename);
	void print_adjacency_modes_summary_csv(const std::string & filename);

	inline double get_maury_scale_factor(){ return this->maury_scale_factor; }
	inline void get_incidence_matrix(Eigen::SparseMatrix<double> & incidence){ incidence = this->Incidence; }
};

//Spectral network functions
template<int RTYPE> void SpectralNetwork<RTYPE>::build_asymmetric_network(const size_t & Ngens, const double & rad0, const double & length0,
			                          const double & scale_factor, const double & asymm_factor)
{
	std::vector<std::shared_ptr<network::Node>> nodes;
	std::vector<std::shared_ptr<network::Edge<network::Node>>> edges;
	std::vector<double> theta_vals;
	nodes.resize(size_t(pow(2,Ngens)));
	theta_vals.resize(nodes.size());
	edges.resize(nodes.size()-1);

	nodes[0] = std::make_shared<network::Node>(0,0,length0);
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
			nodes[node_count] = std::make_shared<network::Node>(nodes[node_in]->get_pos() + dir*length);
			edges[node_count-1] = std::make_shared<network::Edge<network::Node>>(nodes[node_in], nodes[node_count], 1, radius);
			node_count++;
		}
	}
	Tree(nodes,edges);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::fill_matrices()
{
	std::cout << "Filling matrices...\n";
	auto start = std::chrono::system_clock::now();
	//fill vector of edge weights
	this->edge_weight = Eigen::VectorXd::Zero(this->count_edges());
	this->maury_scale_factor = 0;
	for(size_t j = 0; j < this->count_edges(); j++)
	{
		double scale =  this->get_edge(j)->branch_count();
		if(this->tree_type == RESISTANCE_NETWORK)
		{
			scale *= M_PI*1E-06 / (8*VISCOSITY);   //conductance in L / cmH20 s, airway geom in mm
			this->edge_weight[j] = edge_weight_calc<RESISTANCE_NETWORK>(this->get_edge(j)->get_geom()->get_inner_radius(),
				                                      this->get_edge(j)->get_geom()->get_length(), scale);
		}
		if(this->tree_type == DIFFUSION_NETWORK)
		{
			scale *= M_PI * DIFFUSIVITY;
			this->edge_weight[j] = edge_weight_calc<DIFFUSION_NETWORK>(this->get_edge(j)->get_geom()->get_inner_radius(),
				                                      this->get_edge(j)->get_geom()->get_length(), scale);
		}
		if(1.0/this->edge_weight[j] > this->maury_scale_factor) this->maury_scale_factor = 1.0/this->edge_weight[j];
	}
	//initialise various matrices
	this->Degree = Eigen::SparseMatrix<double>(this->count_nodes(), this->count_nodes());
	this->Adjacency = Eigen::SparseMatrix<double>(this->count_nodes(), this->count_nodes());
	this->Incidence = Eigen::SparseMatrix<int>(this->count_nodes(), this->count_edges());
	std::vector<Eigen::Triplet<double>> degree_fill, adjacency_fill;
	std::vector<Eigen::Triplet<int>> incidence_fill;
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
			incidence_fill.push_back(Eigen::Triplet<int>(int(k), int(j), -1)); //incidence matrix

			//adjacency contribution
			size_t ki = this->get_node_in_index(j);
			adjacency_fill.push_back(Eigen::Triplet<double>(int(k), int(ki), this->edge_weight[j]));

			deg += this->edge_weight[j];    //count degree
		}
		//loop over edges out
		for(size_t jo = 0; jo < this->count_edges_out(k); jo++)
		{
			size_t j = this->get_edge_out_index(k,jo);
			incidence_fill.push_back(Eigen::Triplet<int>(int(k), int(j), 1));  //incidence matrix

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
	////*******Fill modified laplacian, just remove first row and col***************//
	this->Lap_11 = this->Laplacian.coeffRef(0,0);
	this->Lap_avec = Eigen::VectorXd::Zero(this->count_nodes() - this->count_term_nodes() - 1);
	this->Modified_Laplacian = Eigen::SparseMatrix<double>(this->count_nodes() - this->count_term_nodes() - 1,
		                                                   this->count_nodes() - this->count_term_nodes() - 1);
	this->Lap_OffDiag = Eigen::SparseMatrix<double>(this->count_nodes() - this->count_term_nodes() - 1,
		                                            this->count_term_nodes());
	this->Lap_Ext = Eigen::SparseMatrix<double>(this->count_term_nodes(), this->count_term_nodes());
	std::vector<Eigen::Triplet<double>> ML_fill, LOD_fill, LE_fill;
	ML_fill.reserve(4*(this->count_nodes()-this->count_term_nodes()-1));
	LOD_fill.reserve(2*(this->count_term_nodes()));
	LE_fill.reserve(this->count_term_nodes());
	//loop over first column
	for (Eigen::SparseMatrix<double>::InnerIterator it(this->Laplacian,0); it; ++it)
	{
		if(it.row() > 0 && it.row() < int(this->get_first_term_index()))  //get inner part of laplacian
		{
			this->Lap_avec(int(it.row()-1)) = it.value();
		}
	}

	for(size_t k = 1; k < this->get_first_term_index(); k++)  //ignore first node -- col ordering
	{
		//loop over rows in each col k
		for (Eigen::SparseMatrix<double>::InnerIterator it(this->Laplacian,k); it; ++it)
		{
			if(it.row() > 0 && it.row() < int(this->get_first_term_index()))  //get inner part of laplacian
			{
				ML_fill.push_back(Eigen::Triplet<double>(int(it.row()-1), int(it.col()-1), it.value()));
			}
		}
	}
	this->Modified_Laplacian.setFromTriplets(ML_fill.begin(),ML_fill.end());

	for(size_t k = this->get_first_term_index(); k < this->count_nodes(); k++)  //loop over term nodes
	{
		//loop over cols
		for (Eigen::SparseMatrix<double>::InnerIterator it(this->Laplacian,k); it; ++it)
		{
			if(it.row() > 0 && it.row() < int(this->get_first_term_index()))  //get off diagonal part
			{
				LOD_fill.push_back(Eigen::Triplet<double>(int(it.row()-1), 
					               int(it.col()-this->get_first_term_index()), it.value()));
			}
			if(it.row() >= int(this->get_first_term_index()))
			{
				LE_fill.push_back(Eigen::Triplet<double>(int(it.row()-this->get_first_term_index()), 
					               int(it.col()-this->get_first_term_index()), it.value()));
			}
		}
	}
	this->Lap_OffDiag.setFromTriplets(LOD_fill.begin(), LOD_fill.end());
	this->Lap_Ext.setFromTriplets(LE_fill.begin(), LE_fill.end());

	this->Degree.makeCompressed();
	this->Adjacency.makeCompressed();
	this->Laplacian.makeCompressed();
	this->Truncated_Laplacian.makeCompressed();
	this->Modified_Laplacian.makeCompressed();
	this->Lap_OffDiag.makeCompressed();
	this->Lap_Ext.makeCompressed();

	//std::cout << this->Truncated_Laplacian.toDense() << "\n\n";

	//this->Laplacian = this->Incidence * Eigen::DiagonalMatrix(this->edge_weight) *  this->Incidence.transpose();
	Eigen::VectorXd Maury_scaled_edge_weights = this->edge_weight;
	for(int j = 0; j < this->edge_weight.size(); j++)
	{
		Maury_scaled_edge_weights[j] *= this->maury_scale_factor;
	}

	this->Maury_matrix = Rmat<RTYPE>(this, Maury_scaled_edge_weights);
	this->fill_fluxes_from_term_vals(this->Nt, Eigen::VectorXd::Ones(this->count_term_nodes()));
	auto end = std::chrono::system_clock::now();
	std::cout << "Matrix fill took " << (std::chrono::duration<double>(end-start)).count() << '\n';
}

template<int RTYPE> void SpectralNetwork<RTYPE>::solve_effective_conductance_problem()
{
	//define matrix
	std::cout << "Solving pressure drop case...\n";
	auto start = std::chrono::system_clock::now();
	//cout << "Populating matrix.\n";
	Eigen::SparseLU<Eigen::SparseMatrix<double>> LUsolve;
	//left hand side is L_int
	LUsolve.compute(this->Modified_Laplacian);
	//RHS is -1(L_OD * -1)
	Eigen::VectorXd e = Eigen::VectorXd::Ones(this->count_term_nodes());
	Eigen::VectorXd Bvec = this->Lap_OffDiag*e;
	//cout << LUsolve.lastErrorMessage() << '\n';
	//cout << "Solving.\n";
	Eigen::VectorXd P_int = LUsolve.solve(Bvec);
	Eigen::VectorXd Q_ext = this->Lap_OffDiag.transpose()*P_int
		                  - this->Lap_Ext*e;
	//cout << "Solved.\n";
	auto end = std::chrono::system_clock::now();
	//extract interior fluxes
	this->fill_fluxes_from_term_vals(this->flux, -Q_ext);
	this->P_int_sqrd = P_int.squaredNorm();
	this->Q_ext_sqrd = Q_ext.squaredNorm();
	this->Ceff = -Q_ext.sum();
	std::cout << "Pressure drop solve took " << (std::chrono::duration<double>(end-start)).count() << '\n';
}

template<int RTYPE> void SpectralNetwork<RTYPE>::fill_fluxes_from_term_vals(Eigen::VectorXd & flux, const Eigen::VectorXd & term_flux)
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

template<int RTYPE> void SpectralNetwork<RTYPE>::get_term_count_edges(Eigen::VectorXd & term_count)
{
	term_count = Eigen::VectorXd::Zero(this->count_edges());
	for(size_t k = this->get_first_term_index(); k < this->count_nodes(); k++)
	{
		//size_t kt = k - this->get_first_term_index();
		size_t j = this->get_edge_in_index(k,0);
		term_count[j] = 1;
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
				term_count[j] += term_count[eo_index];
			}
		}
	}
	

}

template<int RTYPE> void SpectralNetwork<RTYPE>::calc_truncated_laplacian_reffs()
{
	std::cout << "Computing effective resistance contributions" << std::endl;
	this->truncated_laplacian_reffs = Eigen::VectorXd::Zero(this->truncated_laplacian_evalues.size());
	size_t TLsize = this->count_nodes() - this->count_term_nodes();
	for(size_t n = 0; n < size_t(this->truncated_laplacian_evalues.size()); n++)
	{
		if(this->truncated_laplacian_evalues[n] > 0)
		{
			this->truncated_laplacian_reffs[n] = (this->truncated_laplacian_evectors[n][0] - this->truncated_laplacian_evectors[n][TLsize])*
			                                 (this->truncated_laplacian_evectors[n][0] - this->truncated_laplacian_evectors[n][TLsize])/
								              this->truncated_laplacian_evalues[n];
		}
		else
		{
			this->truncated_laplacian_reffs[n] = 0;
		}
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

template<int RTYPE> void SpectralNetwork<RTYPE>::calc_modified_laplacian_ceffs()
{
	std::cout << "Computing effective conductance contributions" << std::endl;
	this->modified_laplacian_ceffs = Eigen::VectorXd::Zero(this->modified_laplacian_evalues.size());
	this->modified_laplacian_ceffs_alt = Eigen::VectorXd::Zero(this->modified_laplacian_evalues.size());
	this->modified_laplacian_sqrd_P_int_contribution = Eigen::VectorXd::Zero(this->modified_laplacian_evalues.size());
	for(size_t n = 0; n < size_t(this->modified_laplacian_evalues.size()); n++)   //loop over evalues
	{
		Eigen::VectorXd uC =  this->modified_laplacian_evectors[n].transpose()*this->Lap_OffDiag;
		double factor = uC.sum();
		this->modified_laplacian_ceffs[n] = -factor*factor/this->modified_laplacian_evalues[n];
		double factor_alt = this->Lap_avec.dot(this->modified_laplacian_evectors[n]);
		this->modified_laplacian_ceffs_alt[n] = -factor_alt*factor_alt/this->modified_laplacian_evalues[n];
		this->modified_laplacian_sqrd_P_int_contribution[n] = factor*factor
			                   / (this->modified_laplacian_evalues[n]*this->modified_laplacian_evalues[n]);
	}
	Eigen::VectorXd e = Eigen::VectorXd::Ones(this->count_term_nodes());
	Eigen::VectorXd Lextsum = this->Lap_Ext*e;
	//double Q = Lextsum.sum() + this->modified_laplacian_ceffs.sum();

	//compute the approximations to pressure and flux
	this->modified_laplacian_pressure = Eigen::VectorXd::Zero(this->count_nodes());
	this->modified_laplacian_flux = Eigen::VectorXd::Zero(this->count_edges());
	Eigen::VectorXd terminal_fluxes = Eigen::VectorXd::Zero(this->count_term_nodes());
	for(size_t n = 0; n < size_t(this->modified_laplacian_evalues.size()); n++)
	{
		Eigen::VectorXd u = this->modified_laplacian_evectors[n];
		Eigen::VectorXd uC =  u.transpose()*this->Lap_OffDiag;
		this->modified_laplacian_pressure.segment(1,u.size()) += (uC.sum()/this->modified_laplacian_evalues[n])*u;
		terminal_fluxes = Lextsum - (uC.sum()/this->modified_laplacian_evalues[n]) * (this->Lap_OffDiag.transpose()*u);
	}
	//fill in rest of fluxes by incompressibility
	this->fill_fluxes_from_term_vals(this->modified_laplacian_flux, terminal_fluxes);
	//fill in other pressures
	//P for node 0 Sum(P_0/r_j) - Sum(P_i/r_j) = Sum(q_j)
	double qsum = 0, Prsum = 0, rsum = 0;
	for(size_t jo = 0; jo < this->count_edges_out(0); jo++) 
	{
		size_t j = this->get_edge_out_index(0,jo);
		size_t kout = this->get_node_out_index(j);
		qsum += this->modified_laplacian_flux[j];
		Prsum += this->modified_laplacian_pressure[kout] / this->edge_weight[j];
		rsum += 1.0 / this->edge_weight[j];
	}
	this->modified_laplacian_pressure[0] = (qsum + Prsum)/rsum;
	//P for terminal nodes
	//Pk = Pkin - rj*qj
	for(size_t k = this->get_first_term_index(); k < this->count_nodes(); k++)
	{
		size_t j = this->get_edge_in_index(k,0);
		size_t kin = this->get_node_in_index(j);
		this->modified_laplacian_pressure[k] = this->modified_laplacian_pressure[kin] 
		                      - this->edge_weight[j]*this->modified_laplacian_flux[j];
	}
}

template<int RTYPE> void SpectralNetwork<RTYPE>::calc_maury_ceffs()
{
	std::cout << "Computing effective conductance contributions" << std::endl;
	this->maury_ceffs = Eigen::VectorXd::Zero(this->maury_evalues.size());
	this->maury_sqrd_Q_ext_contribution = Eigen::VectorXd::Zero(this->maury_evalues.size());
	Eigen::VectorXd term_fluxes = Eigen::VectorXd::Zero(this->count_term_nodes());
	for(size_t n = 0; n < size_t(this->maury_evalues.size()); n++)
	{
		double vec_sum = this->maury_evectors[n].sum();
		this->maury_ceffs[n] = (vec_sum *vec_sum)/(this->maury_evalues[n]);
		this->maury_sqrd_Q_ext_contribution[n] = (vec_sum *vec_sum)/(this->maury_evalues[n]*this->maury_evalues[n]);
		term_fluxes = term_fluxes + (vec_sum/(this->maury_evalues[n]))*this->maury_evectors[n];
	}
	//double Q = this->maury_ceffs.sum();
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

template<int RTYPE> void SpectralNetwork<RTYPE>::print_evectors_vtk(const std::string & filename, const size_t & Nvecs, const Eigen::VectorXd & evalues,
																  const Eigen::VectorXd & contributions, const std::vector<Eigen::VectorXd> & evectors,
																  const std::vector<Eigen::VectorXd> & solution_contributions, const int & sort_option)
{
	std::unordered_map<std::string, std::vector<double>> extra_vals;
	size_t to_print = std::min(Nvecs, size_t(evectors.size()));
	for(size_t n = 0; n < to_print; n++)
	{
		std::stringstream name;
		switch(sort_option)
		{
		case DOMINANT_SORT:
			{
				name << "Dominant_mode_";
			} break;

		case LARGEST_SORT:
			{
				name << "Largest_mode_";
			} break;

		case SMALLEST_SORT:
			{
				name << "Smallest_mode_";
			} break;

		default:
			{
				name << "Mode_";
			} break;
		}
		name << n << "_eigenvalue_" << evalues[n];

		extra_vals[name.str()] = std::vector<double>(evectors[n].data(), evectors[n].data() + evectors[n].size());

		name.clear();
		name.str("");
		if(evectors[0].size() == this->count_nodes())
		{
			name << "Pressure_contribution_mode_"; 
		}
		else
		{
			if(evectors[0].size() == this->count_edges())
			{
				name << "Flux_contribution_mode_";
			}
			else
			{
				name << "Unknown_contribution_mode_";
			}
		}		
		name << n << "_eigenvalue_" << evalues[n];
		extra_vals[name.str()] = std::vector<double>(solution_contributions[n].data(), solution_contributions[n].data() + solution_contributions[n].size());
	}
	this->print_vtk(filename, 1.0, extra_vals);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_evectors_csv(const std::string & filename, const size_t & Nvecs, 
		const Eigen::VectorXd & evalues, const Eigen::VectorXd & contributions, const std::vector<Eigen::VectorXd> & evectors,
		const int & sort_option)
{
	std::vector<std::vector<double>> data;
	std::vector<std::string> headers;
	headers.resize(evectors[0].size()+2);
	headers[0]=std::string("Eigenvalue");
	headers[1]=std::string("Contribution");
	for(int j = 0; j < evectors[0].size(); j++)
	{
		headers[j+2]=std::string("Uentry");
	}
	size_t to_print = std::min(Nvecs, size_t(evalues.size()));
	std::vector<size_t> indices(size_t(evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	switch(sort_option)
	{
		case LARGEST_SORT:
			{
				std::sort(indices.begin(), indices.end(), [&evalues](size_t i1, size_t i2){ return evalues[i1] > evalues[i2]; });
			} break;
		case SMALLEST_SORT:
			{
				std::sort(indices.begin(), indices.end(), [&evalues](size_t i1, size_t i2){ return evalues[i1] < evalues[i2]; });
			} break;
		case DOMINANT_SORT:
			{
				std::sort(indices.begin(), indices.end(), [&contributions](size_t i1, size_t i2){ return contributions[i1] > contributions[i2]; });
			} break;
		default:
			break;
	}
	data.resize(to_print);
	for(size_t n = 0; n < to_print; n++)
	{
		data[n].resize(evectors[indices[n]].size()+2);
		data[n][0] = evalues[indices[n]];
		data[n][1] = contributions[indices[n]];
		for(int j = 0; j < evectors[indices[n]].size(); j++)
		{
			data[n][j+2] = evectors[indices[n]][j];
		}
	}
	write_csv_file(filename, headers, data);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_graph_laplacian_evectors_vtk(const std::string & filename, 
															 const size_t & Nvecs, const int & sort_option)
{
	size_t to_print = std::min(Nvecs, size_t(this->graph_laplacian_evalues.size()));
	std::vector<size_t> indices(size_t(this->graph_laplacian_evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	switch(sort_option)
	{
	case LARGEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->graph_laplacian_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	case SMALLEST_SORT:
	case DOMINANT_SORT:
		{
			Eigen::VectorXd sort_vec = this->graph_laplacian_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) < sort_vec(i2)); });
		} break;

	default:
		break;
	}

	//copy truncated eigenvectors to full network
	std::vector<Eigen::VectorXd> node_evecs;
	Eigen::VectorXd evalues = Eigen::VectorXd::Zero(to_print);
	Eigen::VectorXd blank = Eigen::VectorXd::Zero(to_print);
	std::vector<Eigen::VectorXd> all_blanks(to_print);
	node_evecs.resize(to_print);
	for(size_t n = 0; n < to_print; n++)
	{
		evalues[n] = this->graph_laplacian_evalues[indices[n]];
		node_evecs[n] = this->truncated_laplacian_evectors[indices[n]];
		all_blanks[n] = Eigen::VectorXd::Zero(node_evecs[n].size());
	}
	//print
	this->print_evectors_vtk(filename, Nvecs, evalues, blank, node_evecs, all_blanks, sort_option);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_truncated_laplacian_evectors_vtk(const std::string & filename, 
															 const size_t & Nvecs, const int & sort_option)
{
	size_t to_print = std::min(Nvecs, size_t(this->truncated_laplacian_evalues.size()));
	std::vector<size_t> indices(size_t(this->truncated_laplacian_evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	switch(sort_option)
	{
	case LARGEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->truncated_laplacian_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	case SMALLEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->truncated_laplacian_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) < sort_vec(i2)); });
		} break;

	case DOMINANT_SORT:
		{
			Eigen::VectorXd sort_vec = this->truncated_laplacian_reffs;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	default:
		break;
	}

	//copy truncated eigenvectors to full network
	std::vector<Eigen::VectorXd> node_evecs, pressure_contributions;
	Eigen::VectorXd evalues = Eigen::VectorXd::Zero(to_print), reffs = Eigen::VectorXd::Zero(to_print);
	node_evecs.resize(to_print);
	pressure_contributions.resize(to_print);
	for(size_t n = 0; n < to_print; n++)
	{
		evalues[n] = this->truncated_laplacian_evalues[indices[n]];
		reffs[n] = this->truncated_laplacian_reffs[indices[n]];

		node_evecs[n] = Eigen::VectorXd::Zero(this->count_nodes());
		pressure_contributions[n] = Eigen::VectorXd::Zero(this->count_nodes());
		//all internal nodes remain the same
		node_evecs[n].head(this->truncated_laplacian_evectors[indices[n]].size()) = this->truncated_laplacian_evectors[indices[n]];
		//all term nodes have same value
		//get flux contributions
		node_evecs[n].tail(this->count_nodes()-this->truncated_laplacian_evectors[indices[n]].size())
			= this->truncated_laplacian_evectors[indices[n]][this->truncated_laplacian_evectors[indices[n]].size()-1]
			* Eigen::VectorXd::Ones(this->count_nodes()-this->truncated_laplacian_evectors[indices[n]].size());
	}
	//print
	this->print_evectors_vtk(filename, Nvecs, evalues, reffs, node_evecs, pressure_contributions, sort_option);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_modified_laplacian_evectors_vtk(const std::string & filename, 
															const size_t & Nvecs, const int & sort_option)
{
		
	size_t to_print = std::min(Nvecs, size_t(this->modified_laplacian_evalues.size()));
	std::vector<size_t> indices(size_t(this->modified_laplacian_evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	switch(sort_option)
	{
	case LARGEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->modified_laplacian_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	case SMALLEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->modified_laplacian_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) < sort_vec(i2)); });
		} break;

	case DOMINANT_SORT:
		{
			Eigen::VectorXd sort_vec = this->modified_laplacian_sqrd_P_int_contribution;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) < sort_vec(i2)); });
		} break;

	default:
		break;
	}

	//copy modified eigenvectors to full network
	std::vector<Eigen::VectorXd> node_evecs, pressure_contributions;
	Eigen::VectorXd evalues = Eigen::VectorXd::Zero(to_print), ceffs = Eigen::VectorXd::Zero(to_print);
	Eigen::VectorXd CP = this->Lap_OffDiag * Eigen::VectorXd::Ones(this->count_term_nodes());
	node_evecs.resize(to_print);
	pressure_contributions.resize(to_print);
	for(size_t n = 0; n < to_print; n++)
	{
		evalues[n] = this->modified_laplacian_evalues[indices[n]];
		ceffs[n] = this->modified_laplacian_ceffs[indices[n]];
		Eigen::VectorXd evecth = this->modified_laplacian_evectors[indices[n]];
		node_evecs[n] = Eigen::VectorXd::Zero(this->count_nodes());
		node_evecs[n].segment(1, evecth.size()) = evecth;
		//contribution to pressure solution on internal nodes assuming unit pressure drop
		pressure_contributions[n] = Eigen::VectorXd::Zero(this->count_nodes());
		double weight = evecth.dot(CP) / evalues[n];
		pressure_contributions[n].segment(1, this->modified_laplacian_evectors[indices[n]].size())
				                            = weight * evecth;
	}

	//print
	this->print_evectors_vtk(filename, Nvecs, evalues, ceffs, node_evecs, pressure_contributions, sort_option);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_maury_evectors_vtk(const std::string & filename, 
											   const size_t & Nvecs, const int & sort_option)
{
	size_t to_print = std::min(Nvecs, size_t(this->maury_evalues.size()));
	std::vector<size_t> indices(size_t(this->maury_evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	switch(sort_option)
	{
	case LARGEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->maury_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	case SMALLEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->maury_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) < sort_vec(i2)); });
		} break;

	case DOMINANT_SORT:
		{
			Eigen::VectorXd sort_vec = this->maury_sqrd_Q_ext_contribution;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	default:
		break;
	}

	//convert maury eigenvectors into edge values for visualisation purposes
	std::vector<Eigen::VectorXd> edge_evecs, flux_contributions;
	Eigen::VectorXd evalues = Eigen::VectorXd::Zero(to_print), ceffs = Eigen::VectorXd::Zero(to_print);
	edge_evecs.resize(to_print);
	flux_contributions.resize(to_print);
	for(size_t n = 0; n < to_print; n++)
	{
		evalues[n] = this->maury_evalues[indices[n]];
		this->fill_fluxes_from_term_vals(edge_evecs[n], this->maury_evectors[indices[n]]);
		edge_evecs[n] = edge_evecs[n].array() / this->Nt.array();
		//get flux contirbutions here
		Eigen::VectorXd term_flux_contrib = this->maury_evectors[indices[n]].sum()*this->maury_evectors[indices[n]]
		                                  / (this->maury_evalues[indices[n]]);
		this->fill_fluxes_from_term_vals(flux_contributions[n], term_flux_contrib);
	}
	//print
	this->print_evectors_vtk(filename, Nvecs, evalues, ceffs, edge_evecs, flux_contributions, sort_option);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_adjacency_evectors_vtk(const std::string & filename, 
															 const size_t & Nvecs, const int & sort_option)
{
	size_t to_print = std::min(Nvecs, size_t(this->adjacency_evalues.size()));
	std::vector<size_t> indices(size_t(this->adjacency_evalues.size()));
	iota(indices.begin(), indices.end(), 0);
	switch(sort_option)
	{
	case LARGEST_SORT:
	case DOMINANT_SORT:
		{
			Eigen::VectorXd sort_vec = this->adjacency_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) > sort_vec(i2)); });
		} break;

	case SMALLEST_SORT:
		{
			Eigen::VectorXd sort_vec = this->adjacency_evalues;
			std::sort(indices.begin(), indices.end(), [&sort_vec](int i1, int i2){ return (sort_vec(i1) < sort_vec(i2)); });
		} break;

	default:
		break;
	}

	//copy truncated eigenvectors to full network
	std::vector<Eigen::VectorXd> node_evecs;
	Eigen::VectorXd evalues = Eigen::VectorXd::Zero(to_print);
	Eigen::VectorXd blank = Eigen::VectorXd::Zero(to_print);
	std::vector<Eigen::VectorXd> all_blanks(to_print);
	node_evecs.resize(to_print);
	for(size_t n = 0; n < to_print; n++)
	{
		evalues[n] = this->adjacency_evalues[indices[n]];
		node_evecs[n] = this->adjacency_evectors[indices[n]];
		all_blanks[n] = Eigen::VectorXd::Zero(node_evecs[n].size());
	}
	//print
	this->print_evectors_vtk(filename, Nvecs, evalues, blank, node_evecs, all_blanks, sort_option);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_unit_pressure_drop_solutions_vtk(const std::string & filename)
{
	std::unordered_map<std::string, std::vector<double>> extra_vals;
	extra_vals["Flux"] = std::vector<double>(this->flux.data(), this->flux.data() + this->flux.size());
	Eigen::VectorXd TN;
	this->get_term_count_edges(TN);

	Eigen::VectorXd Flux_PTN = this->flux.array() / TN.array();
	//get number of term nodes in an edge
	extra_vals["Flux_per_TN"] = std::vector<double>(Flux_PTN.data(), Flux_PTN.data() + Flux_PTN.size());
	extra_vals["Pressure"] = std::vector<double>(this->pressure.data(), this->pressure.data() + this->pressure.size());
	if(modified_laplacian_evalues.size() > 0)
	{
		extra_vals["Laplace_flux_approx"] = std::vector<double>(this->modified_laplacian_flux.data(), 
			                              this->modified_laplacian_flux.data() + this->modified_laplacian_flux.size());
		extra_vals["Laplace_pressure_approx"] = std::vector<double>(this->modified_laplacian_pressure.data(),
			                      this->modified_laplacian_pressure.data() + this->modified_laplacian_pressure.size());
	}
	if(maury_evalues.size() > 0)
	{
		extra_vals["Maury_flux_approx"] = std::vector<double>(this->maury_flux.data(), this->maury_flux.data() + this->maury_flux.size());
		extra_vals["Maury_pressure_approx"] = std::vector<double>(this->maury_pressure.data(), this->maury_pressure.data() + this->maury_pressure.size());
	}

	this->print_vtk(filename, 1.0, extra_vals);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_graph_laplacian_modes_summary_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(1);
	headers[0] = std::string("Eigenvalue");
	int nevals = int(this->graph_laplacian_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals);
	for(size_t n = 0; n < nevals; n++)
	{
		data[n].resize(1);
		data[n][0] = this->graph_laplacian_evalues[n];
	}
	write_csv_file<double>(filename, headers, data);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_truncated_laplacian_modes_summary_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(3);
	headers[0] = std::string("Eigenvalue");
	headers[1] = std::string("Absolute_effective_resistance");
	headers[2] = std::string("Relative_effective_resistance");
	int nevals = int(this->truncated_laplacian_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals);
	for(size_t n = 0; n < nevals; n++)
	{
		data[n].resize(3);
		data[n][0] = this->truncated_laplacian_evalues[n];
		data[n][1] = this->truncated_laplacian_reffs[n];
		data[n][2] = this->truncated_laplacian_reffs[n] * this->Ceff;
	}
	write_csv_file<double>(filename, headers, data);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_modified_laplacian_modes_summary_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(5);
	headers[0] = std::string("Eigenvalue");
	headers[1] = std::string("Absolute_effective_conductance");
	headers[2] = std::string("Relative_effective_conductance");
	headers[3] = std::string("Relative_effective_conductance_alt");
	headers[4] = std::string("Relative_sqrd_P_int_frac");
	int nevals = int(this->modified_laplacian_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals+1);
	//add starting point (0 modes) to account for Pext contribution
	data[0].resize(5);
	data[0][0] = -1;
	data[0][1] = this->Lap_Ext.sum();
	data[0][2] = data[0][1] / this->Ceff;
	data[0][3] = this->Lap_11 / this->Ceff;
	data[0][4] = 0;  //does not contribute to Pint
	for(size_t n = 0; n < nevals; n++)
	{
		data[n+1].resize(5);
		data[n+1][0] = this->modified_laplacian_evalues[n];
		data[n+1][1] = this->modified_laplacian_ceffs[n];
		data[n+1][2] = this->modified_laplacian_ceffs[n] / this->Ceff;
		data[n+1][3] = this->modified_laplacian_ceffs_alt[n] / this->Ceff;
		data[n+1][4] = this->modified_laplacian_sqrd_P_int_contribution[n] / this->P_int_sqrd;
	}


	write_csv_file<double>(filename, headers, data);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_maury_modes_summary_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(4);
	headers[0] = std::string("Eigenvalue");
	headers[1] = std::string("Absolute_effective_conductance");
	headers[2] = std::string("Relative_effective_conductance");
	headers[3] = std::string("Relative_sqrd_Q_ext_frac");
	int nevals = int(this->maury_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals);
	for(size_t n = 0; n < nevals; n++)
	{
		data[n].resize(4);
		data[n][0] = this->maury_evalues[n];
		data[n][1] = this->maury_ceffs[n];
		data[n][2] = this->maury_ceffs[n] / this->Ceff;
		data[n][3] = this->maury_sqrd_Q_ext_contribution[n]/this->Q_ext_sqrd;
	}

	write_csv_file<double>(filename, headers, data);
}

template<int RTYPE> void SpectralNetwork<RTYPE>::print_adjacency_modes_summary_csv(const std::string & filename)
{
	std::vector<std::string> headers;
	headers.resize(1);
	headers[0] = std::string("Eigenvalue");
	int nevals = int(this->adjacency_evalues.size());
	std::vector<std::vector<double>> data;
	data.resize(nevals);
	for(size_t n = 0; n < nevals; n++)
	{
		data[n].resize(1);
		data[n][0] = this->adjacency_evalues[n];
	}
	write_csv_file<double>(filename, headers, data);
}

#endif
