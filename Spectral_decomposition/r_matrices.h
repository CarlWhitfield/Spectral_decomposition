#ifndef R_MATRICES_H
#define R_MATRICES_H

#include<network_3D.h>
#include<Eigen/Dense>
#include<Eigen/Sparse>

//choices for Rmat OPT
#define RMAT1 0  //using T matrix
#define RMAT2 1  //using tree directly

//base class for Rmatrix
template<int OPT> class Rbase
{
protected:
	typedef network::Network<network::Node,network::Edge<network::Node>> ANetwork;
	ANetwork * tree;
	size_t Ntermnodes, Nedges;
	Eigen::VectorXd rvec;
	double shift;
public:
	Rbase(){};
	Rbase(ANetwork *t, Eigen::VectorXd conductances)
	{
		this->tree = t;
		this->Ntermnodes = t->count_term_nodes();
		this->Nedges = t->count_edges();
		this->rvec = Eigen::VectorXd::Ones(conductances.size()).array() / conductances.array();
		this->shift = 0;
	}
	virtual void set_shift(const double & s){ this->shift = s; }
	inline int rows() { return int(this->Ntermnodes); }
    inline int cols() { return int(this->Ntermnodes); }

	inline ANetwork * return_tree_pointer() const { return this->tree; }
	inline const Eigen::VectorXd& return_rvec() const { return this->rvec; }

	virtual void perform_op(double *x_in, double *y_out){};
};

//direct R matrix object
template<int OPT> class Rmat: public Rbase<OPT>  //class to perform multiplication by R matrix using T decomposition
{
protected:
	Eigen::SparseMatrix<double> T, T_transpose;
	void initialise(ANetwork *tree,
		            const Eigen::VectorXd & conductance);
public:
	Rmat():Rbase<OPT>(){};
	Rmat(ANetwork * t, const Eigen::VectorXd & conductance)
		:Rbase<OPT>(t, conductance)
	{
		this->initialise(t, conductance);
	}

	void perform_op(double *x_in, double *y_out);
	inline void print_Tmat()
	{
		std::cout << "T: " << T.toDense() << '\n';
	}
	inline void print_diag_r()
	{
		std::cout << "diag(r): " << this->diag_r.toDense() << '\n';
	}

};

//inverse of R matrix object
template<int OPT> class Rmatinv: public Rbase<OPT>
{
private:
	Eigen::SparseMatrix<double> full_mat;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	Eigen::VectorXd rvec_raw;
	int full_mat_dim;
	double shift;
	void reinitialise();
public:
	Rmatinv():Rbase<OPT>(){};
	Rmatinv(ANetwork *tree, 
	        const Eigen::VectorXd & conductance):Rbase<OPT>(tree, conductance)
	{
		this->rvec_raw = this->rvec;
		this->reinitialise();
	}
	
	void set_shift(const double & shift);
	inline int rows() { return int(this->Ntermnodes); }
    inline int cols() { return int(this->Ntermnodes); }
	void perform_op(double *x_in, double *y_out);
};

//functions for Rmatinv
template<int OPT> void Rmatinv<OPT>::reinitialise()
{
	//do not store Rmat, instead make pressure drop solver
	this->full_mat_dim = int(this->tree->count_nodes() + this->tree->count_edges());    //pressure inside + term_flux + flux in
	this->full_mat = Eigen::SparseMatrix<double>(this->full_mat_dim, this->full_mat_dim);
	std::vector<Eigen::Triplet<double>> A_fill;
	A_fill.reserve(3*this->full_mat_dim);
	//solution will consist of (q, P)
	//rhs is (0, -P_in)
	//loop over edges rij qij + Pj - Pi = 0
	for(size_t j = 0; j < this->tree->count_edges(); j++)
	{
		size_t ki = this->tree->get_node_in_index(j);
		size_t ko = this->tree->get_node_out_index(j);
		A_fill.push_back(Eigen::Triplet<double>(int(j),int(j),this->rvec[j]));
		A_fill.push_back(Eigen::Triplet<double>(int(j),int(this->tree->count_edges() + ko), 1.0));
		A_fill.push_back(Eigen::Triplet<double>(int(j),int(this->tree->count_edges() + ki), -1.0));
	}
	//first node set  P_0 = 0
	A_fill.push_back(Eigen::Triplet<double>(int(this->tree->count_edges()), int(this->tree->count_edges()), 1.0));
	//interior nodes (flux conservation) qik - qkm
	for(size_t k = 1; k < this->tree->get_first_term_index(); k++)
	{
		for(size_t ji = 0; ji < this->tree->count_edges_in(k); ji++)
		{
			size_t j = this->tree->get_edge_in_index(k,ji);
			A_fill.push_back(Eigen::Triplet<double>(int(this->tree->count_edges() + k), int(j), 1.0));
		}
		for(size_t jo = 0; jo < this->tree->count_edges_out(k); jo++)
		{
			size_t j = this->tree->get_edge_out_index(k,jo);
			A_fill.push_back(Eigen::Triplet<double>(int(this->tree->count_edges() + k), int(j), -1.0));
		}
	}
	//terminal nodes pressure = -vector in
	for(size_t k = tree->get_first_term_index(); k < tree->count_nodes(); k++)
	{
		A_fill.push_back(Eigen::Triplet<double>(int(this->tree->count_edges() + k), int(this->tree->count_edges() + k), -1.0));
	}
	this->full_mat.setFromTriplets(A_fill.begin(), A_fill.end());
	this->solver.compute(this->full_mat);
}

template<int OPT> void Rmatinv<OPT>::set_shift(const double & s)
{
	this->shift = s;
	for(size_t k = this->tree->get_first_term_index(); k <this->tree->count_nodes(); k++)
	{
		size_t j = this->tree->get_edge_in_index(k,0);
		this->rvec[j] = this->rvec_raw[j] - s;
	}
	this->reinitialise();
}

template<int OPT> void Rmatinv<OPT>::perform_op(double *x_in, double *y_out)
{
	Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(x_in, this->Ntermnodes);
	Eigen::VectorXd b(this->full_mat.rows());
	b << Eigen::VectorXd::Zero(this->full_mat.rows() - this->Ntermnodes), x;
	Eigen::VectorXd y = this->solver.solve(b);
	for(size_t kt = this->tree->get_first_term_index(); kt < this->tree->count_nodes(); kt++)
	{
		size_t j = this->tree->get_edge_in_index(kt,0);
		y_out[kt - this->tree->get_first_term_index()] = y[j];
	}
}

#endif