#include<network_3D.h>

//choices for Rmat OPT
#define RMAT1 0  //using T matrix
#define RMAT2 1  //using tree directly

template<int OPT> class Rbase
{
protected:
	network::Network<network::Node,network::Edge<network::Node>> *tree;
	size_t Ntermnodes, Nedges;
	Eigen::VectorXd rvec;
	double shift;
public:
	Rbase(){};
	Rbase(network::Network<network::Node,network::Edge<network::Node>> *t, Eigen::VectorXd conductances)
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

	inline network::Network<network::Node,network::Edge<network::Node>>* return_tree_pointer() const { return this->tree; }
	inline const Eigen::VectorXd& return_rvec() const { return this->rvec; }

	virtual void perform_op(double *x_in, double *y_out){};
};

template<int OPT> class Rmat: public Rbase<OPT>  //class to perform multiplication by R matrix using T decomposition
{
private:
	Eigen::SparseMatrix<double> T, T_transpose;
public:
	Rmat():Rbase<OPT>(){};
	Rmat(network::Network<network::Node,network::Edge<network::Node>> *tree, const Eigen::VectorXd & conductance):Rbase<RMAT2>(tree, conductance_weight){};
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

template<> Rmat<RMAT1>::Rmat(network::Network<network::Node,network::Edge<network::Node>> *tree, const Eigen::VectorXd & conductance_weight):Rbase<RMAT1>(tree, conductance_weight)
{
	//Eigen uses column ordering, so calc T transpose first
	this->T_transpose = Eigen::SparseMatrix<double>(this->Ntermnodes, this->Nedges);

	std::vector<Eigen::Triplet<double>> TT_fill;
	TT_fill.reserve(tree->count_horsfield_orders() * this->Ntermnodes);

	//loop over horsfield order 0 (term nodes)
	std::vector<std::vector<int>> ktvec;
	ktvec.resize(tree->count_edges()); 
	for(size_t k = tree->get_first_term_index(); k < tree->count_nodes(); k++)
	{
		size_t k_term = k - tree->get_first_term_index();
		size_t j = tree->get_edge_in_index(k,0);
		TT_fill.push_back(Eigen::Triplet<double>(int(k_term), int(j), 1.0));
		ktvec[j].push_back(int(k_term));  //termnial nodes identified
	}

	for(size_t h = 1; h < tree->count_horsfield_orders(); h++)
	{
		std::vector<Eigen::VectorXd> ktvec_next;
		ktvec_next.resize(tree->count_edges_in_horsfield_order(h));
		//loop over edges in horsfield order
		for(size_t i_edge = 0; i_edge < tree->count_edges_in_horsfield_order(h); i_edge++)
		{
			size_t j = tree->get_edge_index_from_horsfield_order(h,i_edge);
			size_t k_out = tree->get_node_out_index(j);
			//loop over edges out
			for(size_t i_out = 0; i_out < tree->count_edges_out(k_out); i_out++)
			{
				//get index of edges descended from j
				size_t j_out = tree->get_edge_out_index(k_out, i_out);
				//loop over column j_out
				if(i_out==0)
				{
					ktvec[j] = ktvec[j_out];
				}
				else
				{
					ktvec[j].insert(ktvec[j].end(), ktvec[j_out].begin(), ktvec[j_out].end());
				}
			}
			for(size_t i_add=0; i_add < ktvec[j].size(); i_add++)  //cannot insert while using matrix iterator
			{
				//add terms to column j
				size_t k = ktvec[j][i_add];
				TT_fill.push_back(Eigen::Triplet<double>(int(k), int(j), 1.0));
			}
		}
	}
	this->T_transpose.setFromTriplets(TT_fill.begin(), TT_fill.end());
	this->T_transpose.makeCompressed();
	this->T = T_transpose.transpose();
}

template<> void Rmat<RMAT1>::perform_op(double *x_in, double *y_out)
{
	Eigen::VectorXd v = Eigen::Map<Eigen::VectorXd>(x_in, this->Ntermnodes);
	//std::cout << v << '\n';
	//create v mapped to edges
	Eigen::VectorXd ve = this->T * v;

	//std::cout << "ve = " << ve << '\n';
	//generate x vals for edges
	Eigen::VectorXd x = this->rvec.array() * ve.array();

	//std::cout << "x = " << x << '\n';
	//calc u and return
	Eigen::VectorXd u = this->T_transpose * x;
	for(size_t k = 0; k < this->Ntermnodes; k++)
	{
		y_out[k] = u[k] - this->shift*v[k];
	}
}

template<> void Rmat<RMAT2>::perform_op(double *x_in, double *y_out)
{
	Eigen::VectorXd v = Eigen::Map<Eigen::VectorXd>(x_in, this->Ntermnodes);

	//ve = v mapped to edges
	Eigen::VectorXd ve = Eigen::VectorXd::Zero(this->Nedges);

	//loop over horsfield order 0 (term nodes)
	for(size_t k = this->tree->get_first_term_index(); k < this->tree->count_nodes(); k++)
	{
		size_t k_term = k - this->tree->get_first_term_index();
		size_t j = this->tree->get_edge_in_index(k,0);
		ve[j] = v[k_term];
	}

	//loop over horsfield orders > 0
	for(size_t h = 1; h < this->tree->count_horsfield_orders(); h++)
	{
		//loop over edges in horsfield order
		for(size_t i_edge = 0; i_edge < this->tree->count_edges_in_horsfield_order(h);
			i_edge++)
		{
			size_t j = this->tree->get_edge_index_from_horsfield_order(h,i_edge);
			size_t k_out = this->tree->get_node_out_index(j);
			for(size_t i_out = 0; i_out < this->tree->count_edges_out(k_out);
				i_out++)
			{
				//get index of edges descended from j
				size_t j_out = this->tree->get_edge_out_index(k_out, i_out);
				ve[j] += ve[j_out];   //sum over descended edges
			}
		}
	}
	//std::cout << "ve = " << ve << '\n';
	//generate x vals for edges
	Eigen::VectorXd x = this->rvec.array() * ve.array();
	//std::cout << "x = " << x << '\n';
	//fill vector of u vals on edges
	Eigen::VectorXd ue = Eigen::VectorXd::Zero(this->Nedges);
	//weibel order 0
	ue[0] = x[0];
	//loop over weibel orders > 0
	for(size_t w = 1; w < this->tree->count_weibel_orders(); w++)
	{
		for(size_t i_edge = 0; i_edge < this->tree->count_edges_in_weibel_order(w);
			i_edge++)
		{
			size_t j = this->tree->get_edge_index_from_weibel_order(w,i_edge);
			size_t j_in = this->tree->get_edge_in_index(tree->get_node_in_index(j),0);
			ue[j] = ue[j_in] + x[j];
		}
	}

	//map edge u vals back onto term nodes
	Eigen::VectorXd u = Eigen::VectorXd::Zero(this->Ntermnodes);
	for(size_t k = this->tree->get_first_term_index(); k < this->tree->count_nodes(); k++)
	{
		size_t k_term = k - this->tree->get_first_term_index();
		size_t j = this->tree->get_edge_in_index(k,0);
		y_out[k_term] = ue[j] - this->shift*v[k_term];
	}
}

template<int OPT> class Rmatinv: public Rbase<OPT>
{
private:
	Eigen::SparseMatrix<double> full_mat;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	Eigen::VectorXd rvec_raw;
	double full_mat_dim, shift;
	void reinitialise();
public:
	Rmatinv():Rbase<OPT>(){};
	Rmatinv(network::Network<network::Node,network::Edge<network::Node>> *tree, 
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

template<int OPT> void Rmatinv<OPT>::reinitialise()
{
	//do not store Rmat, instead make pressure drop solver
	this->full_mat_dim = this->tree->count_nodes() + this->tree->count_edges();    //pressure inside + term_flux + flux in
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
			A_fill.push_back(Eigen::Triplet<double>(int(this->tree->count_edges() + k), j, 1.0));
		}
		for(size_t jo = 0; jo < this->tree->count_edges_out(k); jo++)
		{
			size_t j = this->tree->get_edge_out_index(k,jo);
			A_fill.push_back(Eigen::Triplet<double>(int(this->tree->count_edges() + k), j, -1.0));
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

//this doesnt work for some reason