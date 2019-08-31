#include<network_3D.h>

//choices for Rmat OPT
#define RMAT1 0  //using T matrix
#define RMAT2 1  //using tree directly

template<int OPT> class Rmat  //class to perform multiplication by R matrix using T decomposition
{
private:
	network::Network<network::Node,network::Edge<network::Node>> *tree;
	Eigen::SparseMatrix<double> T, T_transpose;
	Eigen::SparseMatrix<double> diag_r;
	size_t Ntermnodes, Nedges;
	double diag_sup;
	bool flipped;
public:
	Rmat(){};
	Rmat(network::Network<network::Node,network::Edge<network::Node>> *tree, const Eigen::VectorXd & conductance);
	Rmat(const Rmat<OPT> & r, const double & diag, const bool & flip = false)
	{
		*(this) = r;
		this->diag_sup = diag;
		this->flipped = flip;
	}
	inline int rows() { return int(this->Ntermnodes); }
    inline int cols() { return int(this->Ntermnodes); }
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

template<> Rmat<RMAT1>::Rmat(network::Network<network::Node,network::Edge<network::Node>> *tree,  const Eigen::VectorXd & conductance_weight)
{
	//initialise given some tree
	this->Nedges = tree->count_edges();
	this->Ntermnodes = tree->count_term_nodes();
	//Eigen uses column ordering, so calc T transpose first
	this->T_transpose = Eigen::SparseMatrix<double>(this->Ntermnodes, this->Nedges);
	this->diag_r = Eigen::SparseMatrix<double>(this->Nedges, this->Nedges);
	std::vector<Eigen::Triplet<double>> TT_fill, diag_r_fill;
	TT_fill.reserve(tree->count_horsfield_orders() * this->Ntermnodes);
	diag_r_fill.reserve(this->Nedges);

	//loop over horsfield order 0 (term nodes)
	std::vector<std::vector<int>> ktvec;
	ktvec.resize(tree->count_edges()); 
	for(size_t k = tree->get_first_term_index(); k < tree->count_nodes(); k++)
	{
		size_t k_term = k - tree->get_first_term_index();
		size_t j = tree->get_edge_in_index(k,0);
		TT_fill.push_back(Eigen::Triplet<double>(int(k_term), int(j), 1.0));
		diag_r_fill.push_back(Eigen::Triplet<double>(int(j), int(j), 1.0/conductance_weight[j]));
		ktvec[j].push_back(k_term);  //termnial nodes identified
	}

	for(size_t h = 1; h < tree->count_horsfield_orders(); h++)
	{
		std::vector<Eigen::VectorXd> ktvec_next;
		ktvec_next.resize(tree->count_edges_in_horsfield_order(h));
		//loop over edges in horsfield order
		for(size_t i_edge = 0; i_edge < tree->count_edges_in_horsfield_order(h); i_edge++)
		{
			size_t j = tree->get_edge_index_from_horsfield_order(h,i_edge);
			diag_r_fill.push_back(Eigen::Triplet<double>(int(j), int(j), 1.0/conductance_weight[j]));
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
	this->diag_r.setFromTriplets(diag_r_fill.begin(), diag_r_fill.end());
	this->T_transpose.makeCompressed();
	this->T = T_transpose.transpose();
	this->diag_sup = 0;
	this->flipped = false;
}

template<> Rmat<RMAT2>::Rmat(network::Network<network::Node,network::Edge<network::Node>> *tree_in, const Eigen::VectorXd & conductance_weight)
{
	this->tree = tree_in;
	this->Nedges = tree_in->count_edges();
	this->Ntermnodes = tree_in->count_term_nodes();
	this->diag_r = Eigen::SparseMatrix<double>(this->Nedges, this->Nedges);
	this->diag_r.reserve(this->Nedges);
	for(size_t j = 0; j < tree_in->count_edges(); j++)
	{
		this->diag_r.insert(j,j) = 1.0/conductance_weight[j];
	}
	this->diag_sup = 0;
	this->flipped = false;
}



template<> void Rmat<RMAT1>::perform_op(double *x_in, double *y_out)
{
	Eigen::VectorXd v = Eigen::Map<Eigen::VectorXd>(x_in, this->Ntermnodes);
	//std::cout << v << '\n';
	//create v mapped to edges
	Eigen::VectorXd ve = this->T * v;

	//std::cout << "ve = " << ve << '\n';
	//generate x vals for edges
	Eigen::VectorXd x = this->diag_r * ve;

	//std::cout << "x = " << x << '\n';
	//calc u and return
	Eigen::VectorXd u = this->T_transpose * x;
	for(size_t k = 0; k < this->Ntermnodes; k++)
	{
		if(this->flipped) y_out[k] = -u[k] + diag_sup*v[k];
		else y_out[k] = u[k] + diag_sup*v[k];
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
	Eigen::VectorXd x = this->diag_r * ve;
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
		if(this->flipped) y_out[k_term] = -ue[j]  + diag_sup*v[k_term];
		else y_out[k_term] = ue[j]  + diag_sup*v[k_term];
	}
}