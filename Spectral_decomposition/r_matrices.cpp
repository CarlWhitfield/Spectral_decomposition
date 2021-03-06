#include"r_matrices.h"

//Rmat functions
template<> void Rmat<RMAT1>::initialise(ANetwork * tree, 
	                                    const Eigen::VectorXd & conductance_weight)          
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