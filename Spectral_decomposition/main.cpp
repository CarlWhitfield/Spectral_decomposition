#include<list_template.h>
#include"spectral_operations.h"
#include"define_spectral_options.h"
#include"define_spectral_params.h"

void read_options(const int &argc, char** argv, SpectralOptionList* o, SpectralParameterList *p)  //constructs lung simulation from command line args
{
	//5 possible arguments in total, all filenames, sort by extensions
	std::vector<std::string> extensions;
	extensions.resize(5);
	extensions[0] = OPTIONS_FILE_EXT;
	extensions[1] = PARAMS_FILE_EXT;
	extensions[2] = NODE_FILE_EXT;
	extensions[3] = TERM_NODE_FILE_EXT;
	extensions[4] = BRANCH_FILE_EXT;

	o->get_filenames_from_args(extensions, argc, argv);
	if( o->filename_exists(OPTIONS_FILE_EXT) ) o->read_file(o->get_filename(OPTIONS_FILE_EXT));  //if options file given, read it in
	if( o->filename_exists(PARAMS_FILE_EXT) ) p->read_file(o->get_filename(PARAMS_FILE_EXT));  //if options file given, read it in
	p->check_and_convert(o);  //run checks and conversions here
}

std::shared_ptr<SpectralNetwork<>> setup_tree(SpectralOptionList* o, SpectralParameterList* p)
{
	int tree_type;
	std::shared_ptr<SpectralNetwork<>> tree;
	switch(o->get_option<char>(EDGE_WEIGHT_KEY)->get_value())
	{
	case RESISTANCE_WEIGHTING:
		{
			tree_type = RESISTANCE_NETWORK;
		} break;

	case DIFFUSION_NETWORK:
		{
			tree_type = DIFFUSION_NETWORK;
		} break;
	default:
		{
			std::cerr << "Edge weight option not recognised, assuming resistance weighting.\n";
			tree_type = RESISTANCE_NETWORK;
		}
	}
	

	switch(o->get_option<char>(TREE_KEY)->get_value())  //type of tree
	{
	case MODEL_FROM_FILE:
		{
			tree = std::make_shared<SpectralNetwork<>>(o->get_filename(NODE_FILE_EXT), o->get_filename(BRANCH_FILE_EXT),
				                                       o->get_filename(TERM_NODE_FILE_EXT), tree_type);
		} break;

	case ASYMM_MODEL:
		{
			tree = std::shared_ptr<SpectralNetwork<>>(new SpectralNetwork<>(
				                                       size_t(p->get_param<int>(N_GENERATIONS_KEY)->get_value()),
						                               p->get_param<double>(SEED_RADIUS_KEY)->get_value(), 
						                               p->get_param<double>(SEED_LENGTH_KEY)->get_value(),
						                               p->get_param<double>(SCALE_FACTOR_KEY)->get_value(), 
						                               p->get_param<double>(ASYMMETRY_KEY)->get_value(),
						                               tree_type));
		} break;

	case PERT_MODEL:
		{
			tree = std::shared_ptr<SpectralNetwork<>>(new SpectralNetwork<>(
				                                       size_t(p->get_param<int>(N_GENERATIONS_KEY)->get_value()),
						                               p->get_param<double>(SEED_RADIUS_KEY)->get_value(), 
						                               p->get_param<double>(SEED_LENGTH_KEY)->get_value(),
						                               p->get_param<double>(SCALE_FACTOR_KEY)->get_value(), 
						                               0.0, p->get_param<double>(PERTURBATION_FRACTION_KEY)->get_value(),
						                               size_t(p->get_param<int>(PERTURBATION_GEN_KEY)->get_value()),
						                               tree_type));
		} break;

	default:
		{
			std::cerr << "Model type not recognised. Aborting. \n";
			abort_on_failure();
			return tree;
		} break;
	}
	return tree;
}

int compute_and_print_graph_laplacian_spectrum(SpectralNetwork<> *tree, SpectralOptionList *o, SpectralParameterList *p)
{
	std::cout << "Computing Laplacian Spectrum...\n";
	switch(o->get_option<char>(CUTOFF_MODE_KEY)->get_value())
	{
	case NO_CUTOFF:
		{
			tree->compute_full_graph_laplacian_spectrum();
		} break;

	case FRAC_CUTOFF:
		{
			double frac = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(frac == 1)   //equivalent to computing full spectrum
			{
				tree->compute_full_graph_laplacian_spectrum();
			}
			else
			{
				if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
				{
					int Nlarge = std::max(int(1),int(frac*tree->count_nodes()));
					tree->compute_partial_graph_laplacian_spectrum(0,Nlarge);
				}
				else
				{
					int Nsmall = std::max(int(1),int(frac*tree->count_nodes()));
					tree->compute_partial_graph_laplacian_spectrum(Nsmall,0);
				}
			}
		} break;

	case VALUE_CUTOFF:
		{
			double scaled_value = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
			{
				tree->compute_all_graph_laplacian_modes_above_cutoff(scaled_value);
			}
			else
			{
				tree->compute_all_graph_laplacian_modes_below_cutoff(scaled_value);
			}
		} break;

	default:
		{
			std::cerr << "Do not recognise cutoff option. Assuming full spectrum.\n";
			tree->compute_full_graph_laplacian_spectrum();
		}
	}

	std::cout << "Outputting...\n";
	std::string filehead = o->generate_output_name();
	std::stringstream ss;
	ss << filehead << "_summary.csv";
	tree->print_graph_laplacian_modes_summary_csv(ss.str());
	if(o->get_option<bool>(PRINT_VTKS)->get_value() || o->get_option<bool>(PRINT_CSVS)->get_value())  //print out modes in vtk file
	{
		//continue here
		int Nprint = p->get_param<int>(N_PRINT_KEY)->get_value();
		ss.clear();
		ss.str("");
		switch(o->get_option<char>(OUTPUT_SORT_KEY)->get_value())
		{
		case PRINT_DOMINANT_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_graph_laplacian_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_graph_laplacian_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
			} break;

		case PRINT_LARGEST_MODES:
			{
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_graph_laplacian_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_graph_laplacian_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
			} break;

		case PRINT_SMALLEST_MODES:
			{
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_graph_laplacian_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_graph_laplacian_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;

		case PRINT_ALL_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_graph_laplacian_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_graph_laplacian_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_graph_laplacian_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_graph_laplacian_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
				
				ss.clear();
				ss.str("");
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_graph_laplacian_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_graph_laplacian_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;
		}
	}
	return 0;
}

int compute_and_print_truncated_laplacian_spectrum(SpectralNetwork<> *tree, SpectralOptionList *o, SpectralParameterList *p)
{
	std::cout << "Computing Truncated  Laplacian Spectrum...\n";
	switch(o->get_option<char>(CUTOFF_MODE_KEY)->get_value())
	{
	case NO_CUTOFF:
		{
			tree->compute_full_truncated_laplacian_spectrum();
		} break;

	case FRAC_CUTOFF:
		{
			double frac = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(frac == 1)   //equivalent to computing full spectrum
			{
				tree->compute_full_truncated_laplacian_spectrum();
			}
			else
			{
				if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
				{
					int Nlarge = std::max(int(1),int(frac*(tree->count_nodes() - tree->count_term_nodes() + 1)));
					tree->compute_partial_truncated_laplacian_spectrum(0,Nlarge);
				}
				else
				{
					int Nsmall = std::max(int(1),int(frac*(tree->count_nodes() - tree->count_term_nodes() + 1)));
					tree->compute_partial_truncated_laplacian_spectrum(Nsmall,0);
				}
			}
		} break;

	case VALUE_CUTOFF:
		{
			double scaled_value = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
			{
				tree->compute_all_truncated_laplacian_modes_above_cutoff(scaled_value);
			}
			else
			{
				tree->compute_all_truncated_laplacian_modes_below_cutoff(scaled_value);
			}
		} break;

	default:
		{
			std::cerr << "Do not recognise cutoff option. Assuming full spectrum.\n";
			tree->compute_full_truncated_laplacian_spectrum();
		}
	}

	std::cout << "Outputting...\n";
	std::string filehead = o->generate_output_name();
	std::stringstream ss;
	ss << filehead << "_summary.csv";
	tree->print_truncated_laplacian_modes_summary_csv(ss.str());
	if(o->get_option<bool>(PRINT_VTKS)->get_value() || o->get_option<bool>(PRINT_CSVS)->get_value())  //print out modes in vtk file
	{
		//continue here
		int Nprint = p->get_param<int>(N_PRINT_KEY)->get_value();
		ss.clear();
		ss.str("");
		switch(o->get_option<char>(OUTPUT_SORT_KEY)->get_value())
		{
		case PRINT_DOMINANT_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_truncated_laplacian_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_truncated_laplacian_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
			} break;

		case PRINT_LARGEST_MODES:
			{
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_truncated_laplacian_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_truncated_laplacian_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}

			} break;

		case PRINT_SMALLEST_MODES:
			{
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_truncated_laplacian_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_truncated_laplacian_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;

		case PRINT_ALL_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_truncated_laplacian_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_truncated_laplacian_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
				
				ss.clear();
				ss.str("");
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_truncated_laplacian_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_truncated_laplacian_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
				
				ss.clear();
				ss.str("");
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_truncated_laplacian_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_truncated_laplacian_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;
		}
	}
	return 0;
}

int compute_and_print_modified_laplacian_spectrum(SpectralNetwork<> *tree, SpectralOptionList *o, SpectralParameterList *p)
{
	std::cout << "Computing Modified Laplacian Spectrum...\n";
	switch(o->get_option<char>(CUTOFF_MODE_KEY)->get_value())
	{
	case NO_CUTOFF:
		{
			tree->compute_full_modified_laplacian_spectrum();
		} break;

	case FRAC_CUTOFF:
		{
			double frac = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(frac == 1)   //equivalent to computing full spectrum
			{
				tree->compute_full_modified_laplacian_spectrum();
			}
			else
			{
				if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
				{
					int Nlarge = std::max(int(1),int(frac*(tree->count_nodes() - tree->count_term_nodes() - 1)));
					tree->compute_partial_modified_laplacian_spectrum(0,Nlarge);
				}
				else
				{
					int Nsmall = std::max(int(1),int(frac*(tree->count_nodes() - tree->count_term_nodes() - 1)));
					tree->compute_partial_modified_laplacian_spectrum(Nsmall,0);
				}
			}
		} break;

	case VALUE_CUTOFF:
		{
			double scaled_value = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
			{
				tree->compute_all_modified_laplacian_modes_above_cutoff(scaled_value);
			}
			else
			{
				tree->compute_all_modified_laplacian_modes_below_cutoff(scaled_value);
			}
		} break;

	default:
		{
			std::cerr << "Do not recognise cutoff option. Assuming full spectrum.\n";
			tree->compute_full_modified_laplacian_spectrum();
		}
	}

	std::cout << "Outputting...\n";
	std::string filehead = o->generate_output_name();
	std::stringstream ss;
	ss << filehead << "_summary.csv";
	tree->print_modified_laplacian_modes_summary_csv(ss.str());
	if(o->get_option<bool>(PRINT_VTKS)->get_value() || o->get_option<bool>(PRINT_CSVS)->get_value())  //print out modes in vtk file
	{
		//continue here
		int Nprint = p->get_param<int>(N_PRINT_KEY)->get_value();
		ss.clear();
		ss.str("");
		switch(o->get_option<char>(OUTPUT_SORT_KEY)->get_value())
		{
		case PRINT_DOMINANT_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{;
					tree->print_modified_laplacian_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_modified_laplacian_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
			} break;

		case PRINT_LARGEST_MODES:
			{
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_modified_laplacian_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_modified_laplacian_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
			} break;

		case PRINT_SMALLEST_MODES:
			{
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_modified_laplacian_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_modified_laplacian_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;

		case PRINT_ALL_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_modified_laplacian_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_modified_laplacian_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_modified_laplacian_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_modified_laplacian_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_modified_laplacian_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_modified_laplacian_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;
		}
	}
	return 0;
}

int compute_and_print_maury_spectrum(SpectralNetwork<> *tree, SpectralOptionList *o, SpectralParameterList *p)
{
	std::cout << "Computing Maury Spectrum...\n";
	switch(o->get_option<char>(CUTOFF_MODE_KEY)->get_value())
	{
	case NO_CUTOFF:
		{
			tree->compute_full_maury_spectrum();
		} break;

	case FRAC_CUTOFF:
		{
			double frac = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(frac == 1)   //equivalent to computing full spectrum
			{
				tree->compute_full_maury_spectrum();
			}
			else
			{
				if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
				{
					int Nlarge = std::max(int(1),int(frac*tree->count_term_nodes()));
					tree->compute_partial_maury_spectrum(0,Nlarge);
				}
				else
				{
					int Nsmall = std::max(int(1),int(frac*tree->count_term_nodes()));
					tree->compute_partial_maury_spectrum(Nsmall,0);
				}
			}
		} break;

	case VALUE_CUTOFF:
		{
			double scaled_value = p->get_param<double>(CUTOFF_POINT_KEY)->get_value()
				                * tree->count_term_nodes() / tree->get_maury_scale_factor();    //for maury matrix, value cutoff is scaled by number of term nodes
			if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
			{
				tree->compute_all_maury_modes_above_cutoff(scaled_value);
			}
			else
			{
				tree->compute_all_maury_modes_below_cutoff(scaled_value);
			}
		} break;

	default:
		{
			std::cerr << "Do not recognise cutoff option. Assuming full spectrum.\n";
			tree->compute_full_maury_spectrum();
		}
	}

	std::cout << "Outputting...\n";
	std::string filehead = o->generate_output_name();
	std::stringstream ss;
	ss << filehead << "_summary.csv";
	tree->print_maury_modes_summary_csv(ss.str());
	if(o->get_option<bool>(PRINT_VTKS)->get_value() || o->get_option<bool>(PRINT_CSVS)->get_value())  //print out modes in vtk file
	{
		int Nprint = p->get_param<int>(N_PRINT_KEY)->get_value();
		ss.clear();
		ss.str("");
		switch(o->get_option<char>(OUTPUT_SORT_KEY)->get_value())
		{
		case PRINT_DOMINANT_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_maury_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_maury_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
			} break;

		case PRINT_LARGEST_MODES:
			{
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_maury_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_maury_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
			} break;

		case PRINT_SMALLEST_MODES:
			{
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_maury_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_maury_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;

		case PRINT_ALL_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_maury_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_maury_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_maury_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_maury_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_maury_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_maury_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;
		}
	}
	return 0;
}

int compute_and_print_adjacency_spectrum(SpectralNetwork<> *tree, SpectralOptionList *o, SpectralParameterList *p)
{
	std::cout << "Computing Adjacency Spectrum...\n";
	switch(o->get_option<char>(CUTOFF_MODE_KEY)->get_value())
	{
	case NO_CUTOFF:
		{
			tree->compute_full_adjacency_spectrum();
		} break;

	case FRAC_CUTOFF:
		{
			double frac = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(frac == 1)   //equivalent to computing full spectrum
			{
				tree->compute_full_adjacency_spectrum();
			}
			else
			{
				if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
				{
					int Nlarge = std::max(int(1),int(frac*(tree->count_nodes())));
					tree->compute_partial_adjacency_spectrum(0,Nlarge);
				}
				else
				{
					int Nsmall = std::max(int(1),int(frac*(tree->count_nodes())));
					tree->compute_partial_adjacency_spectrum(Nsmall,0);
				}
			}
		} break;

	case VALUE_CUTOFF:
		{
			double scaled_value = p->get_param<double>(CUTOFF_POINT_KEY)->get_value();
			if(o->get_option<char>(ORDER_KEY)->get_value() == LARGEST_MODES)
			{
				tree->compute_all_adjacency_modes_above_cutoff(scaled_value);
			}
			else
			{
				tree->compute_all_adjacency_modes_below_cutoff(scaled_value);
			}
		} break;

	default:
		{
			std::cerr << "Do not recognise cutoff option. Assuming full spectrum.\n";
			tree->compute_full_adjacency_spectrum();
		}
	}

	std::cout << "Outputting...\n";
	std::string filehead = o->generate_output_name();
	std::stringstream ss;
	ss << filehead << "_summary.csv";
	tree->print_adjacency_modes_summary_csv(ss.str());
	if(o->get_option<bool>(PRINT_VTKS)->get_value() || o->get_option<bool>(PRINT_CSVS)->get_value())  //print out modes in vtk file
	{
		int Nprint = p->get_param<int>(N_PRINT_KEY)->get_value();
		ss.clear();
		ss.str("");
		switch(o->get_option<char>(OUTPUT_SORT_KEY)->get_value())
		{
		case PRINT_DOMINANT_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_adjacency_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_adjacency_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
			} break;

		case PRINT_LARGEST_MODES:
			{
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_adjacency_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_adjacency_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
			} break;

		case PRINT_SMALLEST_MODES:
			{
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_adjacency_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_adjacency_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;

		case PRINT_ALL_MODES:
			{
				ss << filehead << "_dominant";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_adjacency_evectors_vtk(ss.str(), Nprint, DOMINANT_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_adjacency_evectors_csv(ss.str(), Nprint, DOMINANT_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_largest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_adjacency_evectors_vtk(ss.str(), Nprint, LARGEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_adjacency_evectors_csv(ss.str(), Nprint, LARGEST_SORT);
				}
				ss.clear();
				ss.str("");
				ss << filehead << "_smallest";
				if(o->get_option<bool>(PRINT_VTKS)->get_value())
				{
					tree->print_adjacency_evectors_vtk(ss.str(), Nprint, SMALLEST_SORT);
				}
				if(o->get_option<bool>(PRINT_CSVS)->get_value())
				{
					ss << ".csv";
					tree->print_adjacency_evectors_csv(ss.str(), Nprint, SMALLEST_SORT);
				}
			} break;
		}
	}
	return 0;
}

int main(int argc, char * argv[])
{
	std::shared_ptr<SpectralOptionList> options = std::make_shared<SpectralOptionList>();
	std::shared_ptr<SpectralParameterList> params = std::make_shared<SpectralParameterList>();

	read_options(argc, argv, options.get(), params.get());
	
	std::cout << "Initialising system...\n";
	
	//set-up tree
	std::shared_ptr<SpectralNetwork<>> tree = setup_tree(options.get(), params.get());

	//choose operator to compute
	switch(options->get_option<char>(OPERATOR_KEY)->get_value())
	{
	case MAURY:
		{
			if(compute_and_print_maury_spectrum(tree.get(), options.get(), params.get())) return 1;
		} break;

	case INTERNAL_LAPLACIAN:
		{
			compute_and_print_modified_laplacian_spectrum(tree.get(), options.get(), params.get());
		} break;

	case TRUNCATED_LAPLACIAN:
		{
			compute_and_print_truncated_laplacian_spectrum(tree.get(), options.get(), params.get());
		} break;

	case FULL_LAPLACIAN:
		{
			compute_and_print_graph_laplacian_spectrum(tree.get(), options.get(), params.get());
		} break;

	case ADJACENCY:
		{
			compute_and_print_adjacency_spectrum(tree.get(), options.get(), params.get());
		} break;

	default:
		{
			std::cerr << "Do not recognise operator option, assuming Maury.\n";
			compute_and_print_maury_spectrum(tree.get(), options.get(), params.get());
		} break;
	}

	return 0;
}