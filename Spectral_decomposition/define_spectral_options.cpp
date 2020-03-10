#include<iostream>
#include "spectral_option_defaults.h"
#include "define_spectral_options.h"
#include "define_spectral_codes.h"
#include <boost/filesystem.hpp>

using namespace inlist;

//sets up default options (all arguments defined in define_spectral_codes.h)
SpectralOptionList::SpectralOptionList():OptionList()
{
	//mutliple choice options
	this->add(TREE_KEY, new Option<char>(TREE_DEFAULT, std::string(TREE_KEY), 
		      Tree_option_list, Tree_option_name_list, TREE_OPTION_COUNT));
	this->add(EDGE_WEIGHT_KEY, new Option<char>(EDGE_WEIGHT_DEFAULT, std::string(EDGE_WEIGHT_KEY), 
		      Edge_weight_option_list, Edge_weight_option_name_list, EDGE_WEIGHT_OPTION_COUNT));
	this->add(OPERATOR_KEY, new Option<char>(OPERATOR_DEFAULT, std::string(OPERATOR_KEY), 
		      Operator_option_list, Operator_option_name_list, OPERATOR_OPTION_COUNT));
	this->add(ORDER_KEY, new Option<char>(ORDER_DEFAULT, std::string(ORDER_KEY), 
		      Order_option_list, Order_option_name_list, ORDER_OPTION_COUNT));
	this->add(CUTOFF_MODE_KEY, new Option<char>(CUTOFF_MODE_DEFAULT, std::string(CUTOFF_MODE_KEY), 
		      Cutoff_option_list, Cutoff_option_name_list, CUTOFF_OPTION_COUNT));
	this->add(OUTPUT_SORT_KEY, new Option<char>(OUTPUT_SORT_DEFAULT, std::string(OUTPUT_SORT_KEY),
			  Output_sort_option_list, Output_sort_option_name_list, OUTPUT_SORT_OPTION_COUNT));
	//boolean options
	this->add(PRINT_VTKS, new Option<bool>(PRINT_VTK_DEFAULT, std::string(PRINT_VTKS)));
	this->add(PRINT_CSVS, new Option<bool>(PRINT_CSV_DEFAULT, std::string(PRINT_CSVS)));
}

std::string SpectralOptionList::generate_output_name() const
{
	std::stringstream ss;
	if(this->get_option<char>(TREE_KEY)->get_value() == MODEL_FROM_FILE)
	{
		if(this->filename_exists(BRANCH_FILE_EXT))
		{
			boost::filesystem::path p1(this->get_filename(BRANCH_FILE_EXT));
			ss << p1.filename().stem().string() << "_";
		}
		else
		{
			ss << this->get_option<char>(TREE_KEY)->get_value_name() << "_";
		}
	}
	else
	{
		if(this->filename_exists(PARAMS_FILE_EXT))
		{ 
			boost::filesystem::path p1(this->get_filename(PARAMS_FILE_EXT));
			ss << p1.filename().stem().string() << "_";
		}
		else
		{
			ss << "Default_params_";
		}
	}
	ss << this->get_option<char>(EDGE_WEIGHT_KEY)->get_value_name() << "_weigths_";
	ss << this->get_option<char>(OPERATOR_KEY)->get_value_name() << "_";
	//ss << this->get_option<char>(ORDER_KEY)->get_value_name() << "_";
	ss << this->get_option<char>(CUTOFF_MODE_KEY)->get_value_name();

	return ss.str();
}