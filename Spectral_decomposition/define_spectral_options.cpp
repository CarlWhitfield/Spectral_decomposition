#include<iostream>
#include "spectral_option_defaults.h"
#include "define_spectral_options.h"
#include "define_spectral_codes.h"
#include <boost/filesystem.hpp>


//sets up default options (all arguments defined in define_spectral_codes.h)
sgt::SpectralOptionList::SpectralOptionList():OptionList()
{
	//mutliple choice options
	using namespace inlist;
	auto tree_opt = std::make_shared<Option<char>>(TREE_DEFAULT, std::string(TREE_KEY), 
												   sgt::Tree_option_list, sgt::Tree_option_name_list,
												   TREE_OPTION_COUNT);
	auto weight_opt = std::make_shared<Option<char>>(EDGE_WEIGHT_DEFAULT, std::string(EDGE_WEIGHT_KEY), 
		      sgt::Edge_weight_option_list, sgt::Edge_weight_option_name_list, EDGE_WEIGHT_OPTION_COUNT);
	auto op_opt = std::make_shared<Option<char>>(OPERATOR_DEFAULT, std::string(OPERATOR_KEY), 
		      sgt::Operator_option_list, sgt::Operator_option_name_list, OPERATOR_OPTION_COUNT);
	auto order_opt = std::make_shared<Option<char>>(ORDER_DEFAULT, std::string(ORDER_KEY), 
		      sgt::Order_option_list, sgt::Order_option_name_list, ORDER_OPTION_COUNT);
	auto cutoff_opt = std::make_shared<Option<char>>(CUTOFF_MODE_DEFAULT, std::string(CUTOFF_MODE_KEY), 
		      sgt::Cutoff_option_list, sgt::Cutoff_option_name_list, CUTOFF_OPTION_COUNT);
	auto out_opt = std::make_shared<Option<char>>(OUTPUT_SORT_DEFAULT, std::string(OUTPUT_SORT_KEY),
			  sgt::Output_sort_option_list, sgt::Output_sort_option_name_list, OUTPUT_SORT_OPTION_COUNT);
	this->add(TREE_KEY, tree_opt);
	this->add(EDGE_WEIGHT_KEY, weight_opt);
	this->add(OPERATOR_KEY, op_opt);
	this->add(ORDER_KEY, order_opt);
	this->add(CUTOFF_MODE_KEY, cutoff_opt);
	this->add(OUTPUT_SORT_KEY, out_opt);

	//boolean options
	auto vtk_opt = std::make_shared<Option<bool>>(PRINT_VTK_DEFAULT, std::string(PRINT_VTKS));
	auto csv_opt = std::make_shared<Option<bool>>(PRINT_CSV_DEFAULT, std::string(PRINT_CSVS));
	this->add(PRINT_VTKS, vtk_opt);
	this->add(PRINT_CSVS, csv_opt);
}

std::string sgt::SpectralOptionList::generate_output_name() const
{
	using namespace inlist;
	std::stringstream ss;
	if(this->get_option_value<char>(TREE_KEY) == MODEL_FROM_FILE)
	{
		if(this->filename_exists(BRANCH_FILE_EXT))
		{
			boost::filesystem::path p1(this->get_filename(BRANCH_FILE_EXT));
			ss << p1.filename().stem().string() << "_";
		}
		else
		{
			ss << this->get_option_name(TREE_KEY) << "_";
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
	ss << this->get_option_name(EDGE_WEIGHT_KEY) << "_weigths_";
	ss << this->get_option_name(OPERATOR_KEY) << "_";
	//ss << this->get_option<char>(ORDER_KEY)->get_value_name() << "_";
	ss << this->get_option_name(CUTOFF_MODE_KEY);

	return ss.str();
}
