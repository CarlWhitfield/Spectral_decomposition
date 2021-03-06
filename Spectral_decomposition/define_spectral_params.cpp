#include "define_spectral_codes.h"
#include "define_spectral_params.h"
#include "spectral_param_defaults.h"

SpectralParameterList::SpectralParameterList():ParameterList()
{
	using namespace inlist;
	//integer parameters
	auto ngens_par = std::make_shared<Parameter<int>>(N_GENERATIONS_DEFAULT, std::string(N_GENERATIONS_KEY));
	auto nprint_par = std::make_shared<Parameter<int>>(N_PRINT_DEFAULT, std::string(N_PRINT_KEY));
	auto pgen_par = std::make_shared<Parameter<int>>(PERTURBATION_GEN_DEFAULT, std::string(PERTURBATION_GEN_KEY));
	this->add(N_GENERATIONS_KEY, ngens_par);
	this->add(N_PRINT_KEY, nprint_par);
	this->add(PERTURBATION_GEN_KEY, pgen_par);

	//double params
	auto asymm_par = std::make_shared<Parameter<double>>(ASYMMETRY_DEFAULT, std::string(ASYMMETRY_KEY));
	auto pert_par = std::make_shared<Parameter<double>>(PERTURBATION_FRACTION_DEFAULT, std::string(PERTURBATION_FRACTION_KEY));
	auto cutoff_par = std::make_shared<Parameter<double>>(CUTOFF_POINT_DEFAULT, std::string(CUTOFF_POINT_KEY));
	auto srad_par = std::make_shared<Parameter<double>>(SEED_RADIUS_DEFAULT, std::string(SEED_RADIUS_KEY));
	auto slen_par = std::make_shared<Parameter<double>>(SEED_LENGTH_DEFAULT, std::string(SEED_LENGTH_KEY));
	auto sf_par = std::make_shared<Parameter<double>>(SCALE_FACTOR_DEFAULT, std::string(SCALE_FACTOR_KEY));
	this->add(ASYMMETRY_KEY, asymm_par);
	this->add(PERTURBATION_FRACTION_KEY, pert_par);
	this->add(CUTOFF_POINT_KEY, cutoff_par);
	this->add(SEED_RADIUS_KEY, srad_par);
	this->add(SEED_LENGTH_KEY, slen_par);
	this->add(SCALE_FACTOR_KEY, sf_par);
}

void SpectralParameterList::check_and_convert(SpectralOptionList *o)
{
	this->check_OK();

	//check for contradictory parameters
	if(o->get_option<char>(TREE_KEY)->get_value() != MODEL_FROM_FILE) //only used if model is not read from file
	{
		if(this->get_param<int>(N_GENERATIONS_KEY)->get_value() < 1 
		   || this->get_param<int>(N_GENERATIONS_KEY)->get_value() > 20)
		{
			std::cerr << "Invalid number of generations selected (min 1, max 20). Setting to default: " 
				      << N_GENERATIONS_DEFAULT << '\n';
			this->get_param<int>(N_GENERATIONS_KEY)->set_to_default();
		}

		if(o->get_option<char>(TREE_KEY)->get_value() == ASYMM_MODEL)    //only matters if building asymm model
		{
			if(this->get_param<double>(ASYMMETRY_KEY)->get_value() < 0 || this->get_param<double>(ASYMMETRY_KEY)->get_value() >= 1)
			{
				std::cerr << "Invalid asymmetry value selected (0 <= A < 1). Setting to default: " 
				      << ASYMMETRY_DEFAULT << '\n';
				this->get_param<int>(ASYMMETRY_KEY)->set_to_default();
			}
		}

		if(o->get_option<char>(TREE_KEY)->get_value() == PERT_MODEL)    //only matters if building perturbaed symmetric model
		{
			if(this->get_param<double>(PERTURBATION_FRACTION_KEY)->get_value() < -1.0)
			{
				std::cerr << "Invalid perturbation value selected (dr >= -1). Setting to default: " 
				      << PERTURBATION_FRACTION_DEFAULT << '\n';
				this->get_param<int>(PERTURBATION_FRACTION_KEY)->set_to_default();
			}
		}
	}

	if(o->get_option<char>(CUTOFF_MODE_KEY)->get_value() != NO_CUTOFF)  //otherwise unimportant
	{
		//fractional cutoff, f between 0 and 1 indicates fraction of spectrum to be computed
		if(o->get_option<char>(CUTOFF_MODE_KEY)->get_value() == FRAC_CUTOFF)
		{
			if(this->get_param<double>(CUTOFF_POINT_KEY)->get_value() <= 0 || this->get_param<double>(CUTOFF_POINT_KEY)->get_value() > 1)
			{
				std::cerr << "Incompatible value for cutoff fraction: f = " << this->get_param<double>(CUTOFF_POINT_KEY)->get_value() << " and 0 < f <= 1, setting to 1\n";
				this->get_param<double>(CUTOFF_POINT_KEY)->set_phys_value(1.0);
			}
		}
	}
}