#ifndef	DEFINE_SPECTRAL_OPTIONS_H
#define DEFINE_SPECTRAL_OPTIONS_H

#include "list_template.h"
#include "spectral_option_defaults.h"
#include "define_spectral_codes.h"
#include <string>

class SpectralOptionList: public inlist::OptionList<char,bool>
{
public:
	SpectralOptionList();

	std::string generate_output_name() const;
};


#endif