#ifndef	DEFINE_SPECTRAL_PARAMS_H
#define DEFINE_SPECTRAL_PARAMS_H

#include "list_template.h"
#include "define_spectral_options.h"
#include "spectral_param_defaults.h"
#include <string>

namespace sgt
{
	class SpectralParameterList: public inlist::ParameterList<int, double>
	{
	public:
		SpectralParameterList();   //constructor for spectral calc
		void check_and_convert(SpectralOptionList *o);
	};
}

#endif
