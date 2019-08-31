#include<list_template.h>
#include"spectral_operations.h"

#define LARGEST 0
#define SMALLEST 1
#define BOTH 2

const size_t Nprint = 10;  //number of dominant modes to print
const double frac = 1.0;  //fraction of modes to compute
const int opt = LARGEST;
//const size_t Ngens = 5;

int main(int argc, char * argv[])
{
	inlist::OptionList<bool,char> *options = new inlist::OptionList<bool,char>();
	std::vector<std::string> extensions;
	extensions.resize(3);
	extensions[0] = NODE_FILE_EXT;
	extensions[1] = BRANCH_FILE_EXT;
	extensions[2] = TERM_NODE_FILE_EXT;
	options->get_filenames_from_args(extensions, argc, argv);
	
	std::cout << "Initialising system...\n";
	SpectralNetwork<RESISTANCE_NETWORK> tree(options->get_filename(NODE_FILE_EXT), 
		   options->get_filename(BRANCH_FILE_EXT), options->get_filename(TERM_NODE_FILE_EXT));
	
	//SpectralNetwork<RESISTANCE_NETWORK> tree(Ngens, 1.0, 6.0, 3.0, 0.1);

	std::cout << "Computing Laplacian Spectrum...\n";
	size_t Nsmall = 0, Nlarge = 0;
	size_t Nmodes = tree.count_nodes()-tree.count_term_nodes()+1;
	if(frac==1)
	{
		tree.compute_full_truncated_laplacian_spectrum();
	}
	else
	{
		switch(opt)
		{
		case LARGEST:
			{
				Nlarge = size_t(frac*Nmodes);
			} break;
		case SMALLEST:
			{
				Nsmall = size_t(frac*Nmodes);
			} break;
		case BOTH:
			{
				Nlarge = size_t(0.5*frac*Nmodes);
				Nsmall = Nlarge;
			} break;
		}
		tree.compute_truncated_laplacian_spectrum(Nsmall, Nlarge);
	}
	

	std::cout << "Computing Maury Spectrum...\n";
	Nmodes = tree.count_term_nodes();
	if(frac==1)
	{
		tree.compute_full_maury_spectrum();
	}
	else
	{
		switch(opt)
		{
		case LARGEST:
			{
				Nlarge = size_t(frac*Nmodes);
			} break;
		case SMALLEST:
			{
				Nsmall = size_t(frac*Nmodes);
			} break;
		case BOTH:
			{
				Nlarge = size_t(0.5*frac*Nmodes);
				Nsmall = Nlarge;
			} break;
		}
		tree.compute_maury_spectrum(Nmall,Nlarge);
	}
	

	std::cout << "Outputting...\n";
	tree.print_dominant_laplacian_evectors_vtk("Laplacian_evecs", Nprint);
	tree.print_dominant_maury_evectors_vtk("Maury_evecs", Nprint);
	tree.print_unit_pressure_drop_solutions_vtk("Fluxes");
	tree.print_laplacian_modes_csv("Laplacian_modes.csv");
	tree.print_maury_modes_csv("Maury_modes.csv");

	return 0;
}