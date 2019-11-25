#include<list_template.h>
#include"spectral_operations.h"

#define LARGEST 0
#define SMALLEST 1
#define BOTH 2

const size_t Nprint = 100;  //number of dominant modes to print
const double frac = 1.0;  //fraction of modes to compute
const int opt = LARGEST;

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

	tree.print_unit_pressure_drop_solutions_vtk("Fluxes");

	std::cout << "Computing Laplacian Spectrum...\n";
	int Nsmall = 0, Nlarge = 0;
	int Nmodes = int(tree.count_nodes()-tree.count_term_nodes()+1);
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
				Nlarge = int(frac*Nmodes);
			} break;
		case SMALLEST:
			{
				Nsmall = int(frac*Nmodes);
			} break;
		case BOTH:
			{
				Nlarge = int(0.5*frac*Nmodes);
				Nsmall = Nlarge;
			} break;
		}
		tree.compute_truncated_laplacian_spectrum(Nsmall, Nlarge);
	}
	std::cout << "Outputting...\n";
	tree.print_laplacian_evectors_vtk("Laplacian_dominant_evecs", Nprint, DOMINANT_SORT);
	tree.print_laplacian_evectors_vtk("Laplacian_largest_evecs", Nprint, LARGEST_SORT);
	tree.print_laplacian_evectors_csv("Laplacian_dominant_evecs.csv", Nprint, DOMINANT_SORT);
	tree.print_laplacian_evectors_csv("Laplacian_largest_evecs.csv", Nprint, LARGEST_SORT);

	std::cout << "Computing Maury Spectrum...\n";
	Nmodes = int(tree.count_term_nodes());
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
				Nlarge = int(frac*Nmodes);
			} break;
		case SMALLEST:
			{
				Nsmall = int(frac*Nmodes);
			} break;
		case BOTH:
			{
				Nlarge = int(0.5*frac*Nmodes);
				Nsmall = Nlarge;
			} break;
		}
		tree.compute_maury_spectrum(Nsmall, Nlarge);
	}
	

	std::cout << "Outputting...\n";
	tree.print_maury_evectors_vtk("Maury_dominant_evecs", Nprint, DOMINANT_SORT);
	tree.print_maury_evectors_vtk("Maury_largest_evecs", Nprint, LARGEST_SORT);
	tree.print_maury_evectors_csv("Maury_dominant_evecs.csv", Nprint, DOMINANT_SORT);
	tree.print_maury_evectors_csv("Maury_largest_evecs.csv", Nprint, LARGEST_SORT);	

	/*std::cout << "Outputting all modes (may crash depending on size)\n";
	tree.print_laplacian_modes_csv("Laplacian_modes.csv");
	tree.print_maury_modes_csv("Maury_modes.csv");*/

	return 0;
}