# Spectral_decomposition

Command-line tool for computing spectral properties of operators associated with a given airway network. Command line arguments are as follows:

- .options file (plain text): Contains code options (see example). If not given default options will be used (see spectral_option_defaults.h).
- .params file (plain text): Contains code parameters (see example). If not given default params will be used (see spectral_param_defaults.h).

If using the Tree f (from file) option, then files defining the tree network need to be provided at the command line too. These are given in plain text format with the following extensions:

- .nodes: minimum 4 columns: first column is the node index (integer), next 3 are (x,y,z) coordinates of the node in mm (floating pt)
- .branches: minimum 4 columns: col 1: branch index (int), node in index (int), node out index (int), radius mm (floating pt)
- .termnodes: min 1 column: list of node indices that terminate the tree (int).

Code options are:
      Tree: Options are 'a' asymmetric, 'p' perturbation, 'f' from file. First two build asymmetric tree with fixed number of gens. Perturbation option changes the radius of a single airway.

      Weighting: Options are 'r' = resistance or 'd' = diffusion. Determines edge weighting, resistance is Poiseuille r = 8 nu l / Pi r^4 (for viscosity nu) edge weight is w = 1/r or diffusion w = Pi D r^2 / l (diffusivity D)

      Operator: Options are 'm' = maury matrix, 'i' = internal laplacian, 't' truncated laplacian, 'f' full laplacian, 'a' adjacency. Determines the operator of which to find the spectrum.

      Order: Options are 'l' largest, 's' smallest. Order in which to compute the modes (from largest or smallest)

      Cutoff: Options are 'n' no cutoff (full spectrum), 'f' fraction of spectrum, 'v' cutoff by eigenvalue. Determines how calculation is terminated

      Sort: Options are 'l' largest evalue, 's' smallest evalue, 'd' dominant, 'a' print all. Determines how eigenmodes to be printed in full are selected.

      Print_vtk: Print full modes in vtk format? 't' of 'f'

      Print_csv: Print full modes in csv format? 't' of 'f'
     
Code parameters are:
      Print: Integer parameter. No. of eigenmodes to print in full
      
      Gens: Integer parameter. If tree option 'a' or 'p', no. of tree generations.
      
      Pert_gen: Integer parameter. If tree option 'p',  generation to apply perturbation, must be less than Gens
      
      A: Floating pt parameter 0 <= A <= 1. If tree option 'a' or 'p', asymmetric branching factor
      
      Pert_frac: Floating pt parameter. If tree option 'p', relative perturbation to airway resistance (e.g. 0.5 is 50% increase in resistance), cannot be <= -1.
      Cutoff: Floating pt parameter. If Cutoff option 'f' - fraction of spectrum to keep 0 <= Cutoff <= 1. If 'v' - eigenvalue to cutoff at.
      Seed_radius: Floating pt parameter. If tree option 'a' or 'p', radius of first airway arbitrary units
      Seed_length: Floating pt parameter. If tree option 'a' or 'p', length of first airway arbitrary units
      Scale_factor: Floating pt parameter. If tree option 'a' or 'p', scale factor for branching, default 3 (Murray's law).  
