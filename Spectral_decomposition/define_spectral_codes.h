#ifndef DEFINE_SPECTRAL_CODES_H
#define DEFINE_SPECTRAL_CODES_H

#include <vector>
#include <string>

//option codes for tree structure
#define TREE_KEY "Tree"
#define TREE_OPTION_COUNT 3
#define ASYMM_MODEL 'a'     //asymmetric network
#define PERT_MODEL 'p'      //symmetric model with single pertubation
#define MODEL_FROM_FILE 'f'     //model from input file
const char Tree_option_list[] = {ASYMM_MODEL, PERT_MODEL, MODEL_FROM_FILE};
const std::string Tree_option_name_list[] = {"Asymm_model", "Perturbed_model", "Model_from_file"};

//option codes for edge weight
#define EDGE_WEIGHT_KEY "Weighting"
#define EDGE_WEIGHT_OPTION_COUNT 2
#define RESISTANCE_WEIGHTING 'r'
#define DIFFUSION_WEIGHTING 'd'
const char Edge_weight_option_list[] = {RESISTANCE_WEIGHTING, DIFFUSION_WEIGHTING};
const std::string Edge_weight_option_name_list[] = {"Resistance", "Diffusion"};

//options for operator to find spectrum of
#define OPERATOR_KEY "Operator"
#define	OPERATOR_OPTION_COUNT 5
#define MAURY 'm'
#define INTERNAL_LAPLACIAN 'i'
#define TRUNCATED_LAPLACIAN 't'
#define FULL_LAPLACIAN 'f'
#define ADJACENCY 'a'
const char Operator_option_list[] = {MAURY, INTERNAL_LAPLACIAN, TRUNCATED_LAPLACIAN, FULL_LAPLACIAN, ADJACENCY};
const std::string Operator_option_name_list[] = {"Maury", "Internal_Laplacian", "Truncated_Laplacian","Full_Laplacian","Adjacency"};

//options for order of calculation
#define ORDER_KEY "Order"
#define ORDER_OPTION_COUNT 2
#define LARGEST_MODES 'l'
#define SMALLEST_MODES 's'
const char Order_option_list[] = {LARGEST_MODES, SMALLEST_MODES};
const std::string Order_option_name_list[] = {"Largest", "Smallest"};

//option codes for cutoff mode
#define CUTOFF_MODE_KEY "Cutoff"
#define CUTOFF_OPTION_COUNT 3
#define NO_CUTOFF 'n'     //compute full specturm  
#define FRAC_CUTOFF 'f'      //cutoff at some fraction of total spectrum
#define VALUE_CUTOFF 'v'      //cutoff at some value
const char Cutoff_option_list[] = {NO_CUTOFF, FRAC_CUTOFF, VALUE_CUTOFF};
const std::string Cutoff_option_name_list[] = {"Complete", "Fraction", "Cutoff"};

//option codes for sorting of detailed outputs
#define OUTPUT_SORT_KEY "Sort"
#define OUTPUT_SORT_OPTION_COUNT 4
#define PRINT_LARGEST_MODES 'l'
#define PRINT_SMALLEST_MODES 's'
#define PRINT_DOMINANT_MODES 'd'
#define PRINT_ALL_MODES 'a'
const char Output_sort_option_list[] = {PRINT_LARGEST_MODES, PRINT_SMALLEST_MODES, 
	                                    PRINT_DOMINANT_MODES, PRINT_ALL_MODES};
const std::string Output_sort_option_name_list[] = {"Largest", "Smallest", "Dominant", "All"};

#define OPTIONS_FILE_EXT "options"
#define PARAMS_FILE_EXT "params"
#define NODE_FILE_EXT "nodes"
#define BRANCH_FILE_EXT "branches"
#define TERM_NODE_FILE_EXT "termnodes"

//parameter keys
#define N_GENERATIONS_KEY "Gens"
#define N_PRINT_KEY "Print"   //print detailed modes for x
#define PERTURBATION_GEN_KEY "Pert_gen"
#define ASYMMETRY_KEY "A"
#define PERTURBATION_FRACTION_KEY "Pert_frac"
#define CUTOFF_POINT_KEY "Cutoff"
#define SEED_RADIUS_KEY "Seed_radius"
#define SEED_LENGTH_KEY "Seed_length"
#define SCALE_FACTOR_KEY "Scale_factor"

//bool option keys
#define PRINT_VTKS "Print_vtk"
#define PRINT_CSVS "Print_csv"

#endif