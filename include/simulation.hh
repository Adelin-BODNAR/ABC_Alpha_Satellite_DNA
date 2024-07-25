#ifndef SIMULATION_H
#define SIMULATION_H

#include <list>            // For std::list
#include <memory>          // For std::shared_ptr
#include <Bpp/Seq/Sequence.h>  // For bpp::Sequence
#include <Bpp/Seq/Alphabet/AlphabetTools.h> // For bpp::Alphabet
#include <Bpp/Phyl/Tree/TreeTemplate.h>  // For bpp::TreeTemplate


void amplification_simulation (int max_size, double amplification_rate, double monomeric_prob, double alpha_amplification_size, double beta_amplification_size, int min_amplification_size, int max_amplification_size, double alpha_HOR_order, double beta_HOR_order, int min_HOR_order, int max_HOR_order, bpp::Sequence ancestor_monomer, std::shared_ptr<const bpp::Alphabet> alpha, double substitution_rate, double transition_transversion_ratio, int seed, int verbose, bpp::TreeTemplate<bpp::Node>** tree_ptr, std::list<bpp::Node*>* monomers_dataset_ptr );


#endif
