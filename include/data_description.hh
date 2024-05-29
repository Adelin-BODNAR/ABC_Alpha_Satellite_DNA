#ifndef DATA_DESCRIPTION_H
#define DATA_DESCRIPTION_H

#include <string> //std::string
#include <vector>

std::string generate_summary_stats (int nb_orders_tested, int nb_similarities_calculated_per_order, int seed, std::vector<std::string>* list_sequences );

#endif