// effspm/load_inst.hpp
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <ctime>           // for clock_t

#include "freq_miner.hpp"  // bring in Pattern

using namespace std;

// Load from disk (file-path) overload
bool Load_instance(string &items_file, double thresh);

// storage & globals used by both C++ CLI and Python binding:
extern vector<vector<int>> items;
extern vector<Pattern>     DFS;
extern vector<int>         item_dic;

extern string              out_file;
extern bool                b_disp, b_write, use_dic, use_list, pre_pro;

extern unsigned int        M, L, time_limit;
extern unsigned long long  N, E, theta;   // E = total entries

extern clock_t             start_time;
