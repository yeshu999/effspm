// load_inst.hpp
#pragma once

#include <vector>
#include <string>
#include <ctime>      // for clock_t
#include <algorithm>  // for std::max, std::abs
#include "freq_miner.hpp"   // for Pattern

// —————————————————————————————————————————————————————
// Globals mirrored from main.cpp / load_inst.cpp
// —————————————————————————————————————————————————————

// the raw dataset (either loaded from file or injected from Python)
extern std::vector<std::vector<int>> items;

// the DFS frontier of partial patterns
extern std::vector<Pattern>           DFS;

// output filename (if b_write==true)
extern std::string    out_file;

// runtime flags
extern bool           b_disp;      // verbose console print
extern bool           b_write;     // append to out_file
extern bool           use_dic;     // use dictionary remapping
extern bool           use_list;    // list‐based optimization
extern bool           pre_pro;     // run the “-preproc” step

// dataset statistics
extern unsigned int   M;           // max sequence length
extern unsigned long long E;       // total number of “entries” across all sequences
extern unsigned int   L;           // number of distinct items
extern unsigned long long N;       // number of sequences

// support threshold (absolute count)
extern unsigned long long theta;

// time‐limit and wall‐clock start
extern unsigned int   time_limit;  // in seconds
extern std::clock_t   start_time;  // for give_time()

// —————————————————————————————————————————————————————
// load_instance: file‐based loader (same signature as before)
// —————————————————————————————————————————————————————
bool Load_instance(std::string &items_file, double thresh);
