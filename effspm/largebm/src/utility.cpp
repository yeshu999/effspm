#include "utility.hpp"
#include "build_mdd.hpp"
#include "load_inst.hpp"

#include <vector>
#include <ctime>
#include <algorithm>

namespace largebm {

bool check_parent(unsigned long long cur_arc,
                  unsigned long long str_pnt,
                  unsigned long long start,
                  std::vector<unsigned long long>& strpnt_vec)
{
    std::vector<unsigned long long> ancestors;
    
    unsigned long long cur_anct = Tree[cur_arc].anct;

    while (Tree[cur_anct].itmset > Tree[str_pnt].itmset) {
        if (Tree[cur_anct].item > 0)
            ancestors.push_back(cur_anct);
        cur_anct = Tree[cur_anct].anct;
    }

    if (Tree[cur_anct].itmset == Tree[str_pnt].itmset)
        return true;
    else {
        for (auto it = ancestors.rbegin(); it != ancestors.rend(); ++it) {
            for (unsigned long long i = start; i < strpnt_vec.size(); ++i) {
                if (strpnt_vec[i] == *it)
                    return true;
            }
        }
    }

    return false;
}

// return elapsed time in seconds
double give_time(std::clock_t ticks) {
    return static_cast<double>(ticks) / CLOCKS_PER_SEC;
}

} // namespace largebm
