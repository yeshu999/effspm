#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include "load_inst.hpp"
#include "build_mdd.hpp"
#include "utility.hpp"
#include "freq_miner.hpp"

namespace btminer {

// These variables are declared in utility.hpp and defined in utility.cpp
extern std::string out_file;
extern bool b_disp, b_write, use_dic, just_build;
extern clock_t start_time;

            // This one is only used in main
       // Local to main
int time_limit = 30 * 3600;       // Local to main
std::string folder;

int main(int argc, char* argv[]) {
    std::string VV, attr;
    double thresh = 0;

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-' || isdigit(argv[i][1]))
            continue;
        else if (strcmp(argv[i], "-thr") == 0)
            thresh = std::stod(argv[i + 1]);
        else if (strcmp(argv[i], "-file") == 0)
            VV = argv[i + 1];
        else if (strcmp(argv[i], "-N_mult") == 0)
            N_mult = std::stoi(argv[i + 1]);
        else if (strcmp(argv[i], "-M_mult") == 0)
            M_mult = std::stoi(argv[i + 1]);
        else if (strcmp(argv[i], "-time") == 0)
            time_limit = std::stoi(argv[i + 1]);
        else if (strcmp(argv[i], "-jbuild") == 0)
            just_build = true;
        else if (strcmp(argv[i], "-folder") == 0)
            folder = argv[i + 1];
        else if (strcmp(argv[i], "-npre") == 0)
            pre_pro = false;
        else if (strcmp(argv[i], "-dic") == 0)
            use_dic = true;
        else if (strcmp(argv[i], "-out") == 0) {
            if (i + 1 == argc || argv[i + 1][0] == '-')
                b_disp = true;
            else if (argv[i + 1][0] == '+') {
                b_disp = true;
                b_write = true;
                if (strlen(argv[i + 1]) > 1) {
                    out_file = argv[i + 1];
                    out_file = out_file.substr(1);
                } else {
                    out_file = VV;
                }
            } else {
                b_write = true;
                out_file = argv[i + 1];
            }
        } else {
            std::cout << "Command " << argv[i] << " not recognized and skipped.\n";
        }
    }

    std::cout << "\n********************** " << VV << " N_mult: " << N_mult << " M_mult: " << M_mult << "**********************\n";

    std::string item_file = folder + VV + ".txt";

    std::cout << "loading instances...\n";
    start_time = clock();

    if (!Load_instance(item_file, thresh)) {
        std::cout << "Files invalid, exiting.\n";
        std::cin.get();
        return 0;
    }

    if (!just_build && give_time(clock() - start_time) < time_limit) {
        Freq_miner();
        if (give_time(clock() - start_time) >= time_limit)
            std::cout << "TIME LIMIT REACHED\n";
        std::cout << "Mining Complete\n\nFound a total of " << num_patt << " patterns\n";
        std::cout << "\nTotal CPU time " << give_time(clock() - start_time) << " seconds\n\n";
    }

    return 0;
}

} // namespace btminer
