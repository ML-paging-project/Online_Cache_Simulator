#include "../include/RAM.h"
#include <set>
#include <fstream>
#include <stdlib.h>
#include <utility>
#include <random>

std::vector<uint64_t> working_set_simulator(uint32_t seed, bool print=false) {
    std::set<uint64_t> unique_pages;

    LRU_RAM *lru = new LRU_RAM(MEM_SIZE, PAGE_SIZE);

    // printf("Representative Workload\n");
    uint64_t working_size  = WORKING_SET * (MEM_SIZE / PAGE_SIZE);
    uint64_t leftover_size = WORKLOAD * (MEM_SIZE / PAGE_SIZE) - working_size;
    std::mt19937 rand(seed ^ 0xDEADBEEF); // create random number generator which has different randomness

    for(int i = 0; i < ACCESSES; i++) {
        double working_chance = rand() / (double) (0xFFFFFFFF);
        uint64_t v_addr = rand();
        if(working_chance <= LOCALITY)
            v_addr %= working_size;
        else
            v_addr = (v_addr % leftover_size) + working_size;

        // printf("Memory access %i\n", i);
        lru->memory_access(v_addr);

        unique_pages.insert(v_addr);
    }
    if (print) {
        printf("Out of %llu memory accesses with %lu unique virtual pages\n", ACCESSES, unique_pages.size());
        lru->printSuccessFunction();
    }
    std::vector<uint64_t> toret = lru->getSuccessFunction();
    delete lru;
    return toret;
}

std::vector<uint32_t> zipfian_simulator(bool print=false) {
    std::ifstream zipf_data(ZIPH_FILE);
    if (!zipf_data.is_open()) {
        printf("Failed to open zipf file\n");
        exit(EXIT_FAILURE);
    }

    std::set<uint64_t> unique_pages;

    LRU_RAM *lru = new LRU_RAM(MEM_SIZE, PAGE_SIZE);
    Clock_RAM *clk = new Clock_RAM(MEM_SIZE, PAGE_SIZE);
    uint64_t accesses = 0;

    printf("Zipfian Workload\n");
    std::string line;
    while(getline(zipf_data, line)) {
        uint64_t v_addr = std::stol(line);
        lru->memory_access(v_addr);
        clk->memory_access(v_addr);
        accesses++;

        unique_pages.insert(v_addr);
    }

    if (print) {
        printf("LRU:   "); lru->print();
        printf("CLOCK: "); clk->print();
        printf("Out of %llu memory accesses with %lu unique virtual pages\n", accesses, unique_pages.size());
    }
    zipf_data.close();

    return {lru->get_faults(), clk->get_faults()};
}


int main() {
    std::vector<uint64_t> lru_total;
    uint32_t trials = 10;

    std::mt19937 rand(SEED);
    for (int i = 0; i < trials; i++) {
        std::vector<uint64_t> result = working_set_simulator(rand());
        lru_total.resize(result.size());

        for (uint32_t page = 1; page < result.size(); page++)
        {
            lru_total[page] += result[page];
        }
        printf("# Trial: %i\r", i); std::fflush(stdout);
    }
    printf("# Final Statistics (%d)\n", lru_total.size());
    for (uint32_t page = 1; page < lru_total.size(); page++)
    {
        lru_total[page] /= trials;
        printf("%u:%lu\n", page, lru_total[page]);
    }
    // zipfian_simulator(true);
}
