#include <stdlib.h>

#include <cinttypes>
#include <fstream>
#include <random>
#include <set>
#include <utility>
#include <chrono>

#include "IncrementAndKill.h"
#include "LruSizesSim.h"

uint64_t get_next_addr(std::mt19937& gen)
{
  // printf("Representative Workload\n");
  uint64_t working_size = WORKING_SET;
  uint64_t leftover_size = WORKLOAD * working_size;
  double working_chance = gen() / (double)(0xFFFFFFFF);
  uint64_t v_addr = gen();
  if (working_chance <= LOCALITY)
    v_addr %= working_size;
  else
    v_addr = (v_addr % leftover_size) + working_size;
  return v_addr;
}


/*
 * Runs a random sequence of page requests that follow a working_set
 * distribution A majority of accesses are made to a somewhat small working_set
 * while the rest are made randomly to a much larger amount of memory seed: the
 * seed to the random number generator print:    if true print out the results
 * of the simulation returns   the success function
 */
std::vector<std::vector<uint64_t>> working_set_simulator(uint32_t seed, bool print = false) {
  std::set<uint64_t> unique_pages;
    
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;

  LruSizesSim *lru = new LruSizesSim();
  IncrementAndKill *iak = new IncrementAndKill();

  // Order Statistic LRU (stack distance)
  std::mt19937 rand(seed);  // create random number generator
  auto start = high_resolution_clock::now();
  for (uint64_t i = 0; i < ACCESSES; i++) {
    lru->memory_access(get_next_addr(rand));
  }
  std::vector<uint64_t> lru_success = lru->get_success_function();
  auto lru_time =  duration<double>(high_resolution_clock::now() - start).count();
  
  // Increment And Kill
  rand.seed(seed);  // create random number generator
  start = high_resolution_clock::now();
  for (uint64_t i = 0; i < ACCESSES; i++) {
    iak->memory_access(get_next_addr(rand));
  }
  std::vector<uint64_t> iak_success = iak->get_success_function();
  auto iak_time = duration<double>(high_resolution_clock::now() - start).count();

  if (print) {
    // Do this for stats
    start = high_resolution_clock::now();
    rand.seed(seed);  // create random number generator
    for (uint64_t i = 0; i < ACCESSES; i++) {
      unique_pages.insert(get_next_addr(rand));
    }
    auto vec_time = duration<double>(high_resolution_clock::now() - start).count();

    // print the success function if requested
    printf("Out of %" PRIu64 " memory accesses with %lu unique virtual pages\n",
        ACCESSES, unique_pages.size());
    lru->print_success_function();
    iak->print_success_function();
    std::cerr << "VEC TIME: " << vec_time << std::endl;
  }
  std::cerr << "LRU TIME: " << lru_time << std::endl;
  std::cerr << "IAK TIME: " << iak_time << std::endl;
  delete lru;
  delete iak;
  return {lru_success, iak_success};
}

// check that the results of two different simulators are the same
// IMPORTANT: vectors may be of different sizes
bool check_equivalent(std::vector<uint64_t> vec_1, std::vector<uint64_t> vec_2) {
  size_t i = 0;
  uint64_t last_elm = vec_1[0];
  while (i < vec_1.size() && i < vec_2.size()) {
    if (vec_1[i] != vec_2[i]) return false;

    last_elm = vec_1[i++];
  }

  // assert that any remaining elements in either vector
  // are identical to last_elm
  for (size_t j = i; j < vec_1.size(); j++)
    if (vec_1[j] != last_elm) return false;

  for (size_t j = i; j < vec_2.size(); j++)
    if (vec_2[j] != last_elm) return false;

  return true;
}

int main() {
  auto results = working_set_simulator(SEED, false);
  auto lru_results = results[0];
  auto iak_results = results[1];
  bool eq = check_equivalent(iak_results, lru_results);
  std::cerr << "Are results equivalent?: " << (eq? "yes" : "no") << std::endl;
  std::cerr << "Sizes (lru, iak): " << lru_results.size() << ", " << iak_results.size() << std::endl;
}
