#pragma once

#include "CacheSim.h"

#include <cassert>

enum OpType {Null, Increment, Kill, Uninit};

// A single operation, such as increment or kill
class Op {
 private:
  OpType type;     // Do we increment or kill?
  uint64_t start;  // Increment start index / Kill index
  uint64_t end;    // Increment End index
  uint64_t r;      // Increment amount

 public:
  // create an increment
  Op(uint64_t start, uint64_t end)
      : type(Increment), start(start), end(end), r(1){};

  // create a kill
  Op(uint64_t kill) : type(Kill), start(kill){};
  
  // Uninitialized. Used to parallelize making a vector of this without push_back
  Op() : type(Uninit){};

  // create a new Op by projecting another one
  Op(const Op& oth_op, uint64_t proj_start, uint64_t proj_end);

  // This enforces correct merging of operators
  Op& operator+=(Op& oth) {
    assert(oth.type == Increment);
    assert(type == Increment);
    assert(start == oth.start);
    assert(end == oth.end);

    r += oth.r;
    return *this;
  }

  // is this operation passive in the projection defined by
  // proj_start and proj_end
  bool is_passive(uint64_t proj_start, uint64_t proj_end) {
    return type == Increment && start <= proj_start && end >= proj_end;
  }
  
  
  OpType get_type() {assert(type != Uninit); return type;}
  uint64_t get_r()  {assert(type != Uninit); return r;}
  uint64_t get_start() {assert(type != Uninit); return start;}
  uint64_t get_end() {assert(type == Increment); return end;}
};

// A sequence of operators defined by a projection
class ProjSequence {
public:
  std::vector<Op> op_seq;
  uint64_t start;
  uint64_t end;

  // Initialize an empty projection with bounds
  ProjSequence(uint64_t start, uint64_t end) : start(start), end(end) {};

  // We project and merge here. 
  void add_op(Op& new_op) {
    Op proj_op = Op(new_op, start, end);
    if (proj_op.get_type() == Null) return;
    if (op_seq.size() == 0) {
      op_seq.push_back(proj_op);
      return;
    }

    Op& last_op = op_seq[op_seq.size() - 1];
    // If either is not passive, then we cannot merge. We must add it.
    if (!proj_op.is_passive(start, end) || !last_op.is_passive(start, end)) {
      op_seq.push_back(proj_op);
      return;
    }
    // merge with the last op in sequence
    last_op += proj_op;
  }
};

// Implements the IncrementAndKill algorithm
class IncrementAndKill: public CacheSim {
  using tuple = std::pair<uint64_t, uint64_t>;
  private:
    // A vector of all requests
    std::vector<tuple> requests;

    /* A vector of tuples for previous and next
     * Previous defines the last instance of a page
     * next defines the next instance of a page
     */
    std::vector<tuple> prevnext;

    /* This converts the requests into the previous and next vectors
     * Requests is copied, not modified.
     * Precondition: requests must be properly populated.
     */
    void calculate_prevnext();

    /* Returns the distance vector calculated from prevnext.
     * Precondition: prevnext must be properly populated.
     */
    std::vector<uint64_t> get_distance_vector();
    // Shortcut to access prev in prevnext.
    uint64_t& prev(uint64_t i) {return prevnext[i].first;}
    // Shortcut to access next in prevnext.
    uint64_t& next(uint64_t i) {return prevnext[i].second;}
    
    /* Helper function to get_distance_vector.
     * Recursively (and in parallel) populates the distance vector if the 
     * projection is small enough, or calls itself with smaller projections otherwise.
     */
    void do_projections(std::vector<uint64_t>& distance_vector, ProjSequence seq);
    
    /* Helper function for solving a projected sequence that has a size > 1
     * Therefore our base case can be larger so we avoid unnecessary overhead
     */
    void brute_force_alg(std::vector<uint64_t>& distance_vector, ProjSequence seq);
  public:
    // Logs a memory access to simulate. The order this function is called in matters.
    void memory_access(uint64_t addr);
    /* Returns the success function.
     * Does *a lot* of work.
     * When calling print_success_function, the answer is re-computed.
     */
    std::vector<uint64_t> get_success_function();
    IncrementAndKill() = default;
    ~IncrementAndKill() = default;
};

/*
Two distributions, one purely random and the other ziphian
Zipfian: Pick ith item with probability 1 / (i ln n) or 1 / (i H(n))
Pick random number between 1 and log n and then pick a number in 2^k-1 and 2^k

*/
