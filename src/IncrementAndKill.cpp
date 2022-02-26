#include "IncrementAndKill.h"

#include <algorithm>
#include <queue>
#include <omp.h>

void IncrementAndKill::memory_access(uint64_t addr) {
  requests.push_back({addr, access_number++});
}

void IncrementAndKill::calculate_prevnext() {
  // put all requests of the same addr next to each other
  // then order those by access_number
  auto requestcopy = requests;

  std::sort(requestcopy.begin(), requestcopy.end());

  prevnext.resize(requestcopy.size() + 1);

#pragma omp parallel for
  for (uint64_t i = 0; i < requestcopy.size(); i++) {
    auto [addr, access_num] = requestcopy[i];
    auto [last_addr, last_access_num] = i == 0 ? tuple(0, 0): requestcopy[i-1];

    // Using last, check if previous sorted access is the same
    if (last_access_num > 0 && addr == last_addr) {
      prev(access_num) = last_access_num;  // Point access to previous
      next(last_access_num) = access_num;  // previous access to this
    } else {
      prev(access_num) = 0;  // last is different so prev = 0
    }

    // Preemptively point this one's next access to the end
    next(access_num) = requestcopy.size() + 1;
  }
}

std::vector<uint64_t> IncrementAndKill::get_distance_vector() {
  std::vector<uint64_t> distance_vector(requests.size() + 1);

  // Generate the list of operations
  std::vector<Op> operations(2*requests.size());
  //TODO: This was probably better using push_back and reserve (single threaded)
#pragma omp parallel for
  for (uint64_t i = 1; i <= requests.size(); i++) {
    operations[2*i-2] = Op(prev(i) + 1, i - 1);        // Increment(prev(i)+1, i-1, 1)
    operations[2*i-1] = Op(prev(i));  // Kill(prev(i))
  }

  // begin the recursive process
  ProjSequence init_seq(1, requests.size());
  init_seq.op_seq = operations;

  // We want to spin up a bunch of threads, but only start with 1.
  // More will be added in by do_projections.
#pragma omp parallel
#pragma omp single
  do_projections(distance_vector, std::move(init_seq));

  return distance_vector;
}

void IncrementAndKill::brute_force_alg(std::vector<uint64_t>& distance_vector, ProjSequence cur) {
  // loop through each operation in the sequence and apply it to the distance vector
  std::vector<bool> killed(cur.end - cur.start + 1);
  for (auto op : cur.op_seq) {
    assert(op.get_type() == Kill || op.get_type() == Increment);
    
    if (op.get_type() == Increment) {
      // For each item in the range increment the distance vector
      for (uint64_t i = op.get_start(); i <= op.get_end(); i++) {
        if (!killed[i - cur.start])
          distance_vector[i] += op.get_r();
      }
    } 
    else {
      // std::cerr << "Kill the op" << std::endl;
      // Kill the item in question
      killed[op.get_start() - cur.start] = true;
    }
  }
}

//recursively (and in parallel) perform all the projections
void IncrementAndKill::do_projections(std::vector<uint64_t>& distance_vector, ProjSequence cur) {
  // base case -- use the brute force algorithm
  // no need for locking as projections have unique bounds
  if (cur.end - cur.start < 64) {
    brute_force_alg(distance_vector, cur);
  }
  else {
    uint64_t dist = cur.end - cur.start;
    uint64_t mid = (dist) / 2 + cur.start;

    // generate projected sequence for first half
    ProjSequence fst_half(cur.start, mid);
    for (uint64_t i = 0; i < cur.op_seq.size(); i++) {
      fst_half.add_op(cur.op_seq[i]);
    }
#pragma omp task shared(distance_vector) mergeable final(dist <= 1024) 
    do_projections(distance_vector, std::move(fst_half));

    // generate projected sequence for second half
    ProjSequence snd_half(mid + 1, cur.end);
    for (uint64_t i = 0; i < cur.op_seq.size(); i++) {
      snd_half.add_op(cur.op_seq[i]);
    }
    do_projections(distance_vector, std::move(snd_half));
  }
}

std::vector<uint64_t> IncrementAndKill::get_success_function() {
  calculate_prevnext();
  auto distances = get_distance_vector();

  // a point representation of successes
  std::vector<uint64_t> success(distances.size());
  for (uint64_t i = 1; i < distances.size(); i++) {
    if (prev(i + 1) != 0) success[distances[prev(i + 1)] + 1]++;
  }
  // integrate
  uint64_t running_count = 0;
  for (uint64_t i = 1; i < success.size(); i++) {
    running_count += success[i];
    success[i] = running_count;
  }
  return success;
}

// Create a new Operation by projecting another
Op::Op(const Op& oth_op, uint64_t proj_start, uint64_t proj_end) {
  // check if Op becomes Null
  // Increments are Null if we end before the start OR have a 'bad' interval
  // (end before our start)
  if (oth_op.type == Increment &&
      (oth_op.end < oth_op.start || oth_op.end < proj_start ||
       oth_op.start > proj_end)) {
    type = Null;
    return;
  }

  // kills do not shrink unless out of bounds
  if (oth_op.type == Kill) {
    if (oth_op.start < proj_start || oth_op.start > proj_end) {
      type = Null;
    } else {
      type = Kill;
      start = oth_op.start;
    }
    return;
  }

  // shrink the operation by the projection
  uint64_t n_start = oth_op.start;
  uint64_t n_end = oth_op.end;
  if (proj_start > n_start) n_start = proj_start;
  if (proj_end < n_end) n_end = proj_end;

  type = Increment;
  start = n_start;
  end = n_end;
  r = oth_op.r;
}
