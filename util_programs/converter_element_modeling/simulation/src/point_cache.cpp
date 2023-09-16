#include <iostream>
#include <vector>
#include <list>
#include <thread>
#include <mutex>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "bin.hpp"
#include "point_cache.hpp"

//PointCache::PointCache() { tl_vecs.reserve(4); } // Guess: 4 worker threads

PointCache::~PointCache() {
  // Acquire the mutex to make sure no one is still using
  Lock();
  // Check if all issued vectors have been returned
  if (issued_vecs.size())
    std::cerr << "WARNING: vector not returned to PointCache!\n";
  Unlock();
}

bool PointCache::is_issued(std::vector<GeantParticle>* v) const {
  return std::find(issued_vecs.begin(), issued_vecs.end(), v) != issued_vecs.end();
}

void PointCache::Lock() {
  vec_mutex.lock();
}
void PointCache::Unlock() { vec_mutex.unlock(); }

std::vector<GeantParticle>* PointCache::GetVec() {
  for (auto& ele : tl_vecs) {
    // Look for an un-issued vector
    if (is_issued(&ele)) continue;
    // ele not found in issued_vecs -> can be issued
    issued_vecs.push_back(&ele);
    return &ele;
  }
  // No un-issued vectors -> make a new one
  tl_vecs.push_back({});
  issued_vecs.push_back(&(tl_vecs.back()));
  return issued_vecs.back();
}

void PointCache::ReturnVec(std::vector<GeantParticle>* v) {
  issued_vecs.erase(std::find(issued_vecs.begin(), issued_vecs.end(), v), issued_vecs.end());
  return;
}

unsigned PointCache::NumIssued() const {
  return issued_vecs.size();
}

void PointCache::DumpData(std::vector<GeantParticle>* out_vec) const {
  // Copies the contents of all the vectors in tl_vecs to
  // the given vector
  size_t total_len = std::accumulate(tl_vecs.begin(), tl_vecs.end(), 0,
      [](const size_t tot, const auto& v) { return tot + v.size(); });
  out_vec->resize(total_len);
  for (const auto& v : tl_vecs)
    std::copy(v.begin(), v.end(), std::back_inserter(*out_vec));
  return;
}

void PointCache::Clear() {
  for (auto& v : tl_vecs) v.clear();
}



// Boilerplate for the iterator class
PointCacheIt::PointCacheIt(PointCache* cache) : p{nullptr}, c{cache}, c_it{nullptr} { }
PointCacheIt::PointCacheIt(std::vector<GeantParticle>::iterator pp, PointCache* pc,
    std::list<std::vector<GeantParticle>>::iterator pc_it) : p(pp), c(pc), c_it(pc_it) {}

PointCacheIt& PointCacheIt::operator++() {
  // If we've reached the end of one of the vectors
  // (other than the last one), shift to the next vector
  // Case: p != end -> increment p
  if (p != c_it->end()) ++p;
  // Case: p == end -> increment c_it until we get to
  // a non-empty vector, or the last vector, then set p to c_it->begin
  while (p == c_it->end() && c_it != std::prev(c->tl_vecs.end())) {
    ++c_it;
    p = c_it->begin();
  }
  return *this;
}

PointCacheIt PointCacheIt::operator++(int) {
  PointCacheIt tmp = *this;
  ++(*this);
  return tmp;
}

PointCacheIt& PointCacheIt::operator--() {
  // If we've reached the beginning of one of the vectors
  // (other than the first one), shift to the previous vector
  // Case: p != begin -> just decrement p
  if (p != c_it->begin()) {
    --p;
    return *this;
  }
  // Case: p == begin -> decrement c_it until we get to
  // a non-empty vector, or the last vector, then set p to c_it->end
  while (p == c_it->begin() && c_it != c->tl_vecs.begin()) {
    --c_it;
    p = c_it->end();
  }
  // Finally, decrement p so that it points to c_it->back
  --p;
  return *this;
}

PointCacheIt PointCacheIt::operator--(int) {
  PointCacheIt tmp = *this;
  --(*this);
  return tmp;
}

bool PointCacheIt::operator==(const PointCacheIt& o) const {
  return p == o.p && c == o.c && c_it == o.c_it;
}
bool PointCacheIt::operator!=(const PointCacheIt& o) const {
  return !this->operator==(o);
}

GeantParticle& PointCacheIt::operator*() { return *p; }
