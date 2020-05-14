#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "bin.hpp"
#include "point_cache.hpp"

PointCacheIt::PointCacheIt(PointCache* cache) : p{nullptr}, c{cache}, c_it{nullptr} { }
PointCacheIt::PointCacheIt(std::vector<DataPoint>::iterator pp, PointCache* pc,
    std::vector<std::vector<DataPoint>>::iterator pc_it) : p(pp), c(pc), c_it(pc_it) {}

PointCacheIt PointCacheIt::operator++() {
  if (p == c_it->end()-1 && c_it != c->tl_vecs.end()-1) {
    ++c_it;
    p = c_it->begin();
    return *this;
  }
  ++p;
  return *this;
}

PointCacheIt PointCacheIt::operator++(int) {
  PointCacheIt tmp = *this;
  ++(*this);
  return tmp;
}

PointCacheIt PointCacheIt::operator--() {
  if (p == c_it->begin() && c_it != c->tl_vecs.begin()) {
    --c_it;
    p = c_it->end()-1;
    return *this;
  }
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

DataPoint& PointCacheIt::operator*() { return *p; }

PointCache::PointCache()  {
  tl_vecs.reserve(4);
}

PointCache::~PointCache() {
  // Acquire the mutex to make sure no one is still using
  std::scoped_lock<std::mutex> g(vec_mutex);
  // Check if all issued vectors have been returned
  if (issued_vecs.size())
    std::cerr << "WARNING: vector not returned to PointCache!\n";
}

bool PointCache::is_issued(std::vector<DataPoint>* v) const {
  return std::find(issued_vecs.begin(), issued_vecs.end(), v) != issued_vecs.end();
}

void PointCache::Lock() {
  vec_mutex.lock();
}
void PointCache::Unlock() { vec_mutex.unlock(); }

std::vector<DataPoint>* PointCache::GetVec() {
  //std::cout << "acquiring lock guard\n";
  //std::scoped_lock<std::mutex> g(vec_mutex);
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

void PointCache::ReturnVec(std::vector<DataPoint>* v) {
  //std::scoped_lock<std::mutex> g(vec_mutex);
  // Remove v from issued_vecs
  issued_vecs.erase(std::find(issued_vecs.begin(), issued_vecs.end(), v), issued_vecs.end());
  return;
}

unsigned PointCache::NumIssued() const {
  return issued_vecs.size();
}

void PointCache::DumpData(std::vector<DataPoint>* out_vec) const {
  // Copies the contents of all the vectors in tl_vecs to
  // the given vector
  //std::scoped_lock<std::mutex> g(vec_mutex);
  size_t total_len = std::accumulate(tl_vecs.begin(), tl_vecs.end(), 0,
      [](const size_t tot, const auto& v) { return tot + v.size(); });
  out_vec->resize(total_len);
  for (const auto& v : tl_vecs)
    std::copy(v.begin(), v.end(), std::back_inserter(*out_vec));
  return;
}

void PointCache::Clear() {
  //std::scoped_lock<std::mutex> g(vec_mutex);
  for (auto& v : tl_vecs) v.clear();
}
