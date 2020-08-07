#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <thread>
#include <mutex>

#include "bin.hpp"

class PointCache;

// Define an iterator for the point cache to enable range-for
// and allow use with STL algorithms
class PointCacheIt {
  private:
    std::vector<GeantParticle>::iterator p;
    PointCache *c;
    std::list<std::vector<GeantParticle>>::iterator c_it;
    PointCacheIt(std::vector<GeantParticle>::iterator, PointCache*,
        std::list<std::vector<GeantParticle>>::iterator);

  public:
    PointCacheIt() = delete;
    PointCacheIt(PointCache*);
    PointCacheIt(const PointCacheIt&) = default;
    PointCacheIt& operator=(const PointCacheIt&) = default;
    PointCacheIt(PointCacheIt&&) = default;
    PointCacheIt& operator=(PointCacheIt&&) = default;
    ~PointCacheIt() = default;

    PointCacheIt& operator++();
    PointCacheIt operator++(int);
    PointCacheIt& operator--();
    PointCacheIt operator--(int);
    bool operator==(const PointCacheIt& o) const;
    bool operator!=(const PointCacheIt& o) const;
    GeantParticle& operator*();

    friend class PointCache;
};


class PointCache {
  // This class implements a data structure that is suitable for use
  // with Geant's approach to multithreading. It allows Geant's
  // worker threads to obtain a std::vector into which they can store
  // their simulation results.  This class provides iterators and a
  // DumpData method that the single threaded binning code can use to
  // obtain Geant's results.

  // WHENEVER YOU TOUCH AN INSTANCE OF THIS CLASS, YOU MUST CALL
  // THE METHOD Lock(). WHEN YOU ARE DONE TOUCHING THE INSTANCE,
  // YOU MUST CALL Unlock().

  // An alterantive aproach would be to have each member function lock the
  // mutex before it exectutes. The user would then not be responsible for
  // manually locking and unlocking. However, this approach has a couple of
  // disadvantages. For one, the Lock() and Unlock() methods would still be
  // necessary to lock the mutex while the iterators are in use. Additionally,
  // the use of separate Lock/Unlock methods allow you to call multiple member
  // functions while you have the mutex locked. If the methods were self-locking,
  // one thread might lock the mutex in between two member function calls
  // of another thread.
  private:
    // tl_vecs must be implemented as a list instead of a vector;
    // This class relies on the vectors it manages to not change their location
    // in memory; std::vector breaks this requirement when it grows.
    // This is why tl_vecs must be a list of vectors rather than a vector
    // of vectors.
    std::list<std::vector<GeantParticle>> tl_vecs; // thread-local vectors
    std::vector<std::vector<GeantParticle>*> issued_vecs;
    std::mutex vec_mutex;

    bool is_issued(std::vector<GeantParticle>* v) const;

  public:
    PointCache() = default;
    PointCache(const PointCache&) = delete;
    PointCache& operator=(const PointCache&) = delete;
    PointCache(PointCache&&) = default;
    PointCache& operator=(PointCache&&) = default;
    ~PointCache();

    void Lock();
    void Unlock();
    std::vector<GeantParticle>* GetVec();
    void ReturnVec(std::vector<GeantParticle>* v);

    unsigned NumIssued() const;

    void DumpData(std::vector<GeantParticle>* out) const;
    void Clear();

    inline PointCacheIt begin() {
      auto first_nonempty_vec = std::find_if(tl_vecs.begin(), tl_vecs.end(), [](const auto& v) { return v.size() != 0; });
      if (first_nonempty_vec == tl_vecs.end()) return end();
      return PointCacheIt(first_nonempty_vec->begin(), this, first_nonempty_vec);
    }

    inline PointCacheIt end() {
      return PointCacheIt(tl_vecs.back().end(), this, std::prev(tl_vecs.end()));
    }

    friend class PointCacheIt;
};


