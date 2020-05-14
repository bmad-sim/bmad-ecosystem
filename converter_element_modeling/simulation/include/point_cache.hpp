#pragma once
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>

#include "bin.hpp"

class PointCache;

class PointCacheIt {
  private:
    std::vector<DataPoint>::iterator p;
    PointCache *c;
    std::vector<std::vector<DataPoint>>::iterator c_it;
    PointCacheIt(std::vector<DataPoint>::iterator, PointCache*,
        std::vector<std::vector<DataPoint>>::iterator);

  public:
    PointCacheIt() = delete;
    PointCacheIt(PointCache*);
    PointCacheIt(const PointCacheIt&) = default;
    PointCacheIt& operator=(const PointCacheIt&) = default;
    PointCacheIt(PointCacheIt&&) = default;
    PointCacheIt& operator=(PointCacheIt&&) = default;
    ~PointCacheIt() = default;

    PointCacheIt operator++();
    PointCacheIt operator++(int);
    PointCacheIt operator--();
    PointCacheIt operator--(int);
    bool operator==(const PointCacheIt& o) const;
    bool operator!=(const PointCacheIt& o) const;
    DataPoint& operator*();

    friend class PointCache;
};


class PointCache {
  private:
    std::vector<std::vector<DataPoint>> tl_vecs;
    std::vector<std::vector<DataPoint>*> issued_vecs;
    std::mutex vec_mutex;

    bool is_issued(std::vector<DataPoint>* v) const;

  public:
    PointCache();
    PointCache(const PointCache&) = delete;
    PointCache& operator=(const PointCache&) = delete;
    PointCache(PointCache&&) = default;
    PointCache& operator=(PointCache&&) = default;
    ~PointCache();

    void Lock();
    void Unlock();
    std::vector<DataPoint>* GetVec();
    void ReturnVec(std::vector<DataPoint>* v);

    unsigned NumIssued() const;

    void DumpData(std::vector<DataPoint>* out) const;
    void Clear();

    PointCacheIt begin() {
      return PointCacheIt(tl_vecs.begin()->begin(), this, tl_vecs.begin());
    }
    PointCacheIt end() {
      return PointCacheIt(tl_vecs.back().end(), this, tl_vecs.end()-1);
    }

    friend class PointCacheIt;
};


