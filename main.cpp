#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include "util.hpp"

using sumCoord = util::sumCoord;
using Coord = util::Coord;
using Indexes = std::pair<uint32_t, uint32_t>;

static bool dumping;
static std::string DUMP_SPLINE;
static uint32_t MAX_VALUE;
static uint32_t DEVIATION;
static uint32_t MAX_N;

// For dumping
bool matrix[50][50];

std::vector<sumCoord> sums;
std::vector<Coord> cdf, spline;

unsigned getIndex(uint32_t sum) 
// extract the index of sum in cdf
{
    auto iter = std::lower_bound(cdf.begin(),
                                 cdf.end(),
                                 sum,
                                 [](auto temp, uint32_t x) -> bool {
                                    return temp.first < x;
                                 });
    assert(iter->first == sum);
    return iter - cdf.begin();
    
}


std::vector<Coord> buildBlocked(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take as sums only the sums (index, index)
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    assert(xs.size() == ys.size());
    
    std::unordered_map<uint32_t, bool> seen;
    
    unsigned limit = static_cast<unsigned>(sqrt(MAX_N));
    for (unsigned index = 0; index < MAX_N; ++index) {
        if (index % limit) continue;
        for (unsigned ptr = 0; ptr < MAX_N; ++ptr) {
            if (ptr % limit == 0) {
                uint32_t sum = xs[index] + ys[ptr];
                if (seen[sum]) continue;
                seen[sum] = true;
                ts.push_back(cdf[getIndex(xs[index] + ys[ptr])]);
            }
        }
    }
    std::sort(ts.begin(), ts.end(), [](auto left, auto right) {
       return left.first < right.first; 
    });
    return ts;
}


std::vector<Coord> buildArtificial(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take as sums only the sums (index, index)
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    assert(xs.size() == ys.size());
    uint32_t currPtr = 0;
    for (unsigned index = 0; index < MAX_N; ++index) {
        unsigned currSum = xs[index] + ys[index];
        while (currPtr != cdf.size() && cdf[currPtr].first < currSum)
            ++currPtr;
        ts.push_back(cdf[currPtr++]);
    }
    return ts;
}

std::vector<Coord> buildAnti(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take as sums only the sums (index, MAX_N - index - 1)
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    assert(xs.size() == ys.size());
    uint32_t currPtr = 0;
    
    
    std::unordered_map<uint32_t, bool> seen;
    
    for (unsigned index = 0; index < MAX_N; ++index) {
        unsigned currSum = xs[index] + ys[MAX_N - index - 1];
        if (seen[currSum]) continue;
        seen[currSum] = true;
        ts.push_back(cdf[getIndex(currSum)]);
    }
    ts.push_back(cdf.front());
    ts.push_back(cdf.back());
    std::sort(ts.begin(), ts.end(), [](auto left, auto right) {
       return left.first < right.first; 
    });
    return ts;
}

std::vector<Coord> buildWithCdf(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take as sums only the sums (index, index)
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    assert(xs.size() == ys.size());
    uint32_t currPtr = 0;
    
    
    auto inlineCdf = [&](const std::vector<uint32_t>& arr) {
        std::vector<Coord> tmp(arr.size());
        unsigned index = 0;
        for (auto v : arr)
            tmp[index++] = std::make_pair(v, index);
        return tmp;
    };

    std::vector<Coord> xs_cdf = inlineCdf(xs);
    std::vector<Coord> ys_cdf = inlineCdf(ys);
    
    unsigned desired = static_cast<unsigned>(sqrt(MAX_N));
    std::vector<Coord> xs_spline = util::compressFunc(xs_cdf, desired);
    std::vector<Coord> ys_spline = util::compressFunc(ys_cdf, desired);
    
    
    std::unordered_map<uint32_t, bool> seen;
    
    for (auto x_elem: xs_spline) {
        for (auto y_elem: ys_spline) {
            uint32_t sum = x_elem.first + y_elem.first;
            if (seen[sum]) continue;
            seen[sum] = true;
            ts.push_back(cdf[getIndex(x_elem.first + y_elem.first)]);
        }
    }
    std::sort(ts.begin(), ts.end(), [](auto left, auto right) {
       return left.first < right.first; 
    });
    return ts;
}

std::vector<Coord> buildWithMerge(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take as sums only the sums (index, index)
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    assert(xs.size() == ys.size());
    uint32_t currPtr = 0;
    
    std::vector<Coord> xs_cdf, ys_cdf;
    
    unsigned index = 0, ptr = 0, currPos = 0;
    while (index < xs.size() && ptr < ys.size()) {
        if (xs[index] < ys[ptr])
            xs_cdf.push_back(std::make_pair(xs[index++], currPos++));
        else
            ys_cdf.push_back(std::make_pair(ys[ptr++], currPos++));
    }
    while (index < xs.size())
        xs_cdf.push_back(std::make_pair(xs[index++], currPos++));
    while (ptr < ys.size())
        ys_cdf.push_back(std::make_pair(ys[ptr++], currPos++));
        
    unsigned desired = static_cast<unsigned>(sqrt(MAX_N));
    std::vector<Coord> xs_spline = util::compressFunc(xs_cdf, desired);
    std::vector<Coord> ys_spline = util::compressFunc(ys_cdf, desired);
    
    
    std::unordered_map<uint32_t, bool> seen;
    
    for (auto x_elem: xs_spline) {
        for (auto y_elem: ys_spline) {
            uint32_t sum = x_elem.first + y_elem.first;
            if (seen[sum]) continue;
            seen[sum] = true;
            ts.push_back(cdf[getIndex(x_elem.first + y_elem.first)]);
        }
    }
    std::sort(ts.begin(), ts.end(), [](auto left, auto right) {
       return left.first < right.first; 
    });
    return ts;
}

std::vector<Coord> buildNatural(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take the natural spline
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    return util::compressFunc(cdf, MAX_N);
}

std::vector<Coord> buildRandom(const char* msg, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Take random indexes and sum them up
{
    std::vector<Coord> ts;
    std::cerr << msg << std::endl;
    assert(xs.size() == ys.size());
    
    std::vector<uint32_t> firstOrder(MAX_N), secondOrder(MAX_N);
    for (unsigned index = 0; index < MAX_N; ++index)
        firstOrder[index] = secondOrder[index] = index;
    
    random_shuffle(firstOrder.begin(), firstOrder.end());
    random_shuffle(secondOrder.begin(), secondOrder.end());
    
    
    std::unordered_map<uint32_t, bool> seen;
    
    unsigned limit = static_cast<unsigned>(sqrt(MAX_N));
    for (unsigned index = 0; index < limit; ++index) {
        for (unsigned ptr = 0; ptr < limit; ++ptr) {
            uint32_t sum = xs[firstOrder[index]] + 
                ys[secondOrder[ptr]];
            if (seen[sum]) continue;
            seen[sum] = true;
            ts.push_back(cdf[getIndex(
                xs[firstOrder[index]] + 
                ys[secondOrder[ptr]]
            )]);
        }
    }
    std::sort(ts.begin(), ts.end(), [](auto left, auto right) {
       return left.first < right.first; 
    });
    return ts;
}

Indexes findSum(uint32_t sum, bool dump = false)
// find in sums the last occurence of sum and return the corresponding indexes
{
    auto iter = std::lower_bound(sums.begin(),
                                 sums.end(),
                                 sum,
                                 [](auto temp, uint32_t x) -> bool {
                                    return temp.first <= x;
                                 }) - 1;
    if (dump) {
        std::cerr << "For sum = " << sum << " get : " << iter->second.first << " & " << iter->second.second << std::endl;
    }
    return iter->second;
}

static std::vector<uint32_t> generator(bool deviate = false)
// generates a sorted array of MAX_N unsigned ints without duplicates
{
    std::vector<uint32_t> ret(MAX_N);
    std::unordered_map<uint32_t, bool> seen;
    for (unsigned index = 0; index < MAX_N; ++index) {
        uint32_t value;
        do {
           value = rand() % MAX_VALUE + 1 + deviate * DEVIATION; 
        } while (seen[value]);
        seen[value] = true;
        ret[index] = value;
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}

std::vector<Coord> choose(std::string name, std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Choose the type of the spline
{
    if (name == "artificial") return buildArtificial("artificial", xs, ys);
    if (name == "natural") return buildNatural("natural", xs, ys);
    if (name == "random") return buildRandom("random", xs, ys);
    if (name == "blocked") return buildBlocked("blocked", xs, ys);
    if (name == "anti") return buildBlocked("anti", xs, ys);
    if (name == "withCdf") return buildWithCdf("withCdf", xs, ys);
    if (name == "withMerge") return buildWithMerge("withMerge", xs, ys);
    return buildNatural("natural", xs, ys);
}

int main(int argc, char** argv) {
    if ((argc != 4) && (argc != 5)) {
        std::cerr << "Usage: " << argv[0] << " [DUMP: only with small MAX_N] [MAX_N] [MAX_VALUE] [optional: DEVIATION - default: 0]" << std::endl;
        return -1;
    }
    
    // Get the values
    DUMP_SPLINE = argv[1];
    dumping = DUMP_SPLINE.size() > 1;
    MAX_N = atoi(argv[2]);
    MAX_VALUE = atoi(argv[3]);
    DEVIATION = argc == 5 ? atoi(argv[4]) : 0;
    
    srand(time(NULL));
    std::vector<uint32_t> xs = generator(true);
    std::vector<uint32_t> ys = generator();
    
    auto printVector = [&](std::string msg, const std::vector<uint32_t>& arr) {
        std::cout << msg << std::endl;
        for (auto v : arr)
            std::cout << v << " ";
        std::cout << std::endl;
    };

#if 0
    printVector("X", xs);
    printVector("Y", ys);
#endif
    
    std::unordered_map<uint32_t, bool> seen;
    sums.reserve(xs.size() * ys.size());
    for (unsigned index = 0; index < xs.size(); ++index) {
        for (unsigned ptr = 0; ptr < ys.size(); ++ptr) {
            uint32_t sum = xs[index] + ys[ptr];
            if (seen[sum])
                continue;
            seen[sum] = true;
            sums.push_back(std::make_pair(xs[index] + ys[ptr], std::make_pair(index, ptr)));
        }
    }
    std::sort(sums.begin(), sums.end(), [](auto left, auto right) {
        return left.first < right.first;
    });
    
    // assert(sums.size() == MAX_N * MAX_N);
    cdf = util::buildCdf<sumCoord>(sums);
    std::cerr << "cdf size = " << cdf.size() << std::endl;
    
    unsigned desiredSize = 
#if 0
    MAX_N * log2(sqrt(MAX_N)); 
#else
    MAX_N;
#endif
    
    
    util::printErrors(cdf, buildArtificial("artificial", xs, ys));
    util::printErrors(cdf, buildNatural("natural", xs, ys));
    util::printErrors(cdf, buildRandom("random", xs, ys));
    util::printErrors(cdf, buildBlocked("blocked", xs, ys));
    util::printErrors(cdf, buildBlocked("anti", xs, ys));
    util::printErrors(cdf, buildWithCdf("withCdf", xs, ys));
    util::printErrors(cdf, buildWithMerge("withMerge", xs, ys));
    
    if (dumping) {
        assert(MAX_N <= 50);
        
        
        spline = choose(DUMP_SPLINE, xs, ys);
        
        util::printErrors(cdf, spline);
        std::cerr << "First sum = " << sums.front().first << " & last sum = " << sums.back().first << std::endl;  
        for (auto elem: spline) {
            uint32_t currSum = static_cast<uint32_t>(elem.first);
            Indexes ixs = findSum(currSum, false);
            
            //std::cerr << "Sum = " << elem.first << " and " << ixs.first << " & " << ixs.second << std::endl;
            matrix[ixs.first][ixs.second] = true;
        }
        for (unsigned index = 0; index < xs.size(); ++index) {
            for (unsigned ptr = 0; ptr < ys.size(); ++ptr) {
                std::cout << (matrix[index][ptr] ? '1' : '.') << " "; 
            }
            std::cout << std::endl;
        }
    }
    return 0;
}
