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
std::vector<uint32_t> final_;

unsigned listSize;

typedef struct {
    uint32_t sum, next;
} cell;

std::vector<uint32_t> adj;
std::vector<cell> list;

void addSum(uint32_t sum, uint32_t index)
// add in the hash-table the sum
{
    if (!(0 <= index && index < sums.size())) {
        std::cerr << index << std::endl;
        assert(0);
    }
    ++listSize;
    list[listSize].sum = sum;
    list[listSize].next = adj[index];
    adj[index] = listSize;
}

std::vector<uint32_t> computeChain(uint32_t index)
// get the length of the chain "index"
{
    assert(0 <= index && index < sums.size());
    std::vector<uint32_t> ret;
    for (uint32_t ptr = adj[index]; ptr; ptr = list[ptr].next)
        ret.push_back(list[ptr].sum);
    return ret;
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

  
double interpolate(uint32_t pos)
// evaluate the linear function found on spline for "pos"
{
    auto iter = lower_bound(spline.begin(), spline.end(), pos, [](const Coord &a, double b) { return a.first<b; });
    if (iter->first == pos)
        return iter->second;

    // Compute slope
    double dx = (iter + 0)->first - (iter - 1)->first;
    double dy = (iter + 0)->second - (iter - 1)->second;

    // Compute offset for the linear function
    double ofs = pos - (iter - 1)->first;
    return (iter - 1)->second + ofs * (dy / dx);
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

void insert(uint32_t pos, uint32_t key)
// insert the key in the final order
{
    uint32_t curr = pos;
    while ((curr > 0) && (final_[curr - 1] > key)) {
        final_[curr] = final_[curr - 1];
        --curr;
    }
    final_[curr] = key;
}

int main(int argc, char** argv) {
    srand(time(NULL));
    
    if ((argc != 3) && (argc != 4)) {
        std::cerr << "Usage: " << argv[0] << " [MAX_N] [MAX_VALUE] [optional: DEVIATION - default: 0]" << std::endl;
        return -1;
    }
    
    MAX_N = atoi(argv[1]);
    MAX_VALUE = atoi(argv[2]);
    DEVIATION = argc == 4 ? atoi(argv[3]) : 0;
    
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
    
#if 1
    std::sort(sums.begin(), sums.end(), [](auto left, auto right) {
        return left.first < right.first;
    });
#endif
    
    cdf = util::buildCdf<sumCoord>(sums);
    spline = buildArtificial("artificial", xs, ys);
    
    util::printErrors(cdf, buildArtificial("artificial", xs, ys));
    
    std::cerr << "First point of spline: (" << spline.front().first << "," << spline.front().second << ") and last point of spline: (" << spline.back().first << "," << spline.back().second << ")" << std::endl; 
    
#if 0
    for (auto elem: spline) {
        std::cerr << elem.first << " " << elem.second << std::endl;
    }
#endif
    
    adj.resize(sums.size(), 0);
    list.resize(sums.size() + 1); // the first buffer is not used
    for (auto elem: sums) {
        uint32_t currSum = elem.first;
        double estimation = interpolate(currSum);
        uint32_t index = static_cast<uint32_t>(round(estimation));
        index = index % sums.size();
        addSum(currSum, index);
    }
    
    std::cerr << "alles gut hier" << std::endl;
    
    unsigned currPos = 0;
    double avgLen = 0, maxLen = 0;
    std::cerr << "before resize with " << sums.size() << std::endl;
    final_.resize(sums.size());
    std::cerr << "after resize" << std::endl;
    for (unsigned index = 0; index < sums.size(); ++index) {
        // std::cerr << "Progress: " << index << std::endl;
        std::vector<uint32_t> currChain = computeChain(index);
        if (!currChain.empty()) {
            for (auto elem: currChain)
                insert(currPos++, elem);
        }
    }
    
    assert(currPos == sums.size());
    
    std::cerr << "Is sorted? " << std::is_sorted(final_.begin(), final_.end()) << std::endl;
    return 0;
}
