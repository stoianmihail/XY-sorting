#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include "util.hpp"
//-------------------------------------------------------------
#define SORTING 0
#define RADIX 1
#define CDF 2
#define TEST 3
//-------------------------------------------------------------
using cdf_coord = util::cdf_coord;
using sumCoord = util::sumCoord;
using Coord = util::Coord;
//-------------------------------------------------------------
static uint32_t OPTION;
static uint32_t MAX_VALUE;
static uint32_t DEVIATION;
static uint32_t MAX_N;
//-------------------------------------------------------------
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
//-------------------------------------------------------------
typedef struct {
    uint32_t sum, next;
} cell;
//-------------------------------------------------------------
std::vector<uint32_t> adj;
std::vector<cell> list;
//-------------------------------------------------------------
void addToChain(uint32_t ptr, uint32_t sum, uint32_t index)
// add in the hash-table the sum
{
    list[ptr].sum = sum;
    list[ptr].next = adj[index];
    adj[index] = ptr;
}
//-------------------------------------------------------------
std::vector<uint32_t> collectChain(uint32_t index) 
// Collect the chain into a vector
{
    std::vector<uint32_t> ret;
    for (uint32_t ptr = adj[index]; ptr; ptr = list[ptr].next)
        ret.push_back(list[ptr].sum);
    return ret;
}
//-------------------------------------------------------------
std::vector<cdf_coord> solveWithCdf(std::vector<uint32_t>& xs, std::vector<uint32_t>& ys)
// Solve with the cdf approximation
{
    // Construct an artificial spline
    assert(xs.size() == ys.size());
    std::vector<cdf_coord> spline(xs.size());
    for (unsigned index = 0; index < xs.size(); ++index)
        spline[index].first = xs[index] + ys[index];
    
    auto forward_pointer = [&spline](uint32_t& knotPtr, uint32_t sum) {
        // Iterate through the spline knots and find the segment of sum
        while ((knotPtr < spline.size()) && (spline[knotPtr].first <= sum))
            ++knotPtr;
        
        // Special case when reaching the last element
        knotPtr -= (knotPtr == spline.size());
    };
    
    // Hash with the frequencies of each sum
    std::unordered_map<uint32_t, uint32_t> count_;
    unsigned sums_counter = 0;
    for (unsigned index = 0; index < xs.size(); ++index) {
        // Reset the pointer to the splien knot at every new line
        unsigned knotPtr = 0;
        for (unsigned ptr = 0; ptr < ys.size(); ++ptr) {
            // TODO: assuming that the value is still a 32-bit unsigned value
            // This serves only to the theoretical part of the algorithm
            uint32_t sum = xs[index] + ys[ptr];
            
            // Forward the pointer in the spline knots
            forward_pointer(knotPtr, sum);
            
            // First occurence of sum? Then update the "Mars Smen"
            if (!count_[sum]) {
                ++spline[knotPtr].second;
                ++sums_counter;
            }
            
            // Increase the frequency of the sum, but keep the first bit free
            count_[sum] += 2;
        }
    }
    
    // Update the "Mars Smen"
    for (unsigned index = 1; index < spline.size(); ++index)
        spline[index].second += spline[index - 1].second;
    assert(spline.back().second == sums_counter);
    
#if 0
    std::cerr << "Final Spline : " << std::endl;
    for (auto elem: spline)
        std::cerr << elem.first << " " << elem.second << std::endl;
#endif
    
    auto interpolate = [&spline](uint32_t knotPtr, uint32_t sum) {
        // Compute the slope
        double dx = spline[knotPtr].first - spline[knotPtr - 1].first;
        double dy = spline[knotPtr].second - spline[knotPtr - 1].second;

        // Compute offset for the linear function
        double ofs = sum - spline[knotPtr - 1].first;
        return spline[knotPtr - 1].second + ofs * (dy / dx);
    };
    
    std::cerr << sums_counter << " out of " << xs.size() * ys.size() << std::endl;
    
    // Compute the CDF approximation
    // Create the hash table
    uint32_t sums_size = sums_counter + 1;
    adj.resize(sums_size);
    list.resize(sums_size + 1);
    
    sums_counter = 0;
    for (unsigned index = 0; index < xs.size(); ++index) {
        unsigned knotPtr = 0;
        for (unsigned ptr = 0; ptr < ys.size(); ++ptr) {
            // TODO: assuming that the value is stil a 32-bit unsigned value
            // This serves only to the theoretical part of the algorithm
            uint32_t sum = xs[index] + ys[ptr];
            
            // Forward the pointer in the spline knots
            forward_pointer(knotPtr, sum);
            
            // First occurence of sum?
            if ((count_[sum] & 1) == 0) {
                double estimation = interpolate(knotPtr, sum);
                uint32_t estimatedIndex = static_cast<uint32_t>(round(estimation));
                
                assert(0 <= estimatedIndex && estimatedIndex < sums_size);
                addToChain(++sums_counter, sum, estimatedIndex);
                
                // Mark the sum as seen
                count_[sum] |= 1;
            }
        }
    }
    
    std::cerr << "Chain done" << std::endl;
    
    // Functions to add the values into the final order
    unsigned currBuffer = 0;
    std::vector<cdf_coord> final_(sums_counter);
    
    auto addKeyToFinal = [&final_, &currBuffer](uint32_t key) {
        cdf_coord tmpPair = std::make_pair(key, 0);
        
        // Find the position where to insert key
        uint32_t curr = currBuffer++;
        while ((curr > 0) && (final_[curr - 1].first > tmpPair.first)) {
            final_[curr] = final_[curr - 1];
            --curr;
        }
        final_[curr] = tmpPair;
    };
    
    auto addChain = [addKeyToFinal](std::vector<uint32_t>& values) {
        for (auto elem: values)
            addKeyToFinal(elem);
    };
    
    // Go through the chains from the hash-table
    std::vector<uint32_t> collector_;
    for (unsigned index = 0; index < sums_size; ++index) {
        // TODO: try also with sorting the chain
        collector_ = collectChain(index);
        // And insert each element of the chain
        addChain(collector_);
    }
#if 0
    std::cerr << currBuffer << " vs " << sums_size << std::endl;
    assert(currBuffer == sums_size);
#endif
    
    std::cerr << "Final done" << std::endl;
    
#if 0
    util::printComputedCdf("final", final_);
#endif
    
    // Determine the cdf by looking up the hash table
    unsigned partialSum = 0;
    for (unsigned index = 0; index < final_.size(); ++index) {
        uint32_t currSum = count_[final_[index].first] >> 1;
        final_[index].second += currSum;
        partialSum += currSum;
    }
    assert(partialSum == xs.size() * ys.size());
    return final_;
}
//-------------------------------------------------------------
std::vector<cdf_coord> chooseAlgorithm(uint32_t option)
// Benchmark the different algorithms, by returning the CDF of the sorting.
{
    // Final result - the CDF of the sorting
    std::vector<cdf_coord> cdf;
    
    // Create the input
    std::vector<uint32_t> xs = generator(true);
    std::vector<uint32_t> ys = generator();
    
    auto printVector = [&](std::string msg, const std::vector<uint32_t>& arr) {
        std::cout << msg << std::endl;
        for (auto v : arr)
            std::cout << v << " ";
        std::cout << std::endl;
    };

    // Analyze options
    switch (option) {
        // Sort with a comparison-based sorting algorithm (N^2 log N)
        case SORTING : {
#if 0
            cdf = solveWithSorting(xs, ys);
#endif
            break;
        }
        // Sort with radix sort (N ^ 2)
        case RADIX : {
#if 0
            cdf = solveWithRadix(xs, ys);
#endif
            break;
        }
        // Sort with the CDF approximation (N ^ 2)
        case CDF : {
            cdf = solveWithCdf(xs, ys);
            break;
        }
        // Sort with all the above options
        case TEST : {
#if 0
            cdf = solveWithSorting(xs, ys);
            cdf = solveWithRadix(xs, ys);
            cdf = solveWithCdf(xs, ys);
#endif
        }
        default : {
            std::cerr << "Option " << option << " not (yet) available" << std::endl;
        }
    }
#if 0
    util::printComputedCdf("Choose algorithm", cdf);
#endif
    return cdf;
}
//-------------------------------------------------------------
int main(int argc, char** argv) {
    srand(time(NULL));
    
    if ((argc != 4) && (argc != 5)) {
        std::cerr << "Usage: " << argv[0] << " [OPTION] [MAX_N] [MAX_VALUE] [optional: DEVIATION - default: 0]" << std::endl;
        return -1;
    }
    
    OPTION = atoi(argv[1]);
    MAX_N = atoi(argv[2]);
    MAX_VALUE = atoi(argv[3]);
    DEVIATION = (argc == 5) ? atoi(argv[4]) : 0;

    std::vector<cdf_coord> ret = chooseAlgorithm(OPTION);
    std::cerr << "Is sorted? " << util::is_sorted(ret) << std::endl;

#if 0
    printVector("X", xs);
    printVector("Y", ys);
#endif

#if 0
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

#endif

#if 0    
    cdf = util::buildCdf<sumCoord>(sums);
    spline = buildArtificial("artificial", xs, ys);
    
    util::printErrors(cdf, buildArtificial("artificial", xs, ys));
    
    std::cerr << "First point of spline: (" << spline.front().first << "," << spline.front().second << ") and last point of spline: (" << spline.back().first << "," << spline.back().second << ")" << std::endl; 
#endif
    
#if 0
    for (auto elem: spline) {
        std::cerr << elem.first << " " << elem.second << std::endl;
    }
#endif
    
#if 0
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
#endif
    return 0;
}
