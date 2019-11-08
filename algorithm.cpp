#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <chrono>
#include "util.hpp"
//-------------------------------------------------------------
#define SORTING 0
#define RADIX 1
#define CDF 2
#define TEST 3
//-------------------------------------------------------------
using namespace std::chrono;
using hash_coord = std::pair<uint64_t, bool>;
using cdf_coord = util::cdf_coord;
using sumCoord = util::sumCoord;
using Coord = util::Coord;
//-------------------------------------------------------------
static uint32_t OPTION;
static uint64_t MAX_VALUE;
static uint64_t DEVIATION;
static uint32_t MAX_N;
static uint64_t MAX_DIST_X;
static uint64_t MAX_DIST_Y;
//-------------------------------------------------------------
static std::vector<uint64_t> generator(uint64_t maxDist)
// generates a sorted array of MAX_N unsigned ints without duplicates
{
    std::vector<uint64_t> ret(MAX_N);
#if 0
    for (unsigned index = 0; index < MAX_N; ++index) {
        uint64_t value;
        do {
           value = rand() % MAX_VALUE + 1 + deviate * DEVIATION; 
        } while (seen[value]);
        seen[value] = true;
        ret[index] = value;
    }
    std::sort(ret.begin(), ret.end());
#else
    ret[0] = rand() % maxDist + 1;
    for (unsigned index = 1; index < MAX_N; ++index) {
        ret[index] = ret[index - 1] + rand() % maxDist + 1;
    }
#endif
    return ret;
}
//-------------------------------------------------------------
void applyFrequencies(std::unordered_map<uint64_t, hash_coord>& freqHashTable, std::vector<cdf_coord>& cdf)
// Update cdf with the frequencies from the hash table
{
    uint64_t partialSum = 0;
    for (unsigned index = 0, limit = cdf.size(); index < limit; ++index) {
        uint64_t currSum = freqHashTable[cdf[index].first].first;
        partialSum += currSum;
        cdf[index].second = partialSum;
    }
}
//-------------------------------------------------------------
std::unordered_map<uint64_t, hash_coord> computeFreqHashTable(std::vector<uint64_t>& xs, std::vector<uint64_t>& ys) 
// Return the hash table with the frequencies of the resulted sums
// This is used by all algorithms
{
    // Hash with the frequencies of each sum
    assert(xs.size() == ys.size());
    std::unordered_map<uint64_t, hash_coord> count_;
    for (unsigned index = 0; index < xs.size(); ++index) {
        for (unsigned ptr = 0; ptr < ys.size(); ++ptr) {
            // TODO: assuming that the value is still a 32-bit unsigned value
            // This serves only to the theoretical part of the algorithm
            uint64_t sum = xs[index] + ys[ptr];
            ++count_[sum].first;
        }
    }
    return count_;
}
//-------------------------------------------------------------
std::vector<cdf_coord> preprocessing_sorting(std::unordered_map<uint64_t, hash_coord>& count_, std::vector<uint64_t>& xs, std::vector<uint64_t>& ys)
// Solve the task with sorting
{
    std::vector<cdf_coord> final_(count_.size());
    unsigned index = 0;
    for (auto elem: count_)
        final_[index++] = std::make_pair(elem.first, 0);
    return final_;
}
//-------------------------------------------------------------
void solveWithSorting(std::vector<cdf_coord>& cdf)
{
    auto start = high_resolution_clock::now();
    std::sort(cdf.begin(), cdf.end(), [](auto left, auto right) {
       return left.first < right.first; 
    });
    auto stop = high_resolution_clock::now();
    std::cout << "Solving with Sorting took: " << duration_cast<milliseconds>(stop - start).count() << " ms" << std::endl;
}
//-------------------------------------------------------------
typedef struct {
    uint64_t sum, next;
} cell;
//-------------------------------------------------------------
std::vector<uint64_t> adj;
std::vector<cell> list;
//-------------------------------------------------------------
void addToChain(uint64_t ptr, uint64_t sum, uint64_t index)
// add in the hash-table the sum
{
    list[ptr].sum = sum;
    list[ptr].next = adj[index];
    adj[index] = ptr;
}
//-------------------------------------------------------------
std::vector<uint64_t> collectChain(uint64_t index) 
// Collect the chain into a vector
{
    std::vector<uint64_t> ret;
    for (uint32_t ptr = adj[index]; ptr; ptr = list[ptr].next)
        ret.push_back(list[ptr].sum);
    return ret;
}
//-------------------------------------------------------------
std::vector<cdf_coord> preprocessing_cdf(std::unordered_map<uint64_t, hash_coord>& count_, std::vector<uint64_t>& xs, std::vector<uint64_t>& ys)
// Solve with the cdf approximation
{
    // Construct an artificial spline
    assert(xs.size() == ys.size());
    std::vector<cdf_coord> spline(xs.size());
    for (unsigned index = 0; index < xs.size(); ++index)
        spline[index].first = xs[index] + ys[index];
    
    auto forward_pointer = [&spline](uint32_t& knotPtr, uint64_t sum) {
        // Iterate through the spline knots and find the segment of sum
        while ((knotPtr < spline.size()) && (spline[knotPtr].first <= sum))
            ++knotPtr;
        
        // Special case when reaching the last element
        knotPtr -= (knotPtr == spline.size());
    };
    
    // Hash with the frequencies of each sum
    unsigned sums_counter = count_.size();
    for (unsigned index = 0; index < xs.size(); ++index) {
        // Reset the pointer to the splien knot at every new line
        unsigned knotPtr = 0;
        for (unsigned ptr = 0; ptr < ys.size(); ++ptr) {
            // TODO: assuming that the value is still a 32-bit unsigned value
            // This serves only to the theoretical part of the algorithm
            uint64_t sum = xs[index] + ys[ptr];
            
            // Forward the pointer in the spline knots
            forward_pointer(knotPtr, sum);
            
            // First occurence of sum? Then update the "Mars Smen"
            if (count_[sum].second == false) {
                ++spline[knotPtr].second;
                count_[sum].second = true;
            }
        }
    }
    
    // Update the "Mars Smen"
    for (unsigned index = 1; index < spline.size(); ++index)
        spline[index].second += spline[index - 1].second;
    
#if 0
    std::cerr << "Final Spline : " << std::endl;
    for (auto elem: spline)
        std::cerr << elem.first << " " << elem.second << std::endl;
#endif
    
    auto interpolate = [&spline](uint32_t knotPtr, uint64_t sum) {
        // Compute the slope
        double dx = spline[knotPtr].first - spline[knotPtr - 1].first;
        double dy = spline[knotPtr].second - spline[knotPtr - 1].second;

        // Compute offset for the linear function
        double ofs = sum - spline[knotPtr - 1].first;
        return spline[knotPtr - 1].second + ofs * (dy / dx);
    };
    
#if 0
    std::cerr << sums_counter << " out of " << xs.size() * ys.size() << std::endl;
#endif
    
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
            uint64_t sum = xs[index] + ys[ptr];
            
            // Forward the pointer in the spline knots
            forward_pointer(knotPtr, sum);
            
            // First occurence of sum? (we already used it before, so it changes the meaning)
            if (count_[sum].second == true) {
                double estimation = interpolate(knotPtr, sum);
                uint32_t estimatedIndex = static_cast<uint32_t>(round(estimation));
                
                addToChain(++sums_counter, sum, estimatedIndex);
                
                // Mark the sum as seen
                count_[sum].second = false;
            }
        }
    }
    
    std::cerr << "Chain done" << std::endl;
    
    // Functions to add the values into the final order
    unsigned currBuffer = 0;
    std::vector<cdf_coord> final_(sums_counter);
    
    auto addKeyToFinal = [&final_, &currBuffer](uint64_t key) {
        final_[currBuffer++] = std::make_pair(key, 0);
    };
    
    auto addChain = [addKeyToFinal](std::vector<uint64_t>& values) {
        for (auto elem: values)
            addKeyToFinal(elem);
    };
    
    // Go through the chains from the hash-table
    std::vector<uint64_t> collector_;
    for (unsigned index = 0; index < sums_size; ++index) {
        // TODO: try also with sorting the chain
        collector_ = collectChain(index);
        // And insert each element of the chain
        addChain(collector_);
    }
    std::cerr << "Sums = " << sums_size << std::endl;
#if 0
    std::cerr << currBuffer << " vs " << sums_size << std::endl;
    assert(currBuffer == sums_size);
#endif
    
    std::cerr << "Final done" << std::endl;
    
#if 0
    util::printComputedCdf("final", final_);
#endif
    return final_;
}
//-------------------------------------------------------------    
void solveWithCdf(std::vector<cdf_coord>& cdf)
// Insertion sort for partially sorted cdf
{
    
    auto start = high_resolution_clock::now();
    int size = cdf.size();
    for (unsigned index = 1, limit = cdf.size(); index < limit; ++index) {
        unsigned ptr = index;
        cdf_coord curr = cdf[index];
        while ((ptr > 0) && (cdf[ptr - 1].first > curr.first)) {
            cdf[ptr] = cdf[ptr - 1];
            --ptr;
        }
        cdf[ptr] = curr;
    }
    
    auto stop = high_resolution_clock::now();
    std::cout << "Solving with CDF approximation took: " << duration_cast<milliseconds>(stop - start).count() << " ms" << std::endl;
}
//-------------------------------------------------------------    
std::vector<cdf_coord> chooseAlgorithm(uint32_t option)
// Benchmark the different algorithms, by returning the CDF of the sorting.
{
    // Final result - the CDF of the sorting
    std::vector<cdf_coord> ret;
    
    // Create the input
    std::vector<uint64_t> xs = generator(MAX_DIST_X);
    std::vector<uint64_t> ys = generator(MAX_DIST_Y);
    
    auto printVector = [&](std::string msg, const std::vector<uint64_t>& arr) {
        std::cout << msg << std::endl;
        for (auto v : arr)
            std::cout << v << " ";
        std::cout << std::endl;
    };

#if 0
    printVector("X", xs);
    printVector("Y", ys);
#endif
    
    std::unordered_map<uint64_t, hash_coord> count_ = computeFreqHashTable(xs, ys);
    
    // Analyze options
    switch (option) {
        // Sort with a comparison-based sorting algorithm (N^2 log N)
        // Sort with radix sort (N ^ 2)
        case RADIX : {
            assert(("Radix not yet supported", 0));
#if 0
            cdf = solveWithRadix(xs, ys);
#endif
            break;
        }
        case TEST : {
            std::cerr << "Testing" << std::endl;
        }
        case SORTING : {
#if 1
            auto start = high_resolution_clock::now();
            ret = preprocessing_sorting(count_, xs, ys);
            solveWithSorting(ret);
            applyFrequencies(count_, ret);
            auto stop = high_resolution_clock::now();
            std::cout << "Sorting (preprocessing + solve) took: " << duration_cast<milliseconds>(stop - start).count() << " ms" << std::endl;
#endif
            if (option != TEST)
                break;
        }
        // Sort with the CDF approximation (N ^ 2)
        case CDF : {
            auto start = high_resolution_clock::now();
            ret = preprocessing_cdf(count_, xs, ys);
            solveWithCdf(ret);
            applyFrequencies(count_, ret);
            auto stop = high_resolution_clock::now();
            std::cout << "CDF approximation (preprocessing + solve) took: " << duration_cast<milliseconds>(stop - start).count() << " ms" << std::endl;
            break;
        }
        default : {
            std::cerr << "Option " << option << " not (yet) available" << std::endl;
        }
    }
#if 0
    util::printComputedCdf("Choose algorithm", cdf);
#endif
    std::cerr << "Is cdf? " << util::is_cdf(xs, ys, ret) << std::endl;
    return ret;
}
//-------------------------------------------------------------
int main(int argc, char** argv) {
    srand(time(NULL));
    
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " [OPTION] [MAX_N] [MAX_DIST_X] [MAX_DIST_Y]" << std::endl;
        return -1;
    }
    
    OPTION = atoi(argv[1]);
    MAX_N = atoi(argv[2]);
    MAX_DIST_X = atoi(argv[3]);
    MAX_DIST_Y = atoi(argv[4]);

    std::vector<cdf_coord> ret = chooseAlgorithm(OPTION);
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
