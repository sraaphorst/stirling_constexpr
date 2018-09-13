#include <iostream>

#include "stirling.h"
using namespace stirling;

int main() {
    constexpr size_t N = 10;
    constexpr size_t K = 10;

    std::cout << "Stirling numbers of the first kind:\n";
    constexpr matrix<N, K> S1s = calc_matrix<Stirling1, N, K>();
    for (size_t i = 0; i < N; ++i) {
        std::copy_n(S1s[i].cbegin(), i+1, std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << '\n';
    }

    std::cout << "\nStirling numbers of the second kind:\n";
    constexpr matrix<N, K> S2s = calc_matrix<Stirling2, N, K>();
    for (size_t i = 0; i < N; ++i) {
        std::copy_n(S2s[i].cbegin(), i+1, std::ostream_iterator<size_t>(std::cout, " "));
        std::cout << '\n';
    }
}