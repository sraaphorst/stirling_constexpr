/**
 * stirling.h
 *
 * By Sebastian Raaphorst, 2018.
 */

#pragma once

#include <algorithm>
#include <array>

namespace stirling {
    /// Convenience aliases.
    template <size_t K>
    using matrix_row = std::array<size_t, K>;

    template <size_t N, size_t K>
    using matrix = std::array<matrix_row<K>, N>;

    namespace details {
        template <template <size_t, size_t> typename Function, size_t K, size_t Row, size_t... Cols>
        constexpr matrix_row<K>
        calc_matrix_row(std::index_sequence<Cols...>) {
            return {{Function<Row, Cols>::value()...}};
        }

        template <template <size_t, size_t> typename Function, size_t N, size_t K, size_t... Rows>
        constexpr matrix<N, K>
        calc_matrix_aux(std::index_sequence<Rows...>) {
            return {{calc_matrix_row<Function, K, Rows>(std::make_index_sequence<K>{})...}};
        }
    }

    /**
     * Create an N x K matrix (std::array<std::array<size_t, K>, N>) from the function that is passed in.
     * Note for Stirling numbers, it doesn't make sense for K > N, so only the lower triangular portion of the
     * matrix is of interest.
     *
     * @tparam Function
     * @tparam N
     * @tparam K
     * @return A constexpr evaluation of the function
     */
    template <template <size_t, size_t> typename Function, size_t N, size_t K>
    constexpr matrix<N, K> calc_matrix() {
        static_assert(N >= K, "N must be at least K.");
        return details::calc_matrix_aux<Function, N, K>(std::make_index_sequence<N>{});
    }

    /**
     * Generate Stirling numbers of the first kind:
     * https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind
     *
     * Often denoted s(N, K), it represents the number of permutations according to their number of cycles,
     * with fixed points counted as single cycles.
     *
     * @tparam N the number of objects to permute
     * @tparam K the number of cycles permitted
     */
    template <size_t N, size_t K>
    struct Stirling1 {
        static constexpr size_t value() {
            if constexpr(N == 0 && K == 0)
                return 1;
            else if constexpr(N == 0 || K == 0)
                return 0;
            else
                return (N - 1) * Stirling1<N - 1, K>::value() + Stirling1<N - 1, K - 1>::value();
        }
    };

    /**
     * Generate Stirling numbers of the second kind:
     * https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind
     *
     * Often denoted S(N, K), it represents the number of ways to partition N objects into K non-empty sets.
     * @tparam N
     * @tparam K
     */
    template<size_t N, size_t K>
    struct Stirling2 {
        static constexpr size_t value() {
            if constexpr(N == 0 && K == 0)
                return 1;
            else if constexpr(N == 0 || K == 0)
                return 0;
            else
                return K * Stirling2<N - 1, K>::value() + Stirling2<N - 1, K - 1>::value();
        }
    };
}