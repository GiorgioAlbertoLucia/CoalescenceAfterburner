#ifndef JACOBITRANSFORM_H
#define JACOBITRANSFORM_H

#include "TVector3.h"

#include <vector>
#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// JacobiTransform
//
// Implements the orthonormal Jacobi transformation for A particles.
//
// The transformation matrix is defined as:
//
//   J_{0i} = 1/sqrt(A)                          for 0 <= i <= A-1  (CoM row)
//
//   J_{ni} = -1/sqrt(n(n+1))                    for 0 <= i <= n-1
//   J_{ni} =  n/sqrt(n(n+1))                    for i = n
//   J_{ni} =  0                                 otherwise
//
// for n = 1, ..., A-1.
//
// Given A lab-frame 3-vectors v_0,...,v_{A-1}, the Jacobi coordinates are:
//
//   xi_n = sum_{i=0}^{A-1} J_{ni} * v_i
//
// Row 0 gives the (normalised) centre of mass.
// Rows 1,...,A-1 give the A-1 relative Jacobi coordinates.
//
// Usage:
//   auto xi = JacobiTransform::transform(positions);  // all A coords
//   auto rel = JacobiTransform::relative(positions);  // rows 1..A-1 only
// ─────────────────────────────────────────────────────────────────────────────

class JacobiTransform {
public:

    // Returns the n-th row of the Jacobi matrix J_{ni} for a system of A particles.
    // n = 0 gives the CoM row; n = 1,...,A-1 give relative rows.
    static double element(int n, int i, int A) {
        if (n < 0 || n >= A || i < 0 || i >= A)
            throw std::out_of_range("JacobiTransform::element: index out of range");

        if (n == 0)
            return 1.0 / std::sqrt(static_cast<double>(A));

        if (i < n)
            return -1.0 / std::sqrt(static_cast<double>(n) * (n + 1));
        if (i == n)
            return static_cast<double>(n) / std::sqrt(static_cast<double>(n) * (n + 1));
        return 0.0;
    }

    // Transform A lab-frame 3-vectors to all A Jacobi coordinates.
    // Returns a vector of A TVector3s: index 0 is the (normalised) CoM,
    // indices 1..A-1 are the relative Jacobi coordinates.
    static std::vector<TVector3> transform(const std::vector<TVector3>& vecs) {
        const int A = static_cast<int>(vecs.size());
        if (A < 2)
            throw std::invalid_argument(
                "JacobiTransform::transform: need at least 2 vectors");

        std::vector<TVector3> result(A);
        for (int n = 0; n < A; ++n) {
            TVector3 xi(0., 0., 0.);
            for (int i = 0; i < A; ++i)
                xi += vecs[i] * element(n, i, A);
            result[n] = xi;
        }
        return result;
    }

    // Returns only the A-1 relative Jacobi coordinates (rows 1..A-1).
    // This is what the Wigner density depends on.
    static std::vector<TVector3> relative(const std::vector<TVector3>& vecs) {
        auto all = transform(vecs);
        return std::vector<TVector3>(all.begin() + 1, all.end());
    }

    // Convenience: extract TVector3s from a container of Particles or
    // any type with a .pos or .mom member — handled by the caller via lambdas.
    // See IntegratorA4 for usage.
};

#endif // JACOBITRANSFORM_H