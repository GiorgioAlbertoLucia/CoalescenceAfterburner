#ifndef JACOBITRANSFORM_H
#define JACOBITRANSFORM_H

#include "TVector3.h"
#include "TMatrixD.h"

#include <vector>
#include <cmath>
#include <stdexcept>

/**
 * Implements the orthonormal Jacobi transformation for A particles.
 *
 * The transformation matrix is defined as:
 *
 *   J_{0i} = 1/sqrt(A)                          for 0 <= i <= A-1  (CoM row)
 *   J_{ni} = -1/sqrt(n(n+1))                    for 0 <= i <= n-1
 *   J_{ni} =  n/sqrt(n(n+1))                    for i = n
 *   J_{ni} =  0                                 otherwise
 *
 * for n = 1, ..., A-1.
 *
 * Given A lab-frame 3-vectors v_0,...,v_{A-1}, the Jacobi coordinates are:
 *
 *   xi_n = sum_{i=0}^{A-1} J_{ni} * v_i
 *
 * Row 0 gives the (normalised) centre of mass.
 * Rows 1,...,A-1 give the A-1 relative Jacobi coordinates.
 *
 * Usage:
 *   JacobiTransform jt(A);              // build once, A fixed for the object's lifetime
 *   auto xi  = jt.transform(positions); // all A coords
 *   auto rel = jt.relative(positions);  // rows 1..A-1 only
 */
class JacobiTransform {
public:

    /**
     * @brief Construct the transform for a fixed number of particles A.
     *        The A x A Jacobi matrix is built once here and reused by
     *        all subsequent calls to transform()/relative().
     * @param A Total number of particles (must be >= 2)
     * @throws std::invalid_argument if A < 2
     */
    explicit JacobiTransform(int A) : fA(A), fJ(A, A) {
        if (A < 2)
            throw std::invalid_argument(
                "JacobiTransform::JacobiTransform: need at least 2 particles");

        for (int n = 0; n < A; ++n)
            for (int i = 0; i < A; ++i)
                fJ(n, i) = element(n, i, A);
    }

    /**
     * @brief Get the element of the Jacobi transformation matrix
     * @param n Row index (0 for CoM, 1..A-1 for relative coordinates)
     * @param i Column index (0..A-1)
     * @param A Total number of particles
     * @return The value of the Jacobi transformation matrix at (n, i)
     * @throws std::out_of_range if n or i are out of bounds
     */
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

    /**
     * @brief Number of particles this instance was built for.
     */
    int A() const { return fA; }

    /**
     * @brief Transform A lab-frame 3-vectors to all A Jacobi coordinates
     * @param vecs A vector of A TVector3s (lab-frame 3-vectors)
     * @return A vector of A TVector3s: index 0 is the (normalised) CoM, indices 1..A-1 are the relative Jacobi coordinates
     * @throws std::invalid_argument if vecs.size() != A
     */
    std::vector<TVector3> transform(const std::vector<TVector3>& vecs) const {
        if (static_cast<int>(vecs.size()) != fA)
            throw std::invalid_argument(
                "JacobiTransform::transform: input size does not match A");

        TMatrixD V(fA, 3);
        for (int i = 0; i < fA; ++i) {
            V(i, 0) = vecs[i].X();
            V(i, 1) = vecs[i].Y();
            V(i, 2) = vecs[i].Z();
        }

        TMatrixD Xi(fA, 3);
        Xi.Mult(fJ, V);

        std::vector<TVector3> result(fA);
        for (int n = 0; n < fA; ++n)
            result[n].SetXYZ(Xi(n, 0), Xi(n, 1), Xi(n, 2));

        return result;
    }

    /**
     * @brief Transform A Particles (positions and momenta) to Jacobi coordinates
     * @param particles A vector of A Particles
     * @return A vector of A Particles with Jacobi-transformed positions and momenta
     * @throws std::invalid_argument if particles.size() != A
     */
    std::vector<Particle> transform(const std::vector<Particle>& particles) const {
        if (static_cast<int>(particles.size()) != fA)
            throw std::invalid_argument(
                "JacobiTransform::transform: input size does not match A");

        std::vector<TVector3> jacobiPositions(fA), jacobiMomenta(fA);
        for (int n = 0; n < fA; ++n) {
            jacobiPositions[n] = particles[n].pos;
            jacobiMomenta[n] = particles[n].mom.Vect();
        }

        auto r_jacobi = transform(jacobiPositions);
        auto k_jacobi = transform(jacobiMomenta);

        std::vector<Particle> result;
        for (int n = 0; n < fA; ++n) {
            TLorentzVector mom(k_jacobi[n], 0.);
            const float mass = particles[n].mom.M();
            mom.SetE(std::sqrt(k_jacobi[n].Mag2() + mass * mass));
            result.emplace_back(particles[n].pdg, mom, r_jacobi[n]);
        }

        return result;
    }

private:
    int fA;
    TMatrixD fJ;   // cached A x A Jacobi matrix, built once in the constructor
};

#endif // JACOBITRANSFORM_H