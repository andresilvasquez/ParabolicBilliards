#ifndef BWM_BOUNDARY_WALL_METHOD_HPP
#define BWM_BOUNDARY_WALL_METHOD_HPP

#include <Eigen/Dense>
#include <boost/math/special_functions/hankel.hpp>
#include <complex>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <omp.h>

namespace bwm {

class BoundaryWallMethod {
public:
    using Complex   = std::complex<double>;
    using Vector2d  = Eigen::Vector2d;
    using VectorXcd = Eigen::VectorXcd;
    using MatrixXcd = Eigen::MatrixXcd;

    struct Spectrum {
        std::vector<double> k_values;
        std::vector<double> norm_values;
        std::vector<double> resonances;
    };

    BoundaryWallMethod(const std::vector<Vector2d>& boundary_points,
                       const std::vector<double>& gamma_vals,
                       double k_, double angle_, int threads_)
      : boundary(boundary_points),
        N(boundary_points.size()),
        k(k_), angle(angle_), threads(threads_)
    {
        if (gamma_vals.size() == 1)
            gamma.assign(N, gamma_vals[0]);
        else
            gamma = gamma_vals;

        kx = k * std::cos(angle);
        ky = k * std::sin(angle);

        segment_lengths = calculateSegmentLengths();
        M_flat          = buildMMatrixFlat();
        T               = buildTMatrix();
    }

    Complex incidentWave(double x, double y) const {
        static const Complex I(0.0, 1.0);
        return std::exp(I * (kx*x + ky*y));
    }

    Spectrum scanSpectrum(double k_min=0.1, double k_max=2.0,
                          int num_points=100,
                          double refine_window=0.05,
                          int refine_iterations=2,
                          double tolerance=1e-4)
    {
        auto [ks, norms] = initialScan(k_min, k_max, num_points);
        auto peaks       = findPeaks(ks, norms, k_min, k_max);
        auto refined     = refineResonances(peaks,
                                            refine_window,
                                            refine_iterations,
                                            tolerance,
                                            num_points/10);
        return { std::move(ks), std::move(norms), std::move(refined) };
    }

    std::vector<Complex> computeScatteredWave(const std::vector<Vector2d>& obs) {
        VectorXcd phi(N);
        for (size_t i = 0; i < N; ++i)
            phi[i] = incidentWave(boundary[i].x(), boundary[i].y());
        VectorXcd Tphi = T * phi;

        std::vector<Complex> psi(obs.size());
        #pragma omp parallel for num_threads(threads)
        for (size_t i = 0; i < obs.size(); ++i) {
            psi[i] = incidentWave(obs[i].x(), obs[i].y());
            for (size_t j = 0; j < N; ++j) {
                Complex G = green(obs[i], boundary[j]);
                psi[i] += G * segment_lengths[j] * Tphi[j];
            }
        }
        return psi;
    }

private:
    size_t idx(size_t i, size_t j) const { return i*N + j; }
    std::pair<size_t,size_t> ij(size_t id) const { return { id/N, id%N }; }

    std::vector<double> calculateSegmentLengths() const {
        std::vector<double> L(N);
        for (size_t i = 0; i < N; ++i) {
            size_t j = (i+1)%N;
            L[i] = (boundary[j] - boundary[i]).norm();
        }
        return L;
    }

    Complex green(const Vector2d& r1, const Vector2d& r2) const {
        double R = (r1 - r2).norm();                                  // Distance between the observation point and the boundary point
        if (R < 1e-10) return handleDiagonal();                       // If the observation point is near the boundary point, avoid infinite values
        static const Complex I(0.0, 1.0);                             // Declare the complex unit
        return I*(1.0/4.0) * boost::math::cyl_hankel_1(0, k*R); // Compute the Green function given the Hankel function of order zero
    }

    Complex handleDiagonal() const {
        double avg = std::accumulate(segment_lengths.begin(),
                                     segment_lengths.end(), 0.0) / N;
        static const Complex I(0.0, 1.0);
        return I*(1.0/4.0)
             * boost::math::cyl_hankel_1(0, k*(avg/2.0))
             * avg;
    }

    std::vector<Complex> buildMMatrixFlat() {
        std::vector<Complex> M(N*N);
        #pragma omp parallel for num_threads(threads)
        for (size_t id = 0; id < N*N; ++id) {
            auto [i,j] = ij(id);
            // If the entry Mij is diagonal call handleDiagonal(), otherwise compute the Green funtion times the segment lenght
            M[id] = (i==j
                     ? handleDiagonal()
                     : green(boundary[i], boundary[j]) * segment_lengths[j]);
        }
        return M;
    }

    MatrixXcd buildTMatrix() {
        MatrixXcd Mmat(N,N);
        for (size_t i=0;i<N;++i)
            for (size_t j=0;j<N;++j)
                Mmat(i,j) = M_flat[i*N+j];

        bool all_inf = std::all_of(gamma.begin(), gamma.end(),
                                   [](double g){ return std::isinf(g); });
        if (all_inf) {
            return -Mmat.inverse();
        } else {
            Eigen::VectorXd gvec = Eigen::Map<const Eigen::VectorXd>(gamma.data(), N);
            Eigen::MatrixXd G    = gvec.asDiagonal();
            MatrixXcd Iden = MatrixXcd::Identity(N,N);
            return G.cast<Complex>() * ((Iden - G.cast<Complex>()*Mmat).inverse());
        }
    }

    std::pair<std::vector<double>,std::vector<double>>
    initialScan(double k_min, double k_max, int np) {
        std::vector<double> ks(np), norms(np);
        for (int i=0;i<np;++i) {
            ks[i] = k_min + i*(k_max - k_min)/(np-1);
            updateK(ks[i]);
            norms[i] = T.array().abs().sum();
        }
        return {ks, norms};
    }

    std::vector<double> findPeaks(const std::vector<double>& ks,
                                  const std::vector<double>& norms,
                                  double k_min, double k_max) const
    {
        std::vector<bool> valid(norms.size(), true);
        std::vector<double> peaks;
        while (true) {
            double mv = -std::numeric_limits<double>::infinity();
            int midx = -1;
            for (size_t i=0;i<norms.size();++i)
                if (valid[i] && norms[i]>mv) {
                    mv = norms[i];
                    midx = i;
                }
            if (midx<0) break;
            valid[midx] = false;
            for (int i=midx+1; i+1<(int)norms.size(); ++i) {
                if (!valid[i]) continue;
                if (norms[i]>norms[i-1]) break;
                valid[i] = false;
            }
            for (int i=midx-1; i>0; --i) {
                if (!valid[i]) continue;
                if (norms[i]>norms[i+1]) break;
                valid[i] = false;
            }
            double kp = ks[midx];
            if (kp!=k_min && kp!=k_max)
                peaks.push_back(kp);
        }
        return peaks;
    }

    std::vector<double> refineResonances(const std::vector<double>& peaks,
                                         double window, int iter,
                                         double tol, int ppw)
    {
        std::vector<double> refined;
        for (double kg : peaks) {
            double ck = kg, w = window;
            for (int it=0; it<iter; ++it) {
                std::vector<double> kws(ppw), vals(ppw);
                for (int j=0;j<ppw;++j) {
                    kws[j] = ck - w + j*(2*w)/(ppw-1);
                    updateK(kws[j]);
                    vals[j] = T.array().abs().sum();
                }
                auto im = std::max_element(vals.begin(), vals.end());
                double nk = kws[std::distance(vals.begin(), im)];
                if (std::abs(nk - ck) < tol) break;
                ck = nk; w *= 0.1;
            }
            refined.push_back(ck);
        }
        return refined;
    }

    void updateK(double newk) {
        k  = newk;
        kx = k*std::cos(angle);
        ky = k*std::sin(angle);
        M_flat = buildMMatrixFlat();
        T      = buildTMatrix();
    }

    std::vector<Vector2d> boundary;
    size_t                N;
    int                   threads;
    std::vector<double>   gamma;
    double                k, angle, kx, ky;
    std::vector<double>   segment_lengths;
    std::vector<Complex>  M_flat;
    MatrixXcd             T;
};

} 

#endif 
