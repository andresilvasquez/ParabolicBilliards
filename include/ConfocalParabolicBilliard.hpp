#ifndef BWM_CONFOCAL_PARABOLIC_BILLIARD_HPP
#define BWM_CONFOCAL_PARABOLIC_BILLIARD_HPP

#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace bwm {

class ConfocalParabolicBilliard {
public:
    using Vector2d = Eigen::Vector2d;

    ConfocalParabolicBilliard(double xi0_, double eta0_, int num_points_)
        : xi0(xi0_), eta0(eta0_), num_points(num_points_)
    {
        if (num_points % 4 != 0)
            throw std::invalid_argument("num_points must be a multiple of 4");
        boundary = generateBoundaryPoints();
    }

    const std::vector<Vector2d>& getBoundary() const { return boundary; }

private:
    double xi0, eta0;
    int    num_points;
    std::vector<Vector2d> boundary;

    static std::vector<double> linspace(double start, double end, int num) {
        std::vector<double> v(num);
        if (num == 1) {
            v[0] = start;
        } else {
            double step = (end - start) / (num - 1);
            for (int i = 0; i < num; ++i)
                v[i] = start + i * step;
        }
        return v;
    }

    std::vector<Vector2d> generateBoundaryPoints() const {
        int    n_half = num_points / 2;
        double sum    = xi0 + eta0;
        int    n_eta  = static_cast<int>(std::round(n_half * xi0  / sum));
        int    n_xi   = static_cast<int>(std::round(n_half * eta0 / sum));
        double step   = sum / n_half;

        // η = const (x<0)
        auto xi_seg1  = linspace(step/2.0, xi0 - step/2.0, n_eta);
        auto eta_seg1 = std::vector<double>(xi_seg1.size(), eta0);

        // ξ = const (x>0)
        auto eta_seg2 = linspace(eta0 - step/2.0, step/2.0, n_xi);
        auto xi_seg2  = std::vector<double>(eta_seg2.size(), xi0);

        double x_disp = (xi_seg2[0]*xi_seg2[0] - eta_seg2[0]*eta_seg2[0]) / 2.0;

        std::vector<double> xi_half, eta_half;
        xi_half.reserve(xi_seg1.size() + xi_seg2.size());
        eta_half.reserve(eta_seg1.size() + eta_seg2.size());
        xi_half.insert(xi_half.end(), xi_seg1.begin(), xi_seg1.end());
        xi_half.insert(xi_half.end(), xi_seg2.begin(), xi_seg2.end());
        eta_half.insert(eta_half.end(), eta_seg1.begin(), eta_seg1.end());
        eta_half.insert(eta_half.end(), eta_seg2.begin(), eta_seg2.end());

        std::vector<double> x_half(xi_half.size()), y_half(xi_half.size());
        for (size_t i = 0; i < xi_half.size(); ++i) {
            x_half[i] = ((xi_half[i]*xi_half[i] - eta_half[i]*eta_half[i]) / 2.0) - x_disp;
            y_half[i] = xi_half[i] * eta_half[i];
        }

        std::vector<double> x_ref = x_half, y_ref = y_half;
        std::reverse(x_ref.begin(), x_ref.end());
        std::reverse(y_ref.begin(), y_ref.end());
        for (auto& yy : y_ref) yy = -yy;

        std::vector<Vector2d> pts;
        pts.reserve(x_half.size()*2);
        for (size_t i = 0; i < x_half.size(); ++i)
            pts.emplace_back(x_half[i], y_half[i]);
        for (size_t i = 0; i < x_ref.size(); ++i)
            pts.emplace_back(x_ref[i], y_ref[i]);

        return pts;
    }
};

}

#endif
