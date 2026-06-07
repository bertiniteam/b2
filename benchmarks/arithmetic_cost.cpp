/**
 * Benchmark to measure the relative arithmetic cost of mpfr_complex (MPC)
 * operations at various precisions, to calibrate ArithmeticCost() in
 * amp_tracker.hpp.
 *
 * Reference: Bates, Hauenstein, Sommese, Wampler, "Stepsize control for
 * adaptive multiprecision path tracking", 2009. Section 3.1 defines C(P):
 *
 *   "an approximation of C(P) was found using an average cost of computation
 *    in MPFR with different precisions compared with IEEE double precision.
 *    At various precisions, we computed the time of common operations used in
 *    homotopy continuation, e.g. straight-line program evaluation, matrix
 *    multiplication and linear solving."
 *
 *   C(P) = { 1,              if P = double precision
 *           { 10.35 + 0.04*P, otherwise    (P in BITS)
 *
 * Bertini 2 stores precision in decimal digits, so the code uses:
 *   0.04_bits * log2(10) ≈ 0.13_digits  →  10.35 + 0.13*P_digits
 * The intercept 10.35 is stale (2009 MPFR on Opteron 250).
 *
 * ------------------------------------------------------------------
 * KEY METHODOLOGICAL NOTES
 * ------------------------------------------------------------------
 * 1. Baseline is std::complex<double>, NOT double, because tracking
 *    is complex-valued.
 *
 * 2. Multiprecision complex uses GNU MPC via boost::multiprecision
 *    (same type as Bertini 2's mpfr_complex).
 *
 * 3. Simple scalar loops are NOT benchmarked — the compiler optimizes
 *    them away (-O2 can constant-fold a scalar multiply chain to 0 ns).
 *    Instead we benchmark composite operations like the paper:
 *
 *      (a) dot product over length-L vector: proxy for SLP evaluation
 *      (b) dense matrix-vector product (N x N)
 *      (c) LU factorization + solve (N x N)
 *
 * 4. Modern hardware vectorizes std::complex<double> (SSE/AVX) but
 *    not MPC. This CORRECTLY reflects the true cost ratio the tracker
 *    faces on modern hardware — i.e., C(P) should be large.
 *
 * Compile:
 *   micromamba run -n b2-ubuntu g++ -O2 -std=c++17 \
 *       -I/home/helper/micromamba/envs/b2-ubuntu/include \
 *       -L/home/helper/micromamba/envs/b2-ubuntu/lib \
 *       -Wl,-rpath,/home/helper/micromamba/envs/b2-ubuntu/lib \
 *       benchmarks/arithmetic_cost.cpp -o benchmarks/arithmetic_cost \
 *       -lmpfr -lgmp -lmpc
 *
 * Run:
 *   ./benchmarks/arithmetic_cost
 */

#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

namespace bmp = boost::multiprecision;
using mpfr_complex = bmp::number<bmp::backends::mpc_complex_backend<0>, bmp::et_off>;
using dbl_complex  = std::complex<double>;

using Clock   = std::chrono::high_resolution_clock;
using Seconds = std::chrono::duration<double>;

constexpr int MAT_N  = 8;    // matrix dimension for matvec and LU
constexpr int DOT_L  = 64;   // vector length for dot product (≈ SLP evaluation)

// ---------------------------------------------------------------------------
// Helper: ADL-friendly abs for pivot selection
// ---------------------------------------------------------------------------
template<typename C>
double mag(C const& x) {
    using std::abs;
    using boost::multiprecision::abs;
    return static_cast<double>(abs(x));
}

// ---------------------------------------------------------------------------
// (a) Dot product: s = sum_i (u[i] * v[i])  over a length-DOT_L vector.
//     Representative of SLP (straight-line program) evaluation.
//     The feed-back into u[0] prevents dead-store elimination.
// ---------------------------------------------------------------------------
double time_dot_dbl(size_t N_reps)
{
    std::vector<dbl_complex> u(DOT_L), v(DOT_L);
    for (int i = 0; i < DOT_L; ++i) {
        u[i] = {double(i % 7 + 1), double(i % 5 - 2)};
        v[i] = {double(i % 11 - 3), double(i % 9 + 1)};
    }

    dbl_complex s{};
    auto t0 = Clock::now();
    for (size_t r = 0; r < N_reps; ++r) {
        s = {};
        for (int i = 0; i < DOT_L; ++i)
            s += u[i] * v[i];
        u[0] = s + dbl_complex{1e-15, 0};  // loop-carried dep, prevents elim
    }
    auto t1 = Clock::now();
    (void)s;
    return Seconds(t1 - t0).count() / double(N_reps);
}

double time_dot_mp(unsigned digits, size_t N_reps)
{
    mpfr_complex::default_precision(digits);
    std::vector<mpfr_complex> u(DOT_L), v(DOT_L);
    for (int i = 0; i < DOT_L; ++i) {
        u[i] = mpfr_complex{double(i % 7 + 1), double(i % 5 - 2)};
        v[i] = mpfr_complex{double(i % 11 - 3), double(i % 9 + 1)};
    }
    mpfr_complex zero{0, 0};
    mpfr_complex tiny{1e-30, 0};

    mpfr_complex s{};
    auto t0 = Clock::now();
    for (size_t r = 0; r < N_reps; ++r) {
        s = zero;
        for (int i = 0; i < DOT_L; ++i)
            s += u[i] * v[i];
        u[0] = s + tiny;
    }
    auto t1 = Clock::now();
    (void)s;
    return Seconds(t1 - t0).count() / double(N_reps);
}

// ---------------------------------------------------------------------------
// (b) Dense matrix-vector product  y = A*x
// ---------------------------------------------------------------------------
double time_matvec_dbl(size_t N_reps)
{
    std::vector<dbl_complex> A(MAT_N * MAT_N), x(MAT_N), y(MAT_N);
    for (int i = 0; i < MAT_N * MAT_N; ++i)
        A[i] = {double(i % 7 + 1), double(i % 5 - 2)};
    for (int i = 0; i < MAT_N; ++i)
        x[i] = {double(i + 1), double(-i)};

    auto t0 = Clock::now();
    for (size_t r = 0; r < N_reps; ++r) {
        for (int i = 0; i < MAT_N; ++i) {
            y[i] = {};
            for (int j = 0; j < MAT_N; ++j)
                y[i] += A[i * MAT_N + j] * x[j];
        }
        x[0] = y[0] + dbl_complex{1e-15, 0};
    }
    auto t1 = Clock::now();
    (void)x;
    return Seconds(t1 - t0).count() / double(N_reps);
}

double time_matvec_mp(unsigned digits, size_t N_reps)
{
    mpfr_complex::default_precision(digits);
    std::vector<mpfr_complex> A(MAT_N * MAT_N), x(MAT_N), y(MAT_N);
    for (int i = 0; i < MAT_N * MAT_N; ++i)
        A[i] = mpfr_complex{double(i % 7 + 1), double(i % 5 - 2)};
    for (int i = 0; i < MAT_N; ++i)
        x[i] = mpfr_complex{double(i + 1), double(-i)};

    mpfr_complex zero{0, 0}, tiny{1e-30, 0};
    auto t0 = Clock::now();
    for (size_t r = 0; r < N_reps; ++r) {
        for (int i = 0; i < MAT_N; ++i) {
            y[i] = zero;
            for (int j = 0; j < MAT_N; ++j)
                y[i] += A[i * MAT_N + j] * x[j];
        }
        x[0] = y[0] + tiny;
    }
    auto t1 = Clock::now();
    (void)x;
    return Seconds(t1 - t0).count() / double(N_reps);
}

// ---------------------------------------------------------------------------
// (c) LU factorization + back-substitution (Gaussian elim, partial pivoting)
// ---------------------------------------------------------------------------
template<typename C>
void lu_solve(std::vector<C>& A, std::vector<C>& b, int n)
{
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double best = 0;
        for (int row = col; row < n; ++row) {
            double v = mag(A[row * n + col]);
            if (v > best) { best = v; pivot = row; }
        }
        if (pivot != col) {
            for (int j = 0; j < n; ++j)
                std::swap(A[col * n + j], A[pivot * n + j]);
            std::swap(b[col], b[pivot]);
        }
        C diag = A[col * n + col];
        for (int row = col + 1; row < n; ++row) {
            C factor = A[row * n + col] / diag;
            for (int j = col; j < n; ++j)
                A[row * n + j] -= factor * A[col * n + j];
            b[row] -= factor * b[col];
        }
    }
    for (int row = n - 1; row >= 0; --row) {
        for (int j = row + 1; j < n; ++j)
            b[row] -= A[row * n + j] * b[j];
        b[row] /= A[row * n + row];
    }
}

double time_lu_dbl(size_t N_reps)
{
    std::vector<dbl_complex> A0(MAT_N * MAT_N), b0(MAT_N);
    for (int i = 0; i < MAT_N * MAT_N; ++i)
        A0[i] = {double((i * 3 + 7) % 11 + 1), double(i % 5 - 2)};
    for (int i = 0; i < MAT_N; ++i)
        b0[i] = {double(i + 1), double(-i)};

    auto A = A0, b = b0;
    auto t0 = Clock::now();
    for (size_t r = 0; r < N_reps; ++r) {
        A = A0; b = b0;
        lu_solve(A, b, MAT_N);
        b0[0] = b[0] + dbl_complex{1e-15, 0};
    }
    auto t1 = Clock::now();
    (void)b;
    return Seconds(t1 - t0).count() / double(N_reps);
}

double time_lu_mp(unsigned digits, size_t N_reps)
{
    mpfr_complex::default_precision(digits);
    std::vector<mpfr_complex> A0(MAT_N * MAT_N), b0(MAT_N);
    for (int i = 0; i < MAT_N * MAT_N; ++i)
        A0[i] = mpfr_complex{double((i * 3 + 7) % 11 + 1), double(i % 5 - 2)};
    for (int i = 0; i < MAT_N; ++i)
        b0[i] = mpfr_complex{double(i + 1), double(-i)};

    mpfr_complex tiny{1e-30, 0};
    auto A = A0, b = b0;
    auto t0 = Clock::now();
    for (size_t r = 0; r < N_reps; ++r) {
        A = A0; b = b0;
        lu_solve(A, b, MAT_N);
        b0[0] = b[0] + tiny;
    }
    auto t1 = Clock::now();
    (void)b;
    return Seconds(t1 - t0).count() / double(N_reps);
}

// ---------------------------------------------------------------------------
// Linear fit  y = a + b*x
// ---------------------------------------------------------------------------
struct LinFit { double a, b; };
LinFit linear_fit(std::vector<double> const& x, std::vector<double> const& y)
{
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = (int)x.size();
    for (int i = 0; i < n; ++i) {
        sx += x[i]; sy += y[i]; sxx += x[i]*x[i]; sxy += x[i]*y[i];
    }
    double denom = n*sxx - sx*sx;
    return { (sy - (n*sxy-sx*sy)/denom * sx) / n,
             (n*sxy - sx*sy) / denom };
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main()
{
    std::vector<unsigned> prec_list;
    prec_list.push_back(16);  // DoublePrecision() in Bertini 2
    for (unsigned p = 20; p <= 300; p += 10)
        prec_list.push_back(p);

    const size_t N_dbl_dot = 1'000'000,  N_mp_dot = 100'000;
    const size_t N_dbl_mv  =   500'000,  N_mp_mv  =  50'000;
    const size_t N_dbl_lu  =   200'000,  N_mp_lu  =  20'000;

    std::cout << "Arithmetic cost benchmark: MPC complex vs std::complex<double>\n";
    std::cout << "DOT_L=" << DOT_L << " (dot/SLP),  MAT_N=" << MAT_N << " (matvec/LU)\n";
    std::cout << "P = decimal digits (Bertini 2 convention)\n";
    std::cout << "Paper formula used bits: C(P) = 10.35 + 0.04*P_bits\n";
    std::cout << "  => in decimal digits: 10.35 + 0.13*P_digits\n\n";

    // Warm up
    for (int w = 0; w < 3; ++w) {
        time_dot_dbl(N_dbl_dot); time_matvec_dbl(N_dbl_mv); time_lu_dbl(N_dbl_lu);
    }
    double base_dot = time_dot_dbl(N_dbl_dot);
    double base_mv  = time_matvec_dbl(N_dbl_mv);
    double base_lu  = time_lu_dbl(N_dbl_lu);

    std::cout << "Double baselines (ns/op):\n";
    std::cout << "  dot product (L=" << DOT_L << "): " << base_dot*1e9 << "\n";
    std::cout << "  matvec (N="  << MAT_N    << "):  " << base_mv*1e9  << "\n";
    std::cout << "  LU solve (N="<< MAT_N    << "): " << base_lu*1e9  << "\n\n";

    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(8)  << "digits"
              << std::setw(10) << "bits"
              << std::setw(12) << "dot_r"
              << std::setw(12) << "matvec_r"
              << std::setw(12) << "lu_r"
              << std::setw(12) << "mean_r"
              << "\n" << std::string(66, '-') << "\n";

    std::vector<double> fit_x, fit_y;

    for (unsigned p : prec_list) {
        double r_dot, r_mv, r_lu;
        if (p == 16) {
            r_dot = r_mv = r_lu = 1.0;
        } else {
            r_dot = time_dot_mp(p, N_mp_dot) / base_dot;
            r_mv  = time_matvec_mp(p, N_mp_mv)  / base_mv;
            r_lu  = time_lu_mp(p, N_mp_lu)       / base_lu;
        }
        double mean = (r_dot + r_mv + r_lu) / 3.0;
        unsigned bits = (p == 16) ? 53 : unsigned(std::round(p * std::log2(10.0)));

        std::cout << std::setw(8)  << p
                  << std::setw(10) << bits
                  << std::setw(12) << r_dot
                  << std::setw(12) << r_mv
                  << std::setw(12) << r_lu
                  << std::setw(12) << mean
                  << "\n";

        if (p > 16) { fit_x.push_back(double(p)); fit_y.push_back(mean); }
    }

    LinFit fit = linear_fit(fit_x, fit_y);

    std::cout << "\n--- Linear fit for P > 16 (P in decimal digits) ---\n";
    std::cout << std::setprecision(4);
    std::cout << "  C(P) = " << fit.a << " + " << fit.b << " * P\n";
    std::cout << "\nEquivalent in bits (paper comparison):\n";
    std::cout << "  a = " << fit.a
              << "  b_bits = " << fit.b / std::log2(10.0) << "\n";
    std::cout << "\nReplace in amp_tracker.hpp ArithmeticCost():\n";
    std::cout << "  return " << fit.a << " + " << fit.b << " * precision;\n";
    std::cout << "\nPaper (2009, Opteron 250, bits): 10.35 + 0.04 * P_bits\n";
    std::cout << "Current code (digits):           10.35 + 0.13 * P_digits\n";

    return 0;
}
