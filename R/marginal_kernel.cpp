// marginal_kernel.cpp
// -----------------------------------------------------------------------
// Unit log marginal likelihoods for the latent capability model, with the
// capability integrated out by Gauss-Legendre quadrature. Same algorithm as
// unit_loglik_marg() in 30_marginal_lib.R, written as a single loop so that
// no n x K temporaries are formed.
//
// Arguments
//   W, S      n x 3 matrices of declared weights and outcomes (NA where fewer
//             than three attempts)
//   mu, om    location and scale per unit
//   xi, kappa tail shape and log scale adjustment (sigma = SIG0 exp(kappa))
//   tau       attempt scale; delta: order effects (delta[0] = 0)
//   zq, q     splice point on the standardised scale and its probability
//   nodes, wts Gauss-Legendre nodes and weights on [-1, 1]
//   reach     half-width added to the attempt range, in kilograms
// Returns a vector of n unit log marginal likelihoods.

// Compiled without OpenMP, so that it builds with the Apple compiler; the
// pipeline runs one single-threaded job per core. On Linux, adding the line
//   // [[Rcpp::plugins(openmp)]]
// here makes the loops below run in parallel (setting WL_THREADS).
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static inline double log_plogis(double x) {
  // log(1 / (1 + exp(-x))) in a stable form
  if (x >= 0) return -log1p(std::exp(-x));
  return x - log1p(std::exp(x));
}

// Bulk: half-normal with scale 1 above the mode at 0, half-t with scale sL
// and nuL degrees of freedom below it, continuous at the mode; the GPD
// starts a fixed zq right-scales above the mode and carries the mass
// ptail = c (1 - qR) / 2, where c = 2 / (1 + sL phi(0) / t_nu(0)) is the
// density level at the mode and qR the quantile of the right half-normal at
// zq. 'q' is passed as ptail here (the mass in the tail), computed in R.
struct Bulk { double c, lc, left_mass, sL, nuL; };
static inline Bulk make_bulk(double sL, double nuL) {
  double t0 = std::exp(R::lgammafn((nuL + 1.0) / 2.0) - R::lgammafn(nuL / 2.0)) / std::sqrt(nuL * M_PI);
  Bulk b; b.sL = sL; b.nuL = nuL;
  b.c = 2.0 / (1.0 + sL * 0.3989422804014327 / t0);
  b.lc = std::log(b.c);
  b.left_mass = b.c * 0.3989422804014327 * sL / t0;
  return b;
}
static inline double log_f_std(double e, double xi, double sig, double zq, double ptail, const Bulk& B) {
  if (e < 0) {
    double u = e / B.sL;
    return B.lc - 0.9189385332046727 - ((B.nuL + 1.0) / 2.0) * log1p(u * u / B.nuL);
  }
  if (e <= zq) return B.lc - 0.5 * e * e - 0.9189385332046727;
  double y = (e - zq) / sig;
  if (std::fabs(xi) < 1e-6) return std::log(ptail) - std::log(sig) - y;
  double a = 1.0 + xi * y;
  if (a <= 0) return -INFINITY;
  return std::log(ptail) - std::log(sig) - (1.0 / xi + 1.0) * std::log(a);
}

static inline double surv_std(double e, double xi, double sig, double zq, double ptail, const Bulk& B) {
  if (e < 0) return 1.0 - B.left_mass * R::pt(e / B.sL, B.nuL, 1, 0);
  if (e <= zq) return B.c * R::pnorm(e, 0.0, 1.0, 0, 0);
  double y = (e - zq) / sig;
  if (std::fabs(xi) < 1e-6) return ptail * std::exp(-y);
  double a = 1.0 + xi * y;
  if (a <= 0) return 0.0;
  return ptail * std::pow(a, -1.0 / xi);
}

// Declaration model. A declared weight w is a noisy signal of capability:
// the relative margin r = (z - w)/w follows a two-piece normal with mode m
// and scales sm (below the mode, shared across situations) and sp (above).
// Situations: 1 opener, 2 raise after a successful attempt, 3 raise after a
// missed attempt. For situations 2 and 3 the declaration is truncated by the
// rule w >= w_prev + 1, which bounds r above by r_max = (z - w_prev - 1) /
// (w_prev + 1), and the density is renormalised by the two-piece normal cdf
// at r_max. dec = (m1, m2, m3, sp1, sp2, sp3, sm); sm <= 0 switches the
// declaration factors off.
static inline double tpn_cdf(double r, double m, double sm, double sp) {
  double c = 2.0 / (sm + sp);
  if (r < m) return c * sm * R::pnorm((r - m) / sm, 0.0, 1.0, 1, 0);
  return c * sm * 0.5 + c * sp * (R::pnorm((r - m) / sp, 0.0, 1.0, 1, 0) - 0.5);
}
static inline double log_declaration(double z, double w, double wprev, int sit, NumericVector& dec) {
  double sm = dec[6];
  if (sm <= 0 || sit == 0) return 0.0;
  double m = dec[sit - 1], sp = dec[3 + sit - 1];
  double r = (z - w) / w;
  double lc = std::log(2.0 / (sm + sp)) - 0.9189385332046727 - std::log(w);
  double u = (r < m) ? (r - m) / sm : (r - m) / sp;
  double v = lc - 0.5 * u * u;
  if (sit >= 2) {
    double rmax = (z - wprev - 1.0) / (wprev + 1.0);
    if (r > rmax) return -INFINITY;   // cannot happen: w >= wprev + 1 by construction
    double Gm = tpn_cdf(rmax, m, sm, sp);
    v -= std::log(Gm > 1e-300 ? Gm : 1e-300);
  }
  return v;
}

// [[Rcpp::export]]
NumericVector unit_ll_cpp(NumericMatrix W, IntegerMatrix S,
                          NumericVector mu, NumericVector om,
                          double xi, double kappa, double tau, NumericVector delta,
                          double zq, double q, double sig0, double sL, double nuL,
                          NumericVector nodes, NumericVector wts, double reach_mult,
                          IntegerMatrix D, NumericVector dec) {
  int n = W.nrow(), G = nodes.size();
  double sig = sig0 * std::exp(kappa);
  double reach = reach_mult * tau;
  Bulk B = make_bulk(sL, nuL);
  NumericVector out(n);
  std::vector<double> buf(2 * G);
  bool opener = dec[6] > 0;
  for (int i = 0; i < n; ++i) {
    // attempt summary
    bool has_miss = false, has_make = false;
    double wmade_min = INFINITY, wmiss_max = -INFINITY, wlast = -INFINITY;
    int m = 0;
    for (int j = 0; j < 3; ++j) {
      if (ISNAN(W(i, j))) continue;
      ++m; wlast = W(i, j);
      if (S(i, j) == 1) { has_make = true; if (W(i, j) < wmade_min) wmade_min = W(i, j); }
      else { has_miss = true; if (W(i, j) > wmiss_max) wmiss_max = W(i, j); }
    }
    double a = has_make ? wmade_min - reach : mu[i] - 7.0 * om[i] * sL;
    double b = has_miss ? wmiss_max + reach : wlast + reach;
    if (opener) {
      // declaration densities are negligible beyond 8 scales from the mode
      bool any_dec = false;
      for (int j = 0; j < 3; ++j) {
        if (ISNAN(W(i, j)) || D(i, j) == 0) continue;
        int sit = D(i, j);
        double m = dec[sit - 1], sp = dec[3 + sit - 1], sm = dec[6];
        double lo_h = W(i, j) * (1.0 + m - 8.0 * sm), hi_h = W(i, j) * (1.0 + m + 8.0 * sp);
        if (lo_h > a) a = lo_h;
        if (hi_h < b || !any_dec) { if (!has_miss && !any_dec) b = hi_h; else if (hi_h < b) b = hi_h; }
        any_dec = true;
      }
    }
    if (b < a + 1e-3) b = a + 1e-3;
    double zs = mu[i] + om[i] * zq;
    double mid = (zs > a && zs < b) ? zs : 0.5 * (a + b);
    double h[2] = {0.5 * (mid - a), 0.5 * (b - mid)};
    double c[2] = {0.5 * (mid + a), 0.5 * (b + mid)};
    double mx = -INFINITY;
    for (int piece = 0; piece < 2; ++piece) {
      for (int g = 0; g < G; ++g) {
        double z = c[piece] + h[piece] * nodes[g];
        double e = (z - mu[i]) / om[i];
        double lf = log_f_std(e, xi, sig, zq, q, B) - std::log(om[i]);
        double la = 0.0;
        for (int j = 0; j < 3; ++j) {
          if (ISNAN(W(i, j))) continue;
          double eta = (z - W(i, j)) / tau + delta[j];
          la += (S(i, j) == 1) ? log_plogis(eta) : log_plogis(-eta);
          if (opener) la += log_declaration(z, W(i, j), j > 0 ? W(i, j - 1) : 0.0, D(i, j), dec);
        }
        double lw = (h[piece] * wts[g] > 0) ? std::log(h[piece] * wts[g]) : -INFINITY;
        double v = lf + la + lw;
        buf[piece * G + g] = v;
        if (v > mx) mx = v;
      }
    }
    double integ = 0.0;
    if (std::isfinite(mx)) {
      for (int k = 0; k < 2 * G; ++k) integ += std::exp(buf[k] - mx);
      integ *= std::exp(mx);
    }
    double tail = (has_miss || opener) ? 0.0 : surv_std((b - mu[i]) / om[i], xi, sig, zq, q, B);
    double L = integ + tail;
    out[i] = std::log(L > 1e-300 ? L : 1e-300);
  }
  return out;
}

// Posterior summaries of the latent capability for each unit, given the
// global parameters: mean, and probability that Z exceeds a threshold
// (used for the capability gap and for the latent-tail diagnostics).
// [[Rcpp::export]]
NumericMatrix unit_latent_cpp(NumericMatrix W, IntegerMatrix S,
                              NumericVector mu, NumericVector om,
                              double xi, double kappa, double tau, NumericVector delta,
                              double zq, double q, double sig0, double sL, double nuL,
                              NumericVector nodes, NumericVector wts, double reach_mult,
                              NumericVector probs, IntegerMatrix D, NumericVector dec) {
  // returns n x (2 + length(probs)): posterior mean of Z, posterior sd, and
  // posterior quantiles of Z at the requested probabilities (by inverting the
  // quadrature cdf; the tail mass beyond b for no-miss units is placed by
  // sampling the conditional survival at its median for the mean)
  int n = W.nrow(), G = nodes.size(), P = probs.size();
  double sig = sig0 * std::exp(kappa);
  double reach = reach_mult * tau;
  Bulk B = make_bulk(sL, nuL);
  NumericMatrix out(n, 2 + P);
  std::vector<double> zz(2 * G), pp(2 * G);
  for (int i = 0; i < n; ++i) {
    bool has_miss = false, has_make = false;
    double wmade_min = INFINITY, wmiss_max = -INFINITY, wlast = -INFINITY;
    for (int j = 0; j < 3; ++j) {
      if (ISNAN(W(i, j))) continue;
      wlast = W(i, j);
      if (S(i, j) == 1) { has_make = true; if (W(i, j) < wmade_min) wmade_min = W(i, j); }
      else { has_miss = true; if (W(i, j) > wmiss_max) wmiss_max = W(i, j); }
    }
    double a = has_make ? wmade_min - reach : mu[i] - 7.0 * om[i] * sL;
    double b = has_miss ? wmiss_max + reach : wlast + reach;
    bool opener = dec[6] > 0;
    if (opener) {
      bool any_dec = false;
      for (int j = 0; j < 3; ++j) {
        if (ISNAN(W(i, j)) || D(i, j) == 0) continue;
        int sit = D(i, j);
        double m = dec[sit - 1], sp = dec[3 + sit - 1], sm = dec[6];
        double lo_h = W(i, j) * (1.0 + m - 8.0 * sm), hi_h = W(i, j) * (1.0 + m + 8.0 * sp);
        if (lo_h > a) a = lo_h;
        if (!has_miss && !any_dec) b = hi_h; else if (hi_h < b) b = hi_h;
        any_dec = true;
      }
    }
    if (b < a + 1e-3) b = a + 1e-3;
    double zs = mu[i] + om[i] * zq;
    double mid = (zs > a && zs < b) ? zs : 0.5 * (a + b);
    double h[2] = {0.5 * (mid - a), 0.5 * (b - mid)};
    double c[2] = {0.5 * (mid + a), 0.5 * (b + mid)};
    double mx = -INFINITY;
    for (int piece = 0; piece < 2; ++piece) {
      for (int g = 0; g < G; ++g) {
        double z = c[piece] + h[piece] * nodes[g];
        double e = (z - mu[i]) / om[i];
        double lf = log_f_std(e, xi, sig, zq, q, B) - std::log(om[i]);
        double la = 0.0;
        for (int j = 0; j < 3; ++j) {
          if (ISNAN(W(i, j))) continue;
          double eta = (z - W(i, j)) / tau + delta[j];
          la += (S(i, j) == 1) ? log_plogis(eta) : log_plogis(-eta);
          if (opener) la += log_declaration(z, W(i, j), j > 0 ? W(i, j - 1) : 0.0, D(i, j), dec);
        }
        double lw = (h[piece] * wts[g] > 0) ? std::log(h[piece] * wts[g]) : -INFINITY;
        double v = lf + la + lw;
        zz[piece * G + g] = z; pp[piece * G + g] = v;
        if (v > mx) mx = v;
      }
    }
    double integ = 0.0, m1 = 0.0, m2 = 0.0;
    for (int k = 0; k < 2 * G; ++k) {
      double w = std::exp(pp[k] - mx);
      integ += w; m1 += w * zz[k]; m2 += w * zz[k] * zz[k];
    }
    double tail = (has_miss || opener) ? 0.0 : surv_std((b - mu[i]) / om[i], xi, sig, zq, q, B) / std::exp(mx);
    // mean of Z beyond b: use the GPD/normal conditional mean numerically by a
    // coarse grid on the survival function (adequate: tail mass is small)
    double tail_m1 = 0.0, tail_m2 = 0.0;
    if (tail > 0) {
      // conditional distribution beyond b: sample 40 quantiles
      double Sb = surv_std((b - mu[i]) / om[i], xi, sig, zq, q, B);
      for (int k = 0; k < 40; ++k) {
        double u = 1.0 - Sb * (1.0 - (k + 0.5) / 40.0);  // quantile level above b
        // invert: standardised quantile
        double eq;
        double pbar = 1.0 - u;
        if (u <= q) eq = R::qnorm(u, 0.0, 1.0, 1, 0);
        else {
          double p = pbar / (1.0 - q);
          eq = (std::fabs(xi) < 1e-6) ? zq - sig * std::log(p) : zq + sig * (std::pow(p, -xi) - 1.0) / xi;
        }
        double z = mu[i] + om[i] * eq;
        tail_m1 += z / 40.0; tail_m2 += z * z / 40.0;
      }
    }
    double L = integ + tail;
    double mean = (m1 + tail * tail_m1) / L;
    double var = (m2 + tail * tail_m2) / L - mean * mean;
    out(i, 0) = mean; out(i, 1) = std::sqrt(var > 0 ? var : 0);
    // quantiles from the discrete cdf over nodes (tail mass appended at b)
    for (int p = 0; p < P; ++p) {
      double target = probs[p] * L, cum = 0.0; double zq_out = b;
      for (int k = 0; k < 2 * G; ++k) {
        cum += std::exp(pp[k] - mx);
        if (cum >= target) { zq_out = zz[k]; break; }
      }
      out(i, 2 + p) = zq_out;
    }
  }
  return out;
}

// -----------------------------------------------------------------------
// Bivariate stage: joint marginal likelihood of the two lifts of an
// athlete-day under a t copula (psi = 1/nu; psi = 0 is Gaussian) on the
// capability margins, with the latent capabilities integrated over the
// product of the two per-lift quadrature grids. Also returns, optionally,
// one posterior draw of (Z_s, Z_c) per unit.

static inline double cdf_std(double e, double xi, double sig, double zq, double q, const Bulk& B) {
  return 1.0 - surv_std(e, xi, sig, zq, q, B);
}

struct LiftGrid {
  std::vector<double> z, v, u;   // node, integrand value (with weight), cdf
};

static void lift_grid(int i, NumericMatrix& W, IntegerMatrix& S, NumericVector& mu, NumericVector& om,
                      double xi, double kappa, double tau, NumericVector& delta,
                      double zq, double q, double sig0, double sL, double nuL,
                      NumericVector& nodes, NumericVector& wts, double reach_mult,
                      IntegerMatrix& D, NumericVector& dec, LiftGrid& out) {
  int G = nodes.size();
  double sig = sig0 * std::exp(kappa);
  double reach = reach_mult * tau;
  Bulk B = make_bulk(sL, nuL);
  bool opener = dec[6] > 0;
  bool has_miss = false, has_make = false;
  double wmade_min = INFINITY, wmiss_max = -INFINITY, wlast = -INFINITY;
  for (int j = 0; j < 3; ++j) {
    if (ISNAN(W(i, j))) continue;
    wlast = W(i, j);
    if (S(i, j) == 1) { has_make = true; if (W(i, j) < wmade_min) wmade_min = W(i, j); }
    else { has_miss = true; if (W(i, j) > wmiss_max) wmiss_max = W(i, j); }
  }
  double a = has_make ? wmade_min - reach : mu[i] - 7.0 * om[i] * sL;
  double b = has_miss ? wmiss_max + reach : wlast + reach;
  if (opener) {
    bool any_dec = false;
    for (int j = 0; j < 3; ++j) {
      if (ISNAN(W(i, j)) || D(i, j) == 0) continue;
      int sit = D(i, j);
      double m = dec[sit - 1], sp = dec[3 + sit - 1], sm = dec[6];
      double lo_h = W(i, j) * (1.0 + m - 8.0 * sm), hi_h = W(i, j) * (1.0 + m + 8.0 * sp);
      if (lo_h > a) a = lo_h;
      if (!has_miss && !any_dec) b = hi_h; else if (hi_h < b) b = hi_h;
      any_dec = true;
    }
  } else if (!has_miss) {
    // without the opener model the mass beyond b is not negligible for
    // no-miss units; extend the range into the tail instead of adding the
    // analytic survival, so that the copula can be evaluated on the grid
    b = wlast + 8.0 * om[i];
  }
  if (b < a + 1e-3) b = a + 1e-3;
  double zs = mu[i] + om[i] * zq;
  double mid = (zs > a && zs < b) ? zs : 0.5 * (a + b);
  double h[2] = {0.5 * (mid - a), 0.5 * (b - mid)};
  double c[2] = {0.5 * (mid + a), 0.5 * (b + mid)};
  out.z.resize(2 * G); out.v.resize(2 * G); out.u.resize(2 * G);
  for (int piece = 0; piece < 2; ++piece) {
    for (int g = 0; g < G; ++g) {
      double z = c[piece] + h[piece] * nodes[g];
      double e = (z - mu[i]) / om[i];
      double lf = log_f_std(e, xi, sig, zq, q, B) - std::log(om[i]);
      double la = 0.0;
      for (int j = 0; j < 3; ++j) {
        if (ISNAN(W(i, j))) continue;
        double eta = (z - W(i, j)) / tau + delta[j];
        la += (S(i, j) == 1) ? log_plogis(eta) : log_plogis(-eta);
      }
      if (opener) {
        for (int j = 0; j < 3; ++j) {
          if (ISNAN(W(i, j))) continue;
          la += log_declaration(z, W(i, j), j > 0 ? W(i, j - 1) : 0.0, D(i, j), dec);
        }
      }
      double lw = (h[piece] * wts[g] > 0) ? std::log(h[piece] * wts[g]) : -INFINITY;
      int k = piece * G + g;
      out.z[k] = z; out.v[k] = lf + la + lw;
      double uu = cdf_std(e, xi, sig, zq, q, B);
      if (uu < 1e-12) uu = 1e-12; if (uu > 1 - 1e-12) uu = 1 - 1e-12;
      out.u[k] = uu;
    }
  }
}

static inline double log_tcop(double x, double y, double rho, double nu, double lconst, double d) {
  double Q = (x * x - 2.0 * rho * x * y + y * y) / d;
  return lconst - 0.5 * std::log(d) - ((nu + 2.0) / 2.0) * log1p(Q / nu)
    + ((nu + 1.0) / 2.0) * (log1p(x * x / nu) + log1p(y * y / nu));
}

// [[Rcpp::export]]
List joint_ll_cpp(NumericMatrix Ws, IntegerMatrix Ss, NumericMatrix Wc, IntegerMatrix Sc,
                  IntegerMatrix Ds, IntegerMatrix Dc,
                  NumericVector mus, NumericVector oms, NumericVector muc, NumericVector omc,
                  NumericVector ps, NumericVector pc,   // (xi, kappa, tau, d2, d3, zq, sig0, sL, ptail, nuL)
                  NumericVector decs, NumericVector decc,
                  double q, NumericVector nodes, NumericVector wts, double reach_mult,
                  double rho, double psi, bool draw) {
  int n = Ws.nrow(), G = nodes.size(), K = 2 * G;
  NumericVector out(n), zs_draw(n), zc_draw(n);
  LiftGrid gs, gc;
  NumericVector ds = NumericVector::create(0.0, ps[3], ps[4]);
  NumericVector dc = NumericVector::create(0.0, pc[3], pc[4]);
  bool gauss = psi < 1e-8;
  double nu = gauss ? 1e9 : 1.0 / psi;
  double d = 1.0 - rho * rho;
  double lconst = gauss ? 0.0 : (R::lgammafn((nu + 2.0) / 2.0) + R::lgammafn(nu / 2.0) - 2.0 * R::lgammafn((nu + 1.0) / 2.0));
  std::vector<double> xs(K), xc(K), cell(K * K);
  for (int i = 0; i < n; ++i) {
    lift_grid(i, Ws, Ss, mus, oms, ps[0], ps[1], ps[2], ds, ps[5], ps[8], ps[6], ps[7], ps[9], nodes, wts, reach_mult, Ds, decs, gs);
    lift_grid(i, Wc, Sc, muc, omc, pc[0], pc[1], pc[2], dc, pc[5], pc[8], pc[6], pc[7], pc[9], nodes, wts, reach_mult, Dc, decc, gc);
    // transform cdf values to the copula's standard scale
    for (int k = 0; k < K; ++k) {
      xs[k] = gauss ? R::qnorm(gs.u[k], 0.0, 1.0, 1, 0) : R::qt(gs.u[k], nu, 1, 0);
      xc[k] = gauss ? R::qnorm(gc.u[k], 0.0, 1.0, 1, 0) : R::qt(gc.u[k], nu, 1, 0);
    }
    double mx = -INFINITY;
    for (int g = 0; g < K; ++g) {
      if (!std::isfinite(gs.v[g])) { for (int h = 0; h < K; ++h) cell[g * K + h] = -INFINITY; continue; }
      for (int h = 0; h < K; ++h) {
        if (!std::isfinite(gc.v[h])) { cell[g * K + h] = -INFINITY; continue; }
        double lc;
        if (gauss) {
          double x = xs[g], y = xc[h];
          lc = -0.5 * std::log(d) - (rho * rho * (x * x + y * y) - 2.0 * rho * x * y) / (2.0 * d);
        } else lc = log_tcop(xs[g], xc[h], rho, nu, lconst, d);
        double v = gs.v[g] + gc.v[h] + lc;
        cell[g * K + h] = v;
        if (v > mx) mx = v;
      }
    }
    double tot = 0.0;
    for (int k = 0; k < K * K; ++k) if (std::isfinite(cell[k])) tot += std::exp(cell[k] - mx);
    out[i] = std::isfinite(mx) ? mx + std::log(tot) : -690.0;
    if (draw && std::isfinite(mx)) {
      double r = R::unif_rand() * tot, cum = 0.0; int pick = K * K - 1;
      for (int k = 0; k < K * K; ++k) {
        if (!std::isfinite(cell[k])) continue;
        cum += std::exp(cell[k] - mx);
        if (cum >= r) { pick = k; break; }
      }
      zs_draw[i] = gs.z[pick / K]; zc_draw[i] = gc.z[pick % K];
    }
  }
  return List::create(Named("ll") = out, Named("zs") = zs_draw, Named("zc") = zc_draw);
}

// =======================================================================
// Extensions
// -----------------------------------------------------------------------
// Declaration margins. Two parametrisations of the margin of capability over
// the declared weight are available (jac):
//   jac = 1  relative margin R = Z / w - 1, density of w given Z equal to
//            f(R) Z / w^2;
//   jac = 0  log margin R = log(Z / w), density f(R) / w.
// In both, R has a two-piece normal density whose mode and upper scale depend
// on the situation (dec = (m1, m2, m3, sp1, sp2, sp3, sm)), truncated by the
// rule w >= w_prev + 1 for raises (situations 2 and 3) and renormalised.
//
// (1) Shared self-assessment error. The athlete's own estimate of capability
//     on the day is Z_est = Z (1 + eta) (relative) or log Z_est = log Z + eta
//     (log), eta ~ N(0, s_eta^2) common to the day's declarations, whose
//     margins are taken relative to Z_est. Given Z the margin of declaration
//     j is a_j eta + b_j, the truncation point c_j eta + d_j; eta is
//     integrated analytically on each segment between the points where a
//     margin or a truncation point crosses its mode (where the two-piece
//     density switches scale), the product of the declaration densities and
//     the prior being an exact Gaussian in eta there, times a smooth factor
//     (the truncation normalisers and, for the relative margin, the Jacobian
//     (1 + eta)^n relative to exp(n eta)) which is expanded to second order
//     about the segment's centre, absorbed into the Gaussian, and integrated
//     by Gauss-Legendre after the substitution u = Phi((eta - mu) / sigma).
// (2) Athlete level. The days of one athlete share a Gaussian factor: the
//     standardised capability error e of every athlete-day keeps its
//     marginal distribution F, and c = Phi^{-1}(F(e)) is N(0,1) with
//     correlation rho_a between days of the same athlete, c = sqrt(rho_a) u +
//     sqrt(1 - rho_a) v. The likelihood of an athlete integrates u by
//     adaptive Gauss-Hermite quadrature, with the inner integrals over each
//     day's capability on the usual grid.

// dec = (m1, m2, m3, sp1, sp2, sp3, sm1, sm2, sm3): modes, upper scales and lower
// scales of the two-piece normal margin density by situation.
struct DecSetup {
  int nj; bool rel;
  double a[3], b[3], c[3], d[3], lcj[3], m[3], sp[3], sm[3]; int sit[3]; bool tr[3];
};
static inline double dec_sm(NumericVector& dec, int sit) { return dec.size() >= 9 ? dec[6 + sit - 1] : dec[6]; }
static inline void dec_setup(double z, const double* w, const double* wprev, const int* sit, int nj,
                             NumericVector& dec, int jac, DecSetup& ds) {
  ds.nj = nj; ds.rel = (jac == 1);
  for (int j = 0; j < nj; ++j) {
    ds.m[j] = dec[sit[j] - 1]; ds.sp[j] = dec[3 + sit[j] - 1]; ds.sm[j] = dec_sm(dec, sit[j]); ds.sit[j] = sit[j]; ds.tr[j] = (sit[j] >= 2);
    double lc = std::log(2.0 / (ds.sm[j] + ds.sp[j])) - 0.9189385332046727;
    if (ds.rel) {
      ds.a[j] = z / w[j]; ds.b[j] = z / w[j] - 1.0; ds.lcj[j] = lc + std::log(z) - 2.0 * std::log(w[j]);
      ds.c[j] = z / (wprev[j] + 1.0); ds.d[j] = z / (wprev[j] + 1.0) - 1.0;
    } else {
      ds.a[j] = 1.0; ds.b[j] = std::log(z / w[j]); ds.lcj[j] = lc - std::log(w[j]);
      ds.c[j] = 1.0; ds.d[j] = std::log(z / (wprev[j] + 1.0));
    }
  }
}

// base declaration factor (no shared error) for one declaration
static inline double log_declaration_j(const DecSetup& ds, int j, bool trunc) {
  double r = ds.b[j];
  double u = (r < ds.m[j]) ? (r - ds.m[j]) / ds.sm[j] : (r - ds.m[j]) / ds.sp[j];
  double v = ds.lcj[j] - 0.5 * u * u;
  if (trunc && ds.tr[j]) {
    double rmax = ds.d[j];
    if (r > rmax + 1e-9) return -INFINITY;
    double Gm = tpn_cdf(rmax, ds.m[j], ds.sm[j], ds.sp[j]);
    v -= std::log(Gm > 1e-300 ? Gm : 1e-300);
  }
  return v;
}

// closed form without truncation (piecewise Gaussian in eta)
static double log_eta_integral(const DecSetup& ds, double s_eta) {
  int nj = ds.nj;
  double etas[3];
  for (int j = 0; j < nj; ++j) etas[j] = (ds.m[j] - ds.b[j]) / ds.a[j];
  double bp[3]; int ord[3];
  for (int j = 0; j < nj; ++j) { bp[j] = etas[j]; ord[j] = j; }
  for (int i = 0; i < nj; ++i) for (int k = i + 1; k < nj; ++k)
    if (bp[ord[k]] < bp[ord[i]]) { int t = ord[i]; ord[i] = ord[k]; ord[k] = t; }
  double total = -INFINITY;
  for (int piece = 0; piece <= nj; ++piece) {
    double lo = (piece == 0) ? -INFINITY : bp[ord[piece - 1]];
    double hi = (piece == nj) ? INFINITY : bp[ord[piece]];
    if (hi <= lo) continue;
    double A = 1.0 / (2.0 * s_eta * s_eta), B = 0.0, C = -std::log(s_eta) - 0.9189385332046727;
    if (ds.rel) { A += 0.5 * nj; B += nj; }   // (1 + eta)^nj as exp(nj eta - nj eta^2 / 2)
    for (int j = 0; j < nj; ++j) {
      double sc = (hi <= etas[j]) ? ds.sm[j] : ds.sp[j];
      double k2 = ds.a[j] * ds.a[j] / (sc * sc);
      A += 0.5 * k2; B += k2 * etas[j]; C += ds.lcj[j] - 0.5 * k2 * etas[j] * etas[j];
    }
    double mu = B / (2.0 * A), sd = std::sqrt(1.0 / (2.0 * A));
    double lpiece = C + B * B / (4.0 * A) + 0.5 * std::log(M_PI / A);
    double plo = std::isfinite(lo) ? R::pnorm((lo - mu) / sd, 0.0, 1.0, 1, 0) : 0.0;
    double phi = std::isfinite(hi) ? R::pnorm((hi - mu) / sd, 0.0, 1.0, 1, 0) : 1.0;
    double mass = phi - plo;
    if (mass <= 0) continue;
    double v = lpiece + std::log(mass);
    if (total == -INFINITY) total = v; else total = std::max(total, v) + log1p(std::exp(-std::fabs(total - v)));
  }
  return total;
}

static inline void tpn_logF_derivs(double r, double m, double sm, double sp, double& lF, double& d1, double& d2) {
  // log F(r), d/dr log F, d^2/dr^2 log F for the two-piece normal
  double c = 2.0 / (sm + sp);
  double F = tpn_cdf(r, m, sm, sp);
  if (F < 1e-300) F = 1e-300;
  double s = (r < m) ? sm : sp;
  double x = (r - m) / s;
  double f = c * R::dnorm(x, 0.0, 1.0, 0);
  double fp = -f * x / s;
  lF = std::log(F); d1 = f / F; d2 = fp / F - d1 * d1;
}

// log Phi(x) + x^2 / 2 and its first two derivatives (smooth and slowly varying
// for x < 0, where log Phi(x) itself is close to -x^2 / 2)
static inline void lowtail_g(double x, double& g, double& g1, double& g2) {
  double lPhi = R::pnorm(x, 0.0, 1.0, 1, 1);
  double lam = std::exp(-0.5 * x * x - 0.9189385332046727 - lPhi);   // phi(x) / Phi(x)
  g = lPhi + 0.5 * x * x; g1 = x + lam; g2 = 1.0 - lam * (x + lam);
}

// with truncation: exact Gaussian segments times a smooth factor.
// Where the truncation point lies below the mode of a raise's density,
// -log F(rmax) = x^2 / 2 - log(c sm) - g(x) with x = (rmax - m) / sm: the
// quadratic part is taken into the Gaussian of the segment exactly (it
// cancels the declaration's own Gaussian for log margins) and only the slowly
// varying remainder goes into the expanded factor. Expanding -log F itself
// there, as an earlier version did, loses most of the mass when the
// declaration scale is small relative to s_eta.
static double log_eta_integral_pt(const DecSetup& ds, double s_eta, const double* glx, const double* glw, int nGL) {
  int nj = ds.nj;
  const int MAXB = 40;
  double etas[3], bnd[MAXB]; int nb = 0;
  for (int j = 0; j < nj; ++j) {
    etas[j] = (ds.m[j] - ds.b[j]) / ds.a[j];
    bnd[nb++] = etas[j];
    if (ds.tr[j]) {
      double et = (ds.m[j] - ds.d[j]) / ds.c[j];
      bnd[nb++] = et;
      // around the point where the truncation point crosses the mode, log F
      // changes over a few lower scales (steeply above it, logarithmically below
      // it): extra breakpoints at geometric distances on both sides keep the
      // remaining factor close to quadratic on every segment
      for (int q = 0; q < 7 && nb < MAXB - 1; ++q) {
        double h = 0.5 * std::ldexp(1.0, q) * ds.sm[j] / ds.c[j];
        if (std::fabs(et + h) < 8.0 * s_eta) bnd[nb++] = et + h;
        if (std::fabs(et - h) < 8.0 * s_eta) bnd[nb++] = et - h;
      }
    }
  }
  for (int i = 0; i < nb; ++i) for (int k = i + 1; k < nb; ++k) if (bnd[k] < bnd[i]) { double t = bnd[i]; bnd[i] = bnd[k]; bnd[k] = t; }
  const int MAXS = MAXB + 1;
  double lmass[MAXS], mu_s[MAXS], sg_s[MAXS], lo_s[MAXS], hi_s[MAXS], lconst[MAXS], e0_s[MAXS], T0_s[MAXS], T1_s[MAXS], T2_s[MAXS];
  bool low_s[MAXS][3];
  double lmax = -INFINITY;
  for (int sgm = 0; sgm <= nb; ++sgm) {
    double lo = (sgm == 0) ? -INFINITY : bnd[sgm - 1];
    double hi = (sgm == nb) ? INFINITY : bnd[sgm];
    lmass[sgm] = -INFINITY;
    if (hi <= lo) continue;
    double A = 1.0 / (2.0 * s_eta * s_eta), B = 0.0, C = -std::log(s_eta) - 0.9189385332046727;
    if (ds.rel) B += nj;   // exp(nj eta); the remainder of (1 + eta)^nj goes into T
    for (int j = 0; j < nj; ++j) {
      double sc = (hi <= etas[j]) ? ds.sm[j] : ds.sp[j];
      double k2 = ds.a[j] * ds.a[j] / (sc * sc);
      A += 0.5 * k2; B += k2 * etas[j]; C += ds.lcj[j] - 0.5 * k2 * etas[j] * etas[j];
    }
    // truncation points below the mode in this segment: quadratic part of -log F into the Gaussian
    double em = (std::isfinite(lo) && std::isfinite(hi)) ? 0.5 * (lo + hi) : (std::isfinite(lo) ? lo + 1.0 : (std::isfinite(hi) ? hi - 1.0 : 0.0));
    for (int j = 0; j < 3; ++j) low_s[sgm][j] = false;
    for (int j = 0; j < nj; ++j) {
      if (!ds.tr[j] || !(ds.c[j] * em + ds.d[j] < ds.m[j])) continue;
      double kt = ds.c[j] * ds.c[j] / (ds.sm[j] * ds.sm[j]), et = (ds.m[j] - ds.d[j]) / ds.c[j];
      if (!(A - 0.5 * kt > 0.0)) continue;   // cannot happen for log margins; keep the plain expansion otherwise
      low_s[sgm][j] = true;
      A -= 0.5 * kt; B -= kt * et; C += 0.5 * kt * et * et;
    }
    double mu = B / (2.0 * A), sg = std::sqrt(1.0 / (2.0 * A));
    double e0 = mu; if (e0 < lo) e0 = lo; if (e0 > hi) e0 = hi;
    if (e0 <= -0.99) e0 = -0.99;
    double T0 = 0.0, T1 = 0.0, T2 = 0.0;
    if (ds.rel) { T0 = nj * (std::log1p(e0) - e0); T1 = nj * (1.0 / (1.0 + e0) - 1.0); T2 = -nj / ((1.0 + e0) * (1.0 + e0)); }
    for (int j = 0; j < nj; ++j) {
      if (!ds.tr[j]) continue;
      double r0 = ds.c[j] * e0 + ds.d[j];
      if (low_s[sgm][j]) {
        double q = ds.c[j] / ds.sm[j], g, g1, g2;
        lowtail_g((r0 - ds.m[j]) / ds.sm[j], g, g1, g2);
        T0 += -std::log(2.0 * ds.sm[j] / (ds.sm[j] + ds.sp[j])) - g; T1 -= q * g1; T2 -= q * q * g2;
        continue;
      }
      double lF, d1, d2;
      tpn_logF_derivs(r0, ds.m[j], ds.sm[j], ds.sp[j], lF, d1, d2);
      T0 -= lF; T1 -= ds.c[j] * d1; T2 -= ds.c[j] * ds.c[j] * d2;
    }
    double A2 = A - 0.5 * T2;
    if (!(A2 > 0.25 * A)) { A2 = A; T2 = 0.0; }
    double B2 = B + T1 - T2 * e0, C2 = C + T0 - T1 * e0 + 0.5 * T2 * e0 * e0;
    mu = B2 / (2.0 * A2); sg = std::sqrt(1.0 / (2.0 * A2));
    double lo2 = lo, hi2 = hi;
    if (lo2 < mu - 7.0 * sg) lo2 = mu - 7.0 * sg;
    if (hi2 > mu + 7.0 * sg) hi2 = mu + 7.0 * sg;
    if (ds.rel && lo2 <= -0.99) lo2 = -0.99;
    if (hi2 <= lo2) continue;
    double l = (lo2 - mu) / sg, h = (hi2 - mu) / sg;
    double mass = (l + h > 0) ? R::pnorm(l, 0.0, 1.0, 0, 0) - R::pnorm(h, 0.0, 1.0, 0, 0)
                              : R::pnorm(h, 0.0, 1.0, 1, 0) - R::pnorm(l, 0.0, 1.0, 1, 0);
    if (mass <= 1e-300) continue;
    lconst[sgm] = C2 + B2 * B2 / (4.0 * A2) + 0.5 * std::log(M_PI / A2);
    lmass[sgm] = lconst[sgm] + std::log(mass);
    mu_s[sgm] = mu; sg_s[sgm] = sg; lo_s[sgm] = lo2; hi_s[sgm] = hi2;
    e0_s[sgm] = e0; T0_s[sgm] = T0; T1_s[sgm] = T1; T2_s[sgm] = T2;
    if (lmass[sgm] > lmax) lmax = lmass[sgm];
  }
  if (!std::isfinite(lmax)) return -INFINITY;
  double tot = -INFINITY;
  for (int sgm = 0; sgm <= nb; ++sgm) {
    if (lmass[sgm] < lmax - 25.0) continue;
    double mu = mu_s[sgm], sg = sg_s[sgm];
    double l = (lo_s[sgm] - mu) / sg, h = (hi_s[sgm] - mu) / sg;
    bool upper = (l + h > 0);
    double u_lo, u_hi;
    if (upper) { u_lo = R::pnorm(h, 0.0, 1.0, 0, 0); u_hi = R::pnorm(l, 0.0, 1.0, 0, 0); }
    else { u_lo = R::pnorm(l, 0.0, 1.0, 1, 0); u_hi = R::pnorm(h, 0.0, 1.0, 1, 0); }
    double hw = 0.5 * (u_hi - u_lo), cm = 0.5 * (u_hi + u_lo);
    if (hw <= 0) continue;
    double e0 = e0_s[sgm], T0 = T0_s[sgm], T1 = T1_s[sgm], T2 = T2_s[sgm];
    double acc = -INFINITY;
    for (int k = 0; k < nGL; ++k) {
      double u = cm + hw * glx[k];
      double t = R::qnorm(u, 0.0, 1.0, 1, 0);
      if (upper) t = -t;
      double eta = mu + sg * t;
      double lT = ds.rel ? nj * (std::log1p(eta) - eta) : 0.0;
      for (int j = 0; j < nj; ++j) {
        if (!ds.tr[j]) continue;
        double rmax = ds.c[j] * eta + ds.d[j];
        double r = ds.a[j] * eta + ds.b[j];
        if (r > rmax + 1e-9) { lT = -INFINITY; break; }
        if (low_s[sgm][j]) {
          double g, g1, g2;
          lowtail_g((rmax - ds.m[j]) / ds.sm[j], g, g1, g2);
          lT += -std::log(2.0 * ds.sm[j] / (ds.sm[j] + ds.sp[j])) - g;
          continue;
        }
        double Gm = tpn_cdf(rmax, ds.m[j], ds.sm[j], ds.sp[j]);
        lT -= std::log(Gm > 1e-300 ? Gm : 1e-300);
      }
      if (!std::isfinite(lT)) continue;
      double de = eta - e0;
      double lR = lT - (T0 + T1 * de + 0.5 * T2 * de * de);
      double lt = std::log(hw * glw[k]) + lR;
      if (acc == -INFINITY) acc = lt; else acc = std::max(acc, lt) + log1p(std::exp(-std::fabs(acc - lt)));
    }
    if (acc == -INFINITY) continue;
    double v = lconst[sgm] + acc;
    if (tot == -INFINITY) tot = v; else tot = std::max(tot, v) + log1p(std::exp(-std::fabs(tot - v)));
  }
  return tot;
}

// Gauss-Legendre nodes and weights on [-1, 1] by Newton iteration on the
// Legendre polynomial (nodes in increasing order), cached by size.
static void gauleg(int n, std::vector<double>& x, std::vector<double>& w) {
  x.assign(n, 0.0); w.assign(n, 0.0);
  int m = (n + 1) / 2;
  for (int i = 0; i < m; ++i) {
    double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), z1, pp;
    do {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 1; j <= n; ++j) { double p3 = p2; p2 = p1; p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j; }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z; z = z1 - p1 / pp;
    } while (std::fabs(z - z1) > 1e-14);
    x[i] = -z; x[n - 1 - i] = z;
    w[i] = 2.0 / ((1.0 - z * z) * pp * pp); w[n - 1 - i] = w[i];
  }
}
struct GLTable {
  std::vector<double> x[128], w[128];
  bool ready[128];
  GLTable() { for (int i = 0; i < 128; ++i) ready[i] = false; }
  int get(int n, const double*& xp, const double*& wp) {   // returns the number of nodes actually provided
    if (n < 1) n = 1; if (n > 127) n = 127;
    if (!ready[n]) {
#ifdef _OPENMP
#pragma omp critical
#endif
      { if (!ready[n]) { gauleg(n, x[n], w[n]); ready[n] = true; } }
    }
    xp = x[n].data(); wp = w[n].data();
    return n;
  }
};
static GLTable GL_TABLE;

// Quadrature layout of one unit: breakpoints at the range ends, at the
// splice point of the body, and at the modes of the declaration densities
// (where their scale switches from sm to sp); each sub-interval takes a
// number of Gauss-Legendre nodes set by its width relative to the narrowest
// declaration scale on the kilogram scale (wide sub-intervals take G nodes).
// Returns the number of nodes; fills zk (nodes) and lwk (log weights).
// largest node spacing in kilograms
static const double SPACING_MAX = 1.0;
static int layout_nodes(double a, double b, double zs, const double* wdec, const int* sdec, int nj, NumericVector& dec,
                        int G, double* zk, double* lwk, int kmax, bool rel) {
  double bp[12]; int nb = 0;
  bp[nb++] = a; bp[nb++] = b;
  if (zs > a && zs < b) bp[nb++] = zs;
  double smin = INFINITY;
  for (int j = 0; j < nj; ++j) { double sp = dec[3 + sdec[j] - 1]; if (sp * wdec[j] < smin) smin = sp * wdec[j]; }
  double spacing = std::isfinite(smin) ? 0.7 * smin : 1.0;
  if (spacing > SPACING_MAX) spacing = SPACING_MAX;
  if (spacing < 0.25) spacing = 0.25;
  double zlo[3], zst[3], sloc[3]; bool band[3];
  for (int j = 0; j < nj; ++j) {
    double m = dec[sdec[j] - 1], sm = dec_sm(dec, sdec[j]);
    double zstar = rel ? wdec[j] * (1.0 + m) : wdec[j] * std::exp(m);
    if (zstar > a && zstar < b) bp[nb++] = zstar;
    // a lower side of the margin density that is narrow compared with the node spacing gets its own
    // segment, from four lower scales below the mode to the mode, with nodes on the lower scale
    band[j] = sm * wdec[j] < 2.0 * spacing;
    double zlow = rel ? wdec[j] * (1.0 + m - 4.0 * sm) : wdec[j] * std::exp(m - 4.0 * sm);
    if (band[j] && zlow > a && zlow < b) bp[nb++] = zlow;
    zlo[j] = zlow; zst[j] = zstar; sloc[j] = std::max(0.7 * sm * wdec[j], 0.25);
  }
  for (int i = 0; i < nb; ++i) for (int k = i + 1; k < nb; ++k) if (bp[k] < bp[i]) { double t = bp[i]; bp[i] = bp[k]; bp[k] = t; }
  int nk = 0;
  for (int sgm = 0; sgm + 1 < nb; ++sgm) {
    double l = bp[sgm], h = bp[sgm + 1], wd = h - l;
    if (wd <= 1e-9) continue;
    double sg = spacing;
    for (int j = 0; j < nj; ++j) if (band[j] && l >= zlo[j] - 1e-9 && h <= zst[j] + 1e-9 && sloc[j] < sg) sg = sloc[j];
    int n;
    if (wd > 10.0) n = G;
    else { n = (int) std::ceil(wd / sg); if (n < 4) n = 4; if (n > 16) n = 16; }
    if (nk + n > kmax) n = kmax - nk;
    if (n <= 0) break;
    const double* gx; const double* gw;
    n = GL_TABLE.get(n, gx, gw);
    double hw = 0.5 * wd, cm = 0.5 * (l + h);
    for (int g = 0; g < n; ++g) { zk[nk] = cm + hw * gx[g]; lwk[nk] = std::log(hw * gw[g]); ++nk; }
  }
  return nk;
}

// Normal-mixture body with a generalized Pareto tail above an estimated
// threshold (the extreme value mixture model of Behrens et al. 2004 and
// Nascimento et al. 2012 with normal in place of gamma components, since
// the standardized error lives on the real line). For e <= u the density
// is H'(e) = sum_k w_k phi((e - m_k)/s_k)/s_k; above u it is
// (1 - H(u)) g(e - u; sigma, xi). bm = (K, w_1..K, m_1..K, s_1..K, u, sigma, xi).
struct Mix { int K; double w[6], m[6], s[6]; double u, sig, xi, l1mHu, S_u, lo; };
static bool make_mix(const NumericVector& bm, Mix& M) {
  if (bm.size() < 7) return false;
  int K = (int) bm[0]; M.K = K;
  for (int k = 0; k < K; ++k) { M.w[k] = bm[1 + k]; M.m[k] = bm[1 + K + k]; M.s[k] = bm[1 + 2 * K + k]; }
  M.u = bm[1 + 3 * K]; M.sig = bm[2 + 3 * K]; M.xi = bm[3 + 3 * K];
  double Su = 0.0, lo = INFINITY;
  for (int k = 0; k < K; ++k) { Su += M.w[k] * R::pnorm((M.u - M.m[k]) / M.s[k], 0.0, 1.0, 0, 0); double l = M.m[k] - 8.0 * M.s[k]; if (l < lo) lo = l; }
  M.S_u = Su; M.l1mHu = std::log(Su > 1e-300 ? Su : 1e-300); M.lo = lo;
  return true;
}
static inline double log_f_mix(double e, const Mix& M) {
  if (e <= M.u) {
    double lt[6], mx = -INFINITY;
    for (int k = 0; k < M.K; ++k) { double z = (e - M.m[k]) / M.s[k]; lt[k] = std::log(M.w[k]) - std::log(M.s[k]) - 0.9189385332046727 - 0.5 * z * z; if (lt[k] > mx) mx = lt[k]; }
    double sm = 0.0; for (int k = 0; k < M.K; ++k) sm += std::exp(lt[k] - mx);
    return mx + std::log(sm);
  }
  double y = (e - M.u) / M.sig;
  if (std::fabs(M.xi) < 1e-8) return M.l1mHu - std::log(M.sig) - y;
  double a = 1.0 + M.xi * y;
  if (a <= 0) return -INFINITY;
  return M.l1mHu - std::log(M.sig) - (1.0 / M.xi + 1.0) * std::log(a);
}
static inline double surv_mix(double e, const Mix& M) {
  if (e <= M.u) { double S = 0.0; for (int k = 0; k < M.K; ++k) S += M.w[k] * R::pnorm((e - M.m[k]) / M.s[k], 0.0, 1.0, 0, 0); return S; }
  double y = (e - M.u) / M.sig;
  if (std::fabs(M.xi) < 1e-8) return M.S_u * std::exp(-y);
  double a = 1.0 + M.xi * y;
  if (a <= 0) return 0.0;
  return M.S_u * std::pow(a, -1.0 / M.xi);
}
// [[Rcpp::export]]
NumericVector log_f_mix_cpp(NumericVector e, NumericVector bm) { Mix M; make_mix(bm, M); NumericVector o(e.size()); for (int i = 0; i < e.size(); ++i) o[i] = log_f_mix(e[i], M); return o; }
// [[Rcpp::export]]
NumericVector surv_mix_cpp(NumericVector e, NumericVector bm) { Mix M; make_mix(bm, M); NumericVector o(e.size()); for (int i = 0; i < e.size(); ++i) o[i] = surv_mix(e[i], M); return o; }

// [[Rcpp::export]]
NumericVector unit_ll_ext_cpp(NumericMatrix W, IntegerMatrix S, IntegerMatrix D,
                              NumericVector mu, NumericVector om,
                              double xi, double kappa, double tau, NumericVector delta,
                              double zq, double ptail, double sig0, double sL, double nuL,
                              NumericVector nodes, NumericVector wts, double reach_mult,
                              NumericVector dec, double s_eta, double rho_a,
                              IntegerVector astart, IntegerVector alen,
                              NumericVector ghx, NumericVector ghw,
                              NumericVector ghx_e, NumericVector ghw_e, double inflate, int trunc, int nthreads, int jac, int rtau,
                              double pi_miss, NumericVector bmix) {
  // bmix: normal-mixture body (see make_mix); empty for the half-t / half-normal body
  // rtau = 1: the attempt scale is relative, P(S = 1 | Z, w) = F(log(Z / w) / tau + delta_j); 0: absolute, F((Z - w) / tau + delta_j)
  // pi_miss in (0, 1): technical misses, P(S = 1 | Z, w) = (1 - pi_miss) F(.); pi_miss <= 0 switches them off
  int G = nodes.size(), K = 3 * G + 5 * 16 + 8, nA = astart.size(), nGH = ghx.size(), nGHe = ghx_e.size();
  bool tr = trunc > 0, rel = (jac == 1), rt = (rtau == 1), floor_on = pi_miss > 0;
  double l1mp = floor_on ? std::log1p(-pi_miss) : 0.0;
  double sig = sig0 * std::exp(kappa);
  double reach = reach_mult * tau;
  Bulk B = make_bulk(sL, nuL);
  Mix M; bool mix = make_mix(bmix, M);
  bool shared = s_eta > 0, copula = rho_a > 0;
  double sr = std::sqrt(rho_a), s1r = std::sqrt(1.0 - rho_a);
  NumericVector out(nA);
  const double* Wp = W.begin(); const int* Sp = S.begin(); const int* Dp = D.begin(); int nrow = W.nrow();
  double* outp = out.begin();
#ifdef _OPENMP
#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for schedule(dynamic, 32)
#endif
  for (int ia = 0; ia < nA; ++ia) {
    int n_d = alen[ia];
    std::vector<double> v((size_t) n_d * K, 0.0), cc((size_t) n_d * K, 0.0), cbar(n_d, 0.0), zk(K), lwk(K);
    std::vector<int> nk(n_d, 0);
    #define W(i, j) Wp[(i) + (size_t) nrow * (j)]
    #define S(i, j) Sp[(i) + (size_t) nrow * (j)]
    #define D(i, j) Dp[(i) + (size_t) nrow * (j)]
    double ll_indep = 0.0;
    for (int dd = 0; dd < n_d; ++dd) {
      int i = astart[ia] + dd;
      bool has_miss = false, has_make = false;
      double wmade_min = INFINITY, wmiss_max = -INFINITY, wlast = -INFINITY;
      double wdec[3], wpdec[3]; int sdec[3]; int nj = 0;
      for (int j = 0; j < 3; ++j) {
        if (ISNAN(W(i, j))) continue;
        wlast = W(i, j);
        if (S(i, j) == 1) { has_make = true; if (W(i, j) < wmade_min) wmade_min = W(i, j); }
        else { has_miss = true; if (W(i, j) > wmiss_max) wmiss_max = W(i, j); }
        if (D(i, j) > 0) { wdec[nj] = W(i, j); wpdec[nj] = j > 0 ? W(i, j - 1) : 0.0; sdec[nj] = D(i, j); ++nj; }
      }
      double reach_lo = rt ? reach * wmade_min : reach, reach_hi = rt ? reach * (has_miss ? wmiss_max : wlast) : reach;
      double a = has_make ? wmade_min - reach_lo : (mix ? mu[i] + om[i] * M.lo : mu[i] - 7.0 * om[i] * sL);
      double b = has_miss ? wmiss_max + reach_hi : wlast + reach_hi;
      // declaration supports, widened for the shared error
      double wid_lo = shared ? (1.0 - 5.0 * s_eta) : 1.0, wid_hi = shared ? (1.0 + 5.0 * s_eta) : 1.0;
      bool any_dec = false;
      for (int j = 0; j < nj; ++j) {
        double m = dec[sdec[j] - 1], sp = dec[3 + sdec[j] - 1], sm = dec_sm(dec, sdec[j]);
        double lo_h = rel ? wdec[j] * (1.0 + m - 8.0 * sm) * wid_lo : wdec[j] * std::exp(m - 8.0 * sm) * wid_lo;
        double hi_h = rel ? wdec[j] * (1.0 + m + 8.0 * sp) * wid_hi : wdec[j] * std::exp(m + 8.0 * sp) * wid_hi;
        // a truncated raise does not bound capability from below: at the smallest allowed raise the
        // density divided by its truncation probability does not decay below the mode
        if (!(tr && sdec[j] >= 2) && lo_h > a) a = lo_h;
        if (!has_miss && !any_dec) b = hi_h; else if (hi_h < b) b = hi_h;
        any_dec = true;
      }
      if (b < a + 1e-3) b = a + 1e-3;
      double zs = mix ? mu[i] + om[i] * M.u : mu[i] + om[i] * zq;
      int nkd = layout_nodes(a, b, zs, wdec, sdec, nj, dec, G, zk.data(), lwk.data(), K, rel);
      nk[dd] = nkd;
      double mx = -INFINITY;
      {
        for (int k = 0; k < nkd; ++k) {
          double z = zk[k];
          double e = (z - mu[i]) / om[i];
          double lf = (mix ? log_f_mix(e, M) : log_f_std(e, xi, sig, zq, ptail, B)) - std::log(om[i]);
          double la = 0.0;
          for (int j = 0; j < 3; ++j) {
            if (ISNAN(W(i, j))) continue;
            double eta = rt ? ((z > 0) ? std::log(z / W(i, j)) / tau : -INFINITY) + delta[j] : (z - W(i, j)) / tau + delta[j];
            if (!floor_on) la += (S(i, j) == 1) ? log_plogis(eta) : log_plogis(-eta);
            else la += (S(i, j) == 1) ? l1mp + log_plogis(eta) : std::log(pi_miss + (1.0 - pi_miss) * std::exp(log_plogis(-eta)));
          }
          double ld = 0.0;
          if (nj > 0) {
            DecSetup ds; dec_setup(z, wdec, wpdec, sdec, nj, dec, jac, ds);
            if (shared) ld = tr ? log_eta_integral_pt(ds, s_eta, ghx_e.begin(), ghw_e.begin(), nGHe) : log_eta_integral(ds, s_eta);
            else for (int j = 0; j < nj; ++j) ld += log_declaration_j(ds, j, tr);
          }
          double lw = lwk[k];
          double val = lf + la + ld + lw;
          v[(size_t) dd * K + k] = val;
          if (val > mx) mx = val;
          if (copula) {
            double Fe = 1.0 - (mix ? surv_mix(e, M) : surv_std(e, xi, sig, zq, ptail, B));
            if (Fe < 1e-12) Fe = 1e-12; if (Fe > 1.0 - 1e-12) Fe = 1.0 - 1e-12;
            cc[(size_t) dd * K + k] = R::qnorm(Fe, 0.0, 1.0, 1, 0);
          }
        }
      }
      // marginal (independent-days) log likelihood of the unit, and the
      // posterior mean of c under it, for centring the outer quadrature
      double sum = 0.0, csum = 0.0;
      for (int k = 0; k < nk[dd]; ++k) {
        double val = v[(size_t) dd * K + k];
        if (!std::isfinite(val)) continue;
        double wgt = std::exp(val - mx);
        sum += wgt; if (copula) csum += wgt * cc[(size_t) dd * K + k];
      }
      double ll_d = std::isfinite(mx) ? mx + std::log(sum) : -690.0;
      ll_indep += ll_d;
      cbar[dd] = (copula && sum > 0) ? csum / sum : 0.0;
    }
    if (!copula || n_d == 1) { outp[ia] = ll_indep; continue; }
    // adaptive Gauss-Hermite over the shared factor u ~ N(0,1)
    double prec = 1.0 + n_d * rho_a / (1.0 - rho_a), msum = 0.0;
    for (int dd = 0; dd < n_d; ++dd) msum += cbar[dd];
    double mu_u = (sr / (1.0 - rho_a)) * msum / prec, sd_u = std::sqrt(1.0 / prec);
    double tot = -INFINITY;
    for (int kk = 0; kk < nGH; ++kk) {
      double u = mu_u + std::sqrt(2.0) * sd_u * ghx[kk];
      double lterm = std::log(ghw[kk]) + ghx[kk] * ghx[kk] + 0.5 * std::log(2.0) + std::log(sd_u)
                     - 0.5 * u * u - 0.9189385332046727;   // GH weight, e^{x^2}, Jacobian, N(0,1) prior
      for (int dd = 0; dd < n_d; ++dd) {
        double mxd = -INFINITY;
        for (int k = 0; k < nk[dd]; ++k) {
          double val = v[(size_t) dd * K + k];
          if (!std::isfinite(val)) continue;
          double cz = cc[(size_t) dd * K + k];
          double r = (cz - sr * u) / s1r;
          double lcond = -0.5 * r * r + 0.5 * cz * cz - std::log(s1r);   // log phi(r) - log phi(c) - log sqrt(1-rho)
          double t = val + lcond;
          if (t > mxd) mxd = t;
        }
        double sd2 = 0.0;
        for (int k = 0; k < nk[dd]; ++k) {
          double val = v[(size_t) dd * K + k];
          if (!std::isfinite(val)) continue;
          double cz = cc[(size_t) dd * K + k];
          double r = (cz - sr * u) / s1r;
          sd2 += std::exp(val - 0.5 * r * r + 0.5 * cz * cz - std::log(s1r) - mxd);
        }
        lterm += std::isfinite(mxd) ? mxd + std::log(sd2) : -690.0;
      }
      if (tot == -INFINITY) tot = lterm; else tot = std::max(tot, lterm) + log1p(std::exp(-std::fabs(tot - lterm)));
    }
    outp[ia] = tot;
    #undef W
    #undef S
    #undef D
  }
  return out;
}


// =======================================================================
// Two-lift copula stage through the extended marginal model: per-unit node
// values and cdf values of each lift on the adaptive layout, combined under
// a t (psi = 1 / nu > 0) or Gaussian (psi = 0) copula. The athlete-level
// factor is not part of this stage. par vectors: (xi, kappa, tau, d2, d3,
// zq, sig0, sL, ptail, nuL, pi_miss); dec vectors of length 9.
struct LiftGridX { std::vector<double> z, v, u; int nk; };
static void lift_grid_ext(int i, const double* Wp, const int* Sp, const int* Dp, int nrow,
                          const double* mu, const double* om, NumericVector& par, NumericVector& dec,
                          int G, double reach_mult, int jac, int rtau, bool trunc, LiftGridX& out, const NumericVector& bmix) {
  #define W(i, j) Wp[(i) + (size_t) nrow * (j)]
  #define S(i, j) Sp[(i) + (size_t) nrow * (j)]
  #define D(i, j) Dp[(i) + (size_t) nrow * (j)]
  double xi = par[0], kappa = par[1], tau = par[2], zq = par[5], sig0 = par[6], sL = par[7], ptail = par[8], nuL = par[9], pi_miss = par[10];
  double delta[3] = {0.0, par[3], par[4]};
  bool rel = (jac == 1), rt = (rtau == 1), floor_on = pi_miss > 0;
  double l1mp = floor_on ? std::log1p(-pi_miss) : 0.0;
  double sig = sig0 * std::exp(kappa), reach = reach_mult * tau;
  Bulk B = make_bulk(sL, nuL);
  Mix M; bool mix = make_mix(bmix, M);
  bool has_miss = false, has_make = false;
  double wmade_min = INFINITY, wmiss_max = -INFINITY, wlast = -INFINITY;
  double wdec[3], wpdec[3]; int sdec[3]; int nj = 0;
  for (int j = 0; j < 3; ++j) {
    if (ISNAN(W(i, j))) continue;
    wlast = W(i, j);
    if (S(i, j) == 1) { has_make = true; if (W(i, j) < wmade_min) wmade_min = W(i, j); }
    else { has_miss = true; if (W(i, j) > wmiss_max) wmiss_max = W(i, j); }
    if (D(i, j) > 0) { wdec[nj] = W(i, j); wpdec[nj] = j > 0 ? W(i, j - 1) : 0.0; sdec[nj] = D(i, j); ++nj; }
  }
  double reach_lo = rt ? reach * wmade_min : reach, reach_hi = rt ? reach * (has_miss ? wmiss_max : wlast) : reach;
  double a = has_make ? wmade_min - reach_lo : (mix ? mu[i] + om[i] * M.lo : mu[i] - 7.0 * om[i] * sL);
  double b = has_miss ? wmiss_max + reach_hi : wlast + reach_hi;
  bool any_dec = false;
  for (int j = 0; j < nj; ++j) {
    double m = dec[sdec[j] - 1], sp = dec[3 + sdec[j] - 1], sm = dec_sm(dec, sdec[j]);
    double lo_h = rel ? wdec[j] * (1.0 + m - 8.0 * sm) : wdec[j] * std::exp(m - 8.0 * sm);
    double hi_h = rel ? wdec[j] * (1.0 + m + 8.0 * sp) : wdec[j] * std::exp(m + 8.0 * sp);
    if (!(trunc && sdec[j] >= 2) && lo_h > a) a = lo_h;   // as in unit_ll_ext_cpp
    if (!has_miss && !any_dec) b = hi_h; else if (hi_h < b) b = hi_h;
    any_dec = true;
  }
  if (b < a + 1e-3) b = a + 1e-3;
  double zs = mix ? mu[i] + om[i] * M.u : mu[i] + om[i] * zq;
  int K = 3 * G + 5 * 16 + 8;
  out.z.resize(K); out.v.resize(K); out.u.resize(K);
  std::vector<double> lwk(K);
  int nk = layout_nodes(a, b, zs, wdec, sdec, nj, dec, G, out.z.data(), lwk.data(), K, rel);
  out.nk = nk;
  for (int k = 0; k < nk; ++k) {
    double z = out.z[k];
    double e = (z - mu[i]) / om[i];
    double lf = (mix ? log_f_mix(e, M) : log_f_std(e, xi, sig, zq, ptail, B)) - std::log(om[i]);
    double la = 0.0;
    for (int j = 0; j < 3; ++j) {
      if (ISNAN(W(i, j))) continue;
      double eta = rt ? ((z > 0) ? std::log(z / W(i, j)) / tau : -INFINITY) + delta[j] : (z - W(i, j)) / tau + delta[j];
      if (!floor_on) la += (S(i, j) == 1) ? log_plogis(eta) : log_plogis(-eta);
      else la += (S(i, j) == 1) ? l1mp + log_plogis(eta) : std::log(pi_miss + (1.0 - pi_miss) * std::exp(log_plogis(-eta)));
    }
    double ld = 0.0;
    if (nj > 0) { DecSetup ds; dec_setup(z, wdec, wpdec, sdec, nj, dec, jac, ds); for (int j = 0; j < nj; ++j) ld += log_declaration_j(ds, j, trunc); }
    out.v[k] = lf + la + ld + lwk[k];
    double uu = 1.0 - (mix ? surv_mix(e, M) : surv_std(e, xi, sig, zq, ptail, B));
    if (uu < 1e-12) uu = 1e-12; if (uu > 1 - 1e-12) uu = 1 - 1e-12;
    out.u[k] = uu;
  }
  #undef W
  #undef S
  #undef D
}

// [[Rcpp::export]]
List joint_ll_ext_cpp(NumericMatrix Ws, IntegerMatrix Ss, IntegerMatrix Ds, NumericMatrix Wc, IntegerMatrix Sc, IntegerMatrix Dc,
                      NumericVector mus, NumericVector oms, NumericVector muc, NumericVector omc,
                      NumericVector ps, NumericVector pc, NumericVector decs, NumericVector decc,
                      int G, double reach_mult, int jac, int rtau, int trunc,
                      double rho, double psi, bool draw, int nthreads, NumericVector bmixs, NumericVector bmixc) {
  int n = Ws.nrow();
  NumericVector out(n), zs_draw(n), zc_draw(n);
  bool gauss = psi < 1e-8;
  double nu = gauss ? 1e9 : 1.0 / psi;
  double d = 1.0 - rho * rho;
  double lconst = gauss ? 0.0 : (R::lgammafn((nu + 2.0) / 2.0) + R::lgammafn(nu / 2.0) - 2.0 * R::lgammafn((nu + 1.0) / 2.0));
  const double* Wsp = Ws.begin(); const int* Ssp = Ss.begin(); const int* Dsp = Ds.begin();
  const double* Wcp = Wc.begin(); const int* Scp = Sc.begin(); const int* Dcp = Dc.begin();
  const double* musp = mus.begin(); const double* omsp = oms.begin(); const double* mucp = muc.begin(); const double* omcp = omc.begin();
  double* outp = out.begin(); double* zsp = zs_draw.begin(); double* zcp = zc_draw.begin();
  // draws use R's generator and must stay serial
  std::vector<double> ur;
  if (draw) { ur.resize(n); for (int i = 0; i < n; ++i) ur[i] = R::unif_rand(); }
#ifdef _OPENMP
#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (int i = 0; i < n; ++i) {
    LiftGridX gs, gc;
    lift_grid_ext(i, Wsp, Ssp, Dsp, n, musp, omsp, ps, decs, G, reach_mult, jac, rtau, trunc > 0, gs, bmixs);
    lift_grid_ext(i, Wcp, Scp, Dcp, n, mucp, omcp, pc, decc, G, reach_mult, jac, rtau, trunc > 0, gc, bmixc);
    int Ks = gs.nk, Kc = gc.nk;
    std::vector<double> xs(Ks), xc(Kc), cell((size_t) Ks * Kc);
    for (int k = 0; k < Ks; ++k) xs[k] = gauss ? R::qnorm(gs.u[k], 0.0, 1.0, 1, 0) : R::qt(gs.u[k], nu, 1, 0);
    for (int k = 0; k < Kc; ++k) xc[k] = gauss ? R::qnorm(gc.u[k], 0.0, 1.0, 1, 0) : R::qt(gc.u[k], nu, 1, 0);
    double mx = -INFINITY;
    for (int g = 0; g < Ks; ++g) {
      if (!std::isfinite(gs.v[g])) { for (int h = 0; h < Kc; ++h) cell[(size_t) g * Kc + h] = -INFINITY; continue; }
      for (int h = 0; h < Kc; ++h) {
        if (!std::isfinite(gc.v[h])) { cell[(size_t) g * Kc + h] = -INFINITY; continue; }
        double lc;
        if (gauss) { double x = xs[g], y = xc[h]; lc = -0.5 * std::log(d) - (rho * rho * (x * x + y * y) - 2.0 * rho * x * y) / (2.0 * d); }
        else lc = log_tcop(xs[g], xc[h], rho, nu, lconst, d);
        double v = gs.v[g] + gc.v[h] + lc;
        cell[(size_t) g * Kc + h] = v;
        if (v > mx) mx = v;
      }
    }
    double tot = 0.0;
    for (size_t k = 0; k < cell.size(); ++k) if (std::isfinite(cell[k])) tot += std::exp(cell[k] - mx);
    outp[i] = std::isfinite(mx) ? mx + std::log(tot) : -690.0;
    if (draw && std::isfinite(mx)) {
      double r = ur[i] * tot, cum = 0.0; size_t pick = cell.size() - 1;
      for (size_t k = 0; k < cell.size(); ++k) {
        if (!std::isfinite(cell[k])) continue;
        cum += std::exp(cell[k] - mx);
        if (cum >= r) { pick = k; break; }
      }
      zsp[i] = gs.z[pick / Kc]; zcp[i] = gc.z[pick % Kc];
    }
  }
  return List::create(Named("ll") = out, Named("zs") = zs_draw, Named("zc") = zc_draw);
}


// Posterior summaries of the latent capability of each unit under the
// extended marginal model with the athlete factor off: posterior mean,
// standard deviation, and probability that Z exceeds thr[i].
// [[Rcpp::export]]
NumericMatrix unit_latent_ext_cpp(NumericMatrix W, IntegerMatrix S, IntegerMatrix D,
                                  NumericVector mu, NumericVector om, NumericVector par, NumericVector dec,
                                  int G, double reach_mult, int jac, int rtau, int trunc, NumericVector thr, int nthreads, NumericVector bmix) {
  int n = W.nrow();
  NumericMatrix out(n, 3);
  const double* Wp = W.begin(); const int* Sp = S.begin(); const int* Dp = D.begin();
  const double* mup = mu.begin(); const double* omp = om.begin(); const double* thrp = thr.begin();
  double* op = out.begin();
#ifdef _OPENMP
#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (int i = 0; i < n; ++i) {
    LiftGridX g;
    lift_grid_ext(i, Wp, Sp, Dp, n, mup, omp, par, dec, G, reach_mult, jac, rtau, trunc > 0, g, bmix);
    double mx = -INFINITY;
    for (int k = 0; k < g.nk; ++k) if (std::isfinite(g.v[k]) && g.v[k] > mx) mx = g.v[k];
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, sa = 0.0;
    for (int k = 0; k < g.nk; ++k) {
      if (!std::isfinite(g.v[k])) continue;
      double w = std::exp(g.v[k] - mx);
      s0 += w; s1 += w * g.z[k]; s2 += w * g.z[k] * g.z[k]; if (g.z[k] > thrp[i]) sa += w;
    }
    double m = s0 > 0 ? s1 / s0 : NA_REAL;
    op[i] = m; op[i + n] = s0 > 0 ? std::sqrt(std::max(s2 / s0 - m * m, 0.0)) : NA_REAL; op[i + 2 * n] = s0 > 0 ? sa / s0 : NA_REAL;
  }
  return out;
}

// test hooks
// [[Rcpp::export]]
double log_eta_integral_cpp(double z, NumericVector w, NumericVector wprev, IntegerVector sit, NumericVector dec, double s_eta, int jac) {
  std::vector<double> ww(w.begin(), w.end()), wp(wprev.begin(), wprev.end()); std::vector<int> ss(sit.begin(), sit.end());
  DecSetup ds; dec_setup(z, ww.data(), wp.data(), ss.data(), w.size(), dec, jac, ds);
  return log_eta_integral(ds, s_eta);
}
// [[Rcpp::export]]
double log_eta_integral_pt_cpp(double z, NumericVector w, NumericVector wprev, IntegerVector sit, NumericVector dec,
                               double s_eta, NumericVector glx, NumericVector glw, int jac) {
  std::vector<double> ww(w.begin(), w.end()), wp(wprev.begin(), wprev.end()); std::vector<int> ss(sit.begin(), sit.end());
  DecSetup ds; dec_setup(z, ww.data(), wp.data(), ss.data(), w.size(), dec, jac, ds);
  return log_eta_integral_pt(ds, s_eta, glx.begin(), glw.begin(), glx.size());
}
// [[Rcpp::export]]
double log_declaration_cpp(double z, NumericVector w, NumericVector wprev, IntegerVector sit, NumericVector dec, int trunc, int jac) {
  std::vector<double> ww(w.begin(), w.end()), wp(wprev.begin(), wprev.end()); std::vector<int> ss(sit.begin(), sit.end());
  DecSetup ds; dec_setup(z, ww.data(), wp.data(), ss.data(), w.size(), dec, jac, ds);
  double v = 0.0; for (int j = 0; j < (int) w.size(); ++j) v += log_declaration_j(ds, j, trunc > 0);
  return v;
}
