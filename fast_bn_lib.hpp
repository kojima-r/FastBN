#pragma once

#include <cmath>

// ===== chi-square 上側確率 Q(k/2, x/2) を計算するための補助 =====
static inline double gammaln(double x) { return std::lgamma(x); }

// 正規化下側不完全ガンマ P(s,x)（級数展開）と
// 正規化上側不完全ガンマ Q(s,x)（連分数）を使い分け
static double regularized_gamma_P(double s, double x) {
    if (x <= 0.0) return 0.0;
    const int MAXIT = 1000;
    const double EPS = 1e-14;

    // 級数: P(s,x) = e^{-x} x^s / Γ(s) * sum_{n=0..∞} (x^n / Γ(s+n+1))
    double sum = 1.0 / s;
    double term = sum;
    for (int n=1; n<MAXIT; ++n) {
        term *= x / (s + n);
        sum += term;
        if (std::fabs(term) < std::fabs(sum) * EPS) break;
    }
    return std::exp(s*std::log(x) - x - gammaln(s)) * sum;
}

static double regularized_gamma_Q(double s, double x) {
    if (x <= 0.0) return 1.0;
    const int MAXIT = 1000;
    const double EPS = 1e-14;

    // 連分数: Q(s,x) を直接計算（Lentz 法）
    double C = 1.0 / 1e-300;
    double D = 0.0;
    double f = 0.0;

    // 初期化
    double a0 = 0.0;
    double b0 = x + 1.0 - s;
    D = 1.0 / std::max(1e-300, b0);
    C = std::max(1e-300, b0);
    f = D;

    for (int i=1; i<MAXIT; ++i) {
        double a = i * (s - i);
        double b = b0 + 2.0 * i;

        // D = 1 / (b + a*D)
        D = 1.0 / std::max(1e-300, b + a*D);
        C = std::max(1e-300, b + a/C);
        double delta = C * D;
        f *= delta;
        if (std::fabs(delta - 1.0) < EPS) break;
    }

    double pref = std::exp(s*std::log(x) - x - gammaln(s));
    return pref * f;
}

// P と Q の切替（数値安定化）
static double gamma_reg_P(double s, double x) {
    if (x < s + 1.0) return regularized_gamma_P(s, x);
    return 1.0 - regularized_gamma_Q(s, x);
}
static double gamma_reg_Q(double s, double x) {
    if (x < s + 1.0) return 1.0 - regularized_gamma_P(s, x);
    return regularized_gamma_Q(s, x);
}

// カイ二乗の上側確率（p 値）: df=k, statistic=chi2
static double chisq_p_upper(double chi2, int df) {
    if (chi2 < 0.0) return 1.0;
    double s = 0.5 * df;
    double x = 0.5 * chi2;
    return gamma_reg_Q(s, x); // p = Q(df/2, chi2/2)
}


