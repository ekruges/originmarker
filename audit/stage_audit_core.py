"""Quantitative checks for the stage-inference claim chain."""
import numpy as np, json
from scipy import stats, optimize

H = 0.168
LADDER = {"bulk": (2e6, 0.013), "TE biopsy": (15, 0.050),
          "single ESC": (2, 0.199), "blastomere": (2, 0.308)}
out = {}

# ---- A. does an independent-per-template failure model reproduce the ladder? ----
# d = f^n  (every one of n template molecules of that allele must fail to prime)
A = {}
for anchor, (n_a, d_a) in LADDER.items():
    f = d_a ** (1.0 / n_a)
    preds = {k: min(1.0, f ** n) for k, (n, d) in LADDER.items()}
    A[anchor] = {"f_implied": f, "pred": preds,
                 "obs": {k: d for k, (n, d) in LADDER.items()}}
out["A_copy_number_model"] = A

# floor + sampling decomposition: d = e0 + (1-e0)*f^n , fit e0,f to the 4 rungs
def resid(p):
    e0, f = p
    return [e0 + (1 - e0) * f ** n - d for k, (n, d) in LADDER.items()]
fit = optimize.least_squares(resid, [0.01, 0.4], bounds=([0, 0], [1, 1]))
e0, f = fit.x
out["A_floor_fit"] = {"e0": e0, "f": f,
                      "pred": {k: e0 + (1 - e0) * f ** n for k, (n, d) in LADDER.items()},
                      "resid_rms": float(np.sqrt(np.mean(np.array(resid(fit.x)) ** 2)))}

# ---- B. confound bias in d_naive = 1 - h/H ----
# true forward model: h = H_true*(1-F)*(1-L)*(1-d) + drop-in term (set 0 here)
def d_naive(h, Hv=H): return 1 - h / Hv
B = []
for label, Htrue, F, L in [("clean", H, 0, 0),
                           ("2nd cousins F=0.0156", H, 0.0156, 0),
                           ("1st cousins F=0.0625", H, 0.0625, 0),
                           ("uncle-niece F=0.125", H, 0.125, 0),
                           ("5% genome LOH/UPD", H, 0, 0.05),
                           ("20% genome LOH", H, 0, 0.20),
                           ("whole-chr14 UPD (2.9% of markers)", H, 0, 0.029),
                           ("anchor 5% too high", H * 0.95, 0, 0),
                           ("anchor 10% too low", H * 1.10, 0, 0)]:
    row = {"scenario": label}
    for st, (n, d) in LADDER.items():
        h = Htrue * (1 - F) * (1 - L) * (1 - d)
        row[st] = {"h": h, "d_naive": d_naive(h), "bias": d_naive(h) - d}
    B.append(row)
out["B_confounds"] = B

# ---- C. replicate-discordance estimator (no anchor needed) ----
# two independent amplifications of one genome; truly-het marker drops with prob d each.
# among markers het in >=1 replicate: discordant fraction phi = 2d(1-d)/(1-d^2) = 2d/(1+d)
# invert: d = phi/(2-phi).  correlated dropout (shared hard failures, corr rho) biases low.
C = {}
for d in [0.013, 0.05, 0.199, 0.308]:
    phi = 2 * d / (1 + d)
    C[str(d)] = {"phi_expected": phi, "d_recovered": phi / (2 - phi)}
# effect of correlated dropout: P(both drop) = d^2 + rho*d*(1-d)
Ccorr = {}
for rho in [0.0, 0.1, 0.25, 0.5]:
    r = {}
    for d in [0.05, 0.199, 0.308]:
        p_both = d * d + rho * d * (1 - d)
        p_one = 2 * (d - p_both)
        phi = p_one / (p_one + (1 - 2 * d + p_both))
        r[str(d)] = {"phi": phi, "d_hat": phi / (2 - phi), "rel_underest": (phi / (2 - phi)) / d}
    Ccorr[str(rho)] = r
out["C_replicate_estimator"] = {"ideal": C, "correlated": Ccorr}

# precision of that estimator vs number of informative (het-in-either) markers
Cprec = {}
for d in [0.05, 0.199, 0.308]:
    phi = 2 * d / (1 + d)
    for N in [100, 200, 1000, 10000, 100000]:
        se_phi = np.sqrt(phi * (1 - phi) / N)
        se_d = se_phi * 2 / (2 - phi) ** 2  # delta method
        Cprec[f"d={d},N={N}"] = {"se_d": se_d}
out["C_precision"] = Cprec

# ---- D. haploid/diploid boundary at h=0.105, with sampling noise ----
BOUND = 0.105
D = []
truth = {"haploid, drop-in only": 0.045,
         "haploid, high end observed": 0.100,
         "blastomere d=0.308, no drop-in": H * (1 - 0.308),
         "blastomere d=0.308, F=0.0625": H * (1 - 0.308) * (1 - 0.0625),
         "blastomere d=0.40": H * (1 - 0.40),
         "single ESC d=0.199": H * (1 - 0.199)}
for lab, h in truth.items():
    row = {"case": lab, "h_true": h}
    for N in [100, 200, 1000, 10000, 825657]:
        se = np.sqrt(h * (1 - h) / N)
        # P(called haploid) = P(h_hat < BOUND); exact binomial
        k = int(np.floor(BOUND * N))
        p_hap = stats.binom.cdf(k, N, h) if h > 0 else 1.0
        row[f"N={N}"] = {"se": se, "P_called_haploid": float(p_hap)}
    D.append(row)
out["D_boundary"] = D

# marker count needed to keep both error rates < 1% for the tightest real pair
def n_needed(h0, h1, alpha=0.01, beta=0.01):
    z0, z1 = stats.norm.isf(alpha), stats.norm.isf(beta)
    s0, s1 = np.sqrt(h0 * (1 - h0)), np.sqrt(h1 * (1 - h1))
    return ((z0 * s0 + z1 * s1) / abs(h1 - h0)) ** 2
out["D_n_needed"] = {
    "haploid0.100_vs_blastomere0.116": n_needed(0.100, H * (1 - 0.308)),
    "haploid0.045_vs_blastomere0.116": n_needed(0.045, H * (1 - 0.308)),
    "bulk_vs_TE": n_needed(H * (1 - 0.013), H * (1 - 0.050)),
    "TE_vs_singleESC": n_needed(H * (1 - 0.050), H * (1 - 0.199)),
    "singleESC_vs_blastomere": n_needed(H * (1 - 0.199), H * (1 - 0.308)),
}

# ---- E. can the stated drop-in coexist with the ladder in one forward model? ----
E = {}
for i in [0.0, 0.0390, 0.0525]:
    E[str(i)] = {k: {"h": H * (1 - d) + (1 - H) * i,
                     "d_naive": d_naive(H * (1 - d) + (1 - H) * i)}
                 for k, (n, d) in LADDER.items()}
out["E_dropin_consistency"] = E
# what per-marker spurious-het rate is compatible with a haploid at h=0.002?
out["E_haploid_floor"] = {"h=0.002 implies i": 0.002, "h=0.100 implies i": 0.100,
                          "measured dropin_range": [0.0390, 0.0525]}

json.dump(out, open("audit_core.json", "w"), indent=1, default=float)

# ---- compact console report ----
print("A. independent-template model d=f^n, calibrated on each rung:")
for a, v in A.items():
    print(f"  anchored {a:<11} f={v['f_implied']:.4g}  pred " +
          " ".join(f"{k}={v['pred'][k]:.3g}" for k in LADDER))
print(f"  obs      {'':<11}         " + " ".join(f"{k}={LADDER[k][1]:.3g}" for k in LADDER))
print(f"  floor+sampling fit: e0={e0:.4f} f={f:.4f} rms={out['A_floor_fit']['resid_rms']:.4f}")
print("  floor-fit pred:", {k: round(v, 4) for k, v in out["A_floor_fit"]["pred"].items()})
print()
print("B. bias of 1-h/H (bias in dropout units), by confound:")
for row in B:
    print(f"  {row['scenario']:<34} " + " ".join(f"{st.split()[0]}{row[st]['bias']:+.3f}" for st in LADDER))
print()
print("C. replicate-discordance estimator: phi -> d")
for d, v in C.items():
    print(f"  d={d}: expected discordance among het-in-either = {v['phi_expected']:.4f} -> recovers {v['d_recovered']:.4f}")
print("  with correlated dropout (rho): d_hat/d_true")
for rho, r in Ccorr.items():
    print(f"   rho={rho}: " + " ".join(f"d={d}->{v['rel_underest']:.2f}x" for d, v in r.items()))
print()
print("D. P(misclassified as haploid) at boundary h<0.105:")
for row in D:
    print(f"  {row['case']:<34} h={row['h_true']:.4f} " +
          " ".join(f"N={N}:{row[f'N={N}']['P_called_haploid']:.3f}" for N in [100, 200, 1000, 10000]))
print()
print("D2. markers needed for 1%/1% error between adjacent classes:")
for k, v in out["D_n_needed"].items():
    print(f"  {k:<38} N>={v:,.0f}")
print()
print("E. forward model with stated drop-in i:")
for i, v in E.items():
    print(f"  i={i}: " + " ".join(f"{k}:h={v[k]['h']:.4f},d_naive={v[k]['d_naive']:+.3f}" for k in LADDER))
