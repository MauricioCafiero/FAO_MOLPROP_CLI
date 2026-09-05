"""
Statistical comparison of zero-shot / few-shot / frag-shot single-call baselines,
the agentic 5x4 loop (5 proposers common to all four), and two non-LLM GA
baselines (results/batches/ga_baseline/5x4 and 5x4_full, 5 replicates each,
no proposer identity).

Two GA baselines, not one: the GA's substituent pool must match whichever
LLM condition it's being held up against, or a score gap could just reflect
one condition having a bigger chemical space to search, not a better search
strategy (see GA_BASELINE_EXPLAINER.md "Two substituent pools: frag10 vs
full"). `ga_frag10` (5x4) is restricted to the exact 10 fragments frag-shot's
prompt showed the LLM -- the intended comparison is ga_frag10 vs frag. `ga_full`
(5x4_full) searches the full ~390-item combined pool -- the intended
comparison is ga_full vs zero, since zero-shot's prompt shows the LLM no
fragment menu at all. Both are still included in every pooled/omnibus test
below for completeness, but only those two pairings are pool-matched by
design; other GA pairwise rows (e.g. ga_frag10 vs zero, ga_full vs frag) are
reported but weren't designed for a fair chemical-space comparison.

Unit-of-analysis discipline:
  - compound-level rows are NOT independent (multiple compounds share a replicate,
    multiple replicates share a proposer) -> naive compound-level tests pseudoreplicate.
  - Primary tests use REPLICATE-level aggregates (n=5 proposers x ~5 reps per condition)
    with proposer as a blocking/random-effect variable.
  - Compound-level tests are reported as a secondary/liberal lens only, explicitly
    flagged as pseudoreplicated.
  - The GA baselines have no proposer dimension (each is one system, not five), so
    they're included in every condition-pooled test (KW, pairwise MW-U, good-lead
    GEE, compound-level KW) but excluded from the proposer-BLOCKED tests (Friedman,
    paired Wilcoxon, the proposer-random-intercept mixed model) -- those require
    the same 5 units observed under every condition, which GA doesn't satisfy.

Outputs a single Markdown report to stdout (redirect to a file).
"""
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Exchangeable
import itertools
import warnings
warnings.filterwarnings("ignore")

ROOT = "results/batches"

GOOD_LEAD_THRESH = -8.0

# condition -> (csv path, set_label filter or None, proposer col name to use)
SOURCES = {
    "zero":  f"{ROOT}/zero_shot/analysis/compounds_zero_shot.csv",
    "few":   f"{ROOT}/few_shot/analysis/compounds_few_shot.csv",
    "frag":  f"{ROOT}/frag_shot/analysis/compounds_frag_shot.csv",
}

GA_SOURCES = {
    "ga_frag10": f"{ROOT}/ga_baseline/5x4/analysis/compounds_ga_5x4.csv",
    "ga_full":   f"{ROOT}/ga_baseline/5x4_full/analysis/compounds_ga_5x4_full.csv",
}

AGENTIC = {
    "openai":     f"{ROOT}/openai_gpt-5.2_vs_anthropic_5x4/analysis/compounds_e2e_test.csv",
    "anthropic":  f"{ROOT}/anthropic-haiku-4-5_vs_openai_5x4/analysis/compounds_e2e_test.csv",
    "gemini":     f"{ROOT}/gemini-3-flash-preview_vs_openai_5x4/analysis_full5/compounds_gemini-3-flash-preview_vs_openai.csv",
    "kimi":       f"{ROOT}/ollama_kimi-k2.6_vs_openai_5x4/analysis/compounds_ollama_kimi-k2.6_vs_openai.csv",
    "deepseek":   f"{ROOT}/ollama_deepseek-v4-pro_vs_openai_5x4/analysis_full5/compounds_ollama_deepseek-v4-pro_vs_openai.csv",
}

# original-5 set_label -> canonical proposer key, for the zero/few/frag files
BASELINE_LABEL_MAP = {
    "openai": "openai",
    "anthropic": "anthropic",
    "gemini": "gemini",
    "kimi-k2.6": "kimi",
    "deepseek-v4-pro": "deepseek",
}

# First-response baseline: each agentic batch's turn-1 compounds (proposer with
# tools, no adversary, no iteration), mined from the sidecar transcripts by
# code/make_first_response_baseline.py and re-scored with the real Vina engine.
# set_label is already the canonical proposer key, and the proposer dimension is
# intact (5 proposers x 5 replicates), so unlike the GAs this condition joins the
# proposer-BLOCKED tests (Friedman, paired Wilcoxon, mixed model).
FIRST_RESPONSE = f"{ROOT}/first_response_5x4/analysis/compounds_first_response_5x4.csv"

def load_first_response():
    df = pd.read_csv(FIRST_RESPONSE)
    df["proposer"] = df["set_label"]
    df["condition"] = "first"
    return df[["proposer", "condition", "replicate", "docking", "docked_in_pocket"]]

def load_baseline(condition):
    df = pd.read_csv(SOURCES[condition])
    df = df[df["set_label"].isin(BASELINE_LABEL_MAP.keys())].copy()
    df["proposer"] = df["set_label"].map(BASELINE_LABEL_MAP)
    df["condition"] = condition
    return df[["proposer", "condition", "replicate", "docking", "docked_in_pocket"]]

def load_agentic():
    frames = []
    for proposer, path in AGENTIC.items():
        df = pd.read_csv(path)
        df["proposer"] = proposer
        df["condition"] = "agentic"
        frames.append(df[["proposer", "condition", "replicate", "docking", "docked_in_pocket"]])
    return pd.concat(frames, ignore_index=True)

def load_ga(condition):
    df = pd.read_csv(GA_SOURCES[condition])
    df["proposer"] = condition
    df["condition"] = condition
    return df[["proposer", "condition", "replicate", "docking", "docked_in_pocket"]]

def good_lead(row):
    return bool(row["docked_in_pocket"]) and row["docking"] <= GOOD_LEAD_THRESH

def main():
    baseline = pd.concat([load_baseline(c) for c in SOURCES], ignore_index=True)
    agentic = load_agentic()
    ga_frag10 = load_ga("ga_frag10")
    ga_full = load_ga("ga_full")
    df = pd.concat([baseline, agentic, ga_frag10, ga_full, load_first_response()], ignore_index=True)
    n_before = len(df)
    df = df[df["docking"].notna()].copy()
    if len(df) != n_before:
        print(f"(dropped {n_before - len(df)} compound row(s) with missing docking score before analysis)\n")
    df["docked_in_pocket"] = df["docked_in_pocket"].astype(str).str.lower().isin(["true", "1"])
    df["good_lead"] = df.apply(good_lead, axis=1).astype(int)
    df["condition"] = pd.Categorical(
        df["condition"], categories=["zero", "frag", "few", "first", "agentic", "ga_frag10", "ga_full"], ordered=True)

    # ---------- replicate-level aggregation ----------
    rep = df.groupby(["proposer", "condition", "replicate"], observed=True).agg(
        n_compounds=("docking", "size"),
        mean_docking=("docking", "mean"),
        best_docking=("docking", "min"),
        good_lead_rate=("good_lead", "mean"),
    ).reset_index()

    print("# Statistical comparison: zero / frag / few-shot vs. agentic 5x4 vs. GA baselines\n")
    print(f"n replicates per condition (pooled over 5 proposers): "
          f"{rep.groupby('condition', observed=True).size().to_dict()}\n")

    print("## 1. Primary test — replicate-level mean docking score, all 7 conditions pooled\n")
    print("Kruskal-Wallis across all 7 conditions (unit = one replicate's mean docking score, "
          f"n={len(rep)} replicates; zero/frag/few/first/agentic pooled over 5 proposers, ga_frag10/ga_full "
          "are each their own 5 replicates with no proposer dimension):\n")
    groups = [g["mean_docking"].values for _, g in rep.groupby("condition", observed=True)]
    h, p = stats.kruskal(*groups)
    print(f"- H = {h:.3f}, p = {p:.2e}\n")

    print("### Pairwise Mann-Whitney U (Holm-corrected), replicate mean docking\n")
    print("Rows in **bold** are the two pool-matched GA comparisons this baseline was designed for "
          "(ga_frag10 vs frag, ga_full vs zero); other GA rows are reported for completeness but the "
          "GA's chemical space wasn't matched to those specific conditions.\n")
    print("| Comparison | n1 | n2 | U | raw p | Holm-adj p | rank-biserial r (effect size) |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    conds = ["zero", "frag", "few", "first", "agentic", "ga_frag10", "ga_full"]
    pairs = list(itertools.combinations(conds, 2))
    raw_ps = []
    stats_cache = []
    for a, b in pairs:
        xa = rep[rep.condition == a]["mean_docking"].values
        xb = rep[rep.condition == b]["mean_docking"].values
        u, p_ = stats.mannwhitneyu(xa, xb, alternative="two-sided")
        r = 1 - (2 * u) / (len(xa) * len(xb))  # rank-biserial
        raw_ps.append(p_)
        stats_cache.append((a, b, len(xa), len(xb), u, p_, r))
    # Holm-Bonferroni
    order = np.argsort(raw_ps)
    m = len(raw_ps)
    adj_ps = [None] * m
    running_max = 0
    for rank, idx in enumerate(order):
        adj = min(1.0, raw_ps[idx] * (m - rank))
        running_max = max(running_max, adj)
        adj_ps[idx] = running_max
    pool_matched = {frozenset(("frag", "ga_frag10")), frozenset(("zero", "ga_full"))}
    for (a, b, n1, n2, u, p_, r), adj in zip(stats_cache, adj_ps):
        label = f"**{a} vs {b}**" if frozenset((a, b)) in pool_matched else f"{a} vs {b}"
        print(f"| {label} | {n1} | {n2} | {u:.1f} | {p_:.2e} | {adj:.2e} | {r:+.2f} |")
    print()

    print("## 2. Design-matched test — Friedman (proposer as block, n=5 proposers)\n")
    print("Restricted to zero/frag/few/first/agentic -- the GA baselines have no proposer dimension "
          "(each is one system run 5 times, not 5 proposers each run once), so neither can be a block "
          "in a proposer-matched design and both are excluded from this test and §3's mixed model. "
          "They're already covered by §1's pooled tests above. first-response shares the agentic "
          "proposers (its compounds are mined from the same transcripts' turn 1), so it joins the "
          "blocked tests.\n")
    print("Uses each proposer's mean-of-replicate-means per condition, so proposer identity "
          "(the strongest confound: proposers differ far more than conditions within a proposer) "
          "is controlled for by construction.\n")
    conds5 = ["zero", "frag", "few", "first", "agentic"]
    prop_means = rep.groupby(["proposer", "condition"], observed=True)["mean_docking"].mean().unstack("condition")
    prop_means = prop_means.loc[["openai", "anthropic", "gemini", "kimi", "deepseek"], conds5]
    print(prop_means.round(3).to_markdown())
    print()
    fr_stat, fr_p = stats.friedmanchisquare(*[prop_means[c].values for c in conds5])
    print(f"Friedman chi2 = {fr_stat:.3f}, p = {fr_p:.4f} (df=4, n=5 blocks — likely underpowered; "
          "treat as indicative, not confirmatory)\n")

    print("### Paired Wilcoxon signed-rank, agentic vs each other LLM condition (n=5 proposers)\n")
    print("| Comparison | W | p (two-sided, exact) | median diff (other − agentic) |")
    print("|---|---:|---:|---:|")
    for c in ["zero", "frag", "few", "first"]:
        d = prop_means[c] - prop_means["agentic"]
        w, p_ = stats.wilcoxon(d, alternative="two-sided", mode="exact")
        print(f"| {c} vs agentic | {w:.1f} | {p_:.4f} | {d.median():+.3f} |")
    print()

    print("## 3. Linear mixed-effects model — replicate-level, proposer as random intercept\n")
    print("Restricted to zero/frag/few/first/agentic for the same reason as §2 (GA baselines have no "
          "proposer to supply a random intercept for; first-response shares the agentic proposers).\n")
    print("`mean_docking ~ C(condition, Treatment('agentic')) + (1 | proposer)`, "
          "unit = replicate (n=25/condition for zero/few/agentic, n=24 for frag [deepseek frag-shot rep3 extraction miss]). "
          "Coefficients are kcal/mol relative to agentic; negative-going coefficient = single-shot condition scored WORSE "
          "(less negative docking) than agentic.\n")
    model_df = rep[~rep.condition.isin(["ga_frag10", "ga_full"])].copy()
    model_df["condition"] = model_df["condition"].astype(str)
    md = smf.mixedlm("mean_docking ~ C(condition, Treatment('agentic'))", model_df, groups=model_df["proposer"])
    fit = md.fit(reml=True)
    print("```")
    print(fit.summary())
    print("```\n")

    print("## 4. Good-lead rate (binary, compound-level) — GEE with exchangeable within-replicate correlation\n")
    print("`good_lead ~ condition`, family=Binomial, clustered by (proposer, replicate) to respect "
          "non-independence of compounds proposed in the same call. Odds ratios relative to agentic. "
          "Both GA baselines included as separate levels -- this test doesn't need a proposer dimension.\n")
    gee_df = df.copy()
    gee_df["cluster"] = gee_df["proposer"].astype(str) + "_" + gee_df["condition"].astype(str) + "_" + gee_df["replicate"].astype(str)
    gee_df["condition"] = pd.Categorical(
        gee_df["condition"], categories=["agentic", "zero", "frag", "few", "first", "ga_frag10", "ga_full"])
    fam = Binomial()
    ind = Exchangeable()
    gee_model = GEE.from_formula("good_lead ~ C(condition, Treatment('agentic'))", groups="cluster", data=gee_df, family=fam, cov_struct=ind)
    gee_fit = gee_model.fit()
    print("```")
    print(gee_fit.summary())
    print("```\n")

    # flag any condition with 0 or 100% good-lead compounds at the compound level:
    # that produces quasi-complete separation, and the corresponding coefficient/CI
    # above is a numerical artifact (drifts to +-inf), not a meaningful effect size.
    lead_counts = gee_df.groupby("condition", observed=True)["good_lead"].agg(["sum", "count"])
    degenerate = lead_counts[(lead_counts["sum"] == 0) | (lead_counts["sum"] == lead_counts["count"])]
    if len(degenerate):
        print("**Note — quasi-complete separation:** the following condition(s) have 0 or 100% "
              "good-lead compounds, which drives their GEE coefficient toward +-infinity above; "
              "treat those specific rows as \"no good leads observed\" / \"all good leads\", not as a "
              "precise odds ratio:\n")
        for cond, row in degenerate.iterrows():
            print(f"- {cond}: {int(row['sum'])}/{int(row['count'])} good-lead compounds")
        print()

    print("Odds ratios (exp(coef)):\n")
    for name, coef in gee_fit.params.items():
        if name == "Intercept":
            continue
        print(f"- {name}: OR = {np.exp(coef):.3f} (95% CI {np.exp(coef - 1.96*gee_fit.bse[name]):.3f}"
              f" - {np.exp(coef + 1.96*gee_fit.bse[name]):.3f})")
    print()

    print("## 5. Secondary / liberal lens — compound-level Kruskal-Wallis (PSEUDOREPLICATED, flagged)\n")
    print(f"n compounds: {df.groupby('condition', observed=True).size().to_dict()}\n")
    print("This ignores clustering (many compounds per replicate/proposer) and will overstate "
          "significance — included only because it's the comparison most papers report by default, "
          "so it's here to be explicitly contrasted with the corrected tests above.\n")
    groups_c = [g["docking"].values for _, g in df.groupby("condition", observed=True)]
    h_c, p_c = stats.kruskal(*groups_c)
    print(f"- H = {h_c:.3f}, p = {p_c:.2e}\n")

    print("## Artifacts\n")
    print("- Replicate-level table used above saved to `results/batches/analysis/shot_vs_agentic_replicate_level.csv`")
    rep.to_csv(f"{ROOT}/analysis/shot_vs_agentic_replicate_level.csv", index=False)

if __name__ == "__main__":
    main()
