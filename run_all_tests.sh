#!/bin/sh
# run_all_tests.sh — Compile, run, cross-validate all algorithm implementations.
#
# Usage:  ./run_all_tests.sh            (run from repo root)
#         ./run_all_tests.sh --quick    (skip large graphs)
#
# Produces:  results/report.md          (markdown summary table)
#            results/raw/               (raw output from each run)
#
# Sequential execution — no parallelism, clean timing.

set -e

# ── locate repo root ─────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$SCRIPT_DIR"
ALGO="$REPO/algorithms"
DATA="$REPO/data"
RESULTS="$REPO/results"

QUICK=0
if [ "$1" = "--quick" ]; then
    QUICK=1
    echo "Quick mode: skipping large graphs."
    echo ""
fi

# ── timeout (macOS coreutils may not have timeout) ───────────────────────
if command -v gtimeout >/dev/null 2>&1; then
    TIMEOUT_CMD="gtimeout"
elif command -v timeout >/dev/null 2>&1; then
    TIMEOUT_CMD="timeout"
else
    TIMEOUT_CMD=""
fi

run_with_timeout() {
    limit="$1"; shift
    if [ -n "$TIMEOUT_CMD" ]; then
        "$TIMEOUT_CMD" "$limit" "$@" 2>&1
    else
        "$@" 2>&1
    fi
}

# ── setup results directory ──────────────────────────────────────────────
rm -rf "$RESULTS"
mkdir -p "$RESULTS/raw"

# ── general matching algorithms ──────────────────────────────────────────
GENERAL_ALGOS="edmonds-blossom-simple edmonds-blossom-optimized gabow-simple gabow-optimized micali-vazirani"
MV_PURE="micali-vazirani-pure"
BIPARTITE_ALGOS="hopcroft-karp"
LANGS="cpp rust python"

# derive source filename from algorithm directory name
src_name() {
    echo "$1" | tr '-' '_'
}

# ── Phase 1: Compile ─────────────────────────────────────────────────────
echo "============================================"
echo "  Phase 1: Compiling all implementations"
echo "============================================"
echo ""

compile_errors=0

for alg in $GENERAL_ALGOS $MV_PURE $BIPARTITE_ALGOS; do
    base="$(src_name "$alg")"
    alg_dir="$ALGO/$alg"

    # C++
    src="$alg_dir/cpp/${base}.cpp"
    bin="$alg_dir/cpp/${base}_cpp"
    if [ -f "$src" ]; then
        printf "  compile  %-45s" "$alg/cpp"
        if g++ -O3 -std=c++17 "$src" -o "$bin" 2>"$RESULTS/raw/${base}_cpp_compile.log"; then
            echo "✓"
        else
            echo "✗  (see results/raw/${base}_cpp_compile.log)"
            compile_errors=$((compile_errors + 1))
        fi
    fi

    # Rust
    src="$alg_dir/rust/${base}.rs"
    bin="$alg_dir/rust/${base}_rust"
    if [ -f "$src" ]; then
        printf "  compile  %-45s" "$alg/rust"
        if rustc -O "$src" -o "$bin" 2>"$RESULTS/raw/${base}_rust_compile.log"; then
            echo "✓"
        else
            echo "✗  (see results/raw/${base}_rust_compile.log)"
            compile_errors=$((compile_errors + 1))
        fi
    fi
done

echo ""
if [ "$compile_errors" -gt 0 ]; then
    echo "WARNING: $compile_errors compilation(s) failed. Continuing with what compiled."
    echo ""
fi

# ── Phase 2: Build test graph list ───────────────────────────────────────
GENERAL_GRAPHS=""
for size in small medium large; do
    dir="$DATA/general-unweighted/$size"
    [ -d "$dir" ] || continue
    if [ "$QUICK" -eq 1 ] && [ "$size" = "large" ]; then continue; fi
    for f in "$dir"/*.txt; do
        [ -f "$f" ] || continue
        GENERAL_GRAPHS="$GENERAL_GRAPHS $f"
    done
done

BIPARTITE_GRAPHS=""
for size in small medium large; do
    dir="$DATA/bipartite-unweigthed/$size"
    [ -d "$dir" ] || continue
    if [ "$QUICK" -eq 1 ] && [ "$size" = "large" ]; then continue; fi
    for f in "$dir"/*.txt; do
        [ -f "$f" ] || continue
        BIPARTITE_GRAPHS="$BIPARTITE_GRAPHS $f"
    done
done

# ── Phase 3: Run all tests ───────────────────────────────────────────────
echo "============================================"
echo "  Phase 2: Running all tests (sequential)"
echo "============================================"
echo ""

# CSV accumulator: algo,graph,lang,size,time_ms,valid,status
CSV="$RESULTS/raw/all_results.csv"
echo "algo,graph,lang,size,time_ms,valid,status" > "$CSV"

run_one() {
    alg="$1"       # e.g. edmonds-blossom-simple
    graph="$2"     # full path to .txt
    lang="$3"      # cpp|rust|python

    base="$(src_name "$alg")"
    gname="$(basename "$graph" .txt)"
    logfile="$RESULTS/raw/${base}_${lang}_${gname}.log"

    case "$lang" in
        cpp)
            bin="$ALGO/$alg/cpp/${base}_cpp"
            [ -x "$bin" ] || { echo "skip" > "$logfile"; return; }
            run_with_timeout 300 "$bin" "$graph" > "$logfile"
            ;;
        rust)
            bin="$ALGO/$alg/rust/${base}_rust"
            [ -x "$bin" ] || { echo "skip" > "$logfile"; return; }
            run_with_timeout 300 "$bin" "$graph" > "$logfile"
            ;;
        python)
            src="$ALGO/$alg/python/${base}.py"
            [ -f "$src" ] || { echo "skip" > "$logfile"; return; }
            run_with_timeout 300 uv run "$src" "$graph" > "$logfile"
            ;;
    esac

    # parse output
    size=$(grep "^Matching size:" "$logfile" | tail -1 | awk '{print $3}')
    tms=$(grep "^Time:" "$logfile" | awk '{print $2}')
    valid=$(grep "VALIDATION" "$logfile" | head -1)

    case "$valid" in
        *PASSED*) vflag="PASS" ;;
        *FAILED*) vflag="FAIL" ;;
        *)        vflag="NONE" ;;
    esac

    [ -z "$size" ] && size="ERR"
    [ -z "$tms" ] && tms="ERR"

    echo "$alg,$gname,$lang,$size,$tms,$vflag" >> "$CSV"
}

total=0
for alg in $GENERAL_ALGOS $MV_PURE; do
    for graph in $GENERAL_GRAPHS; do
        gname="$(basename "$graph" .txt)"
        for lang in $LANGS; do
            total=$((total + 1))
            printf "  [%3d] %-35s %-6s %-30s" "$total" "$alg" "$lang" "$gname"
            run_one "$alg" "$graph" "$lang"
            # re-read last csv line for display
            last="$(tail -1 "$CSV")"
            sz="$(echo "$last" | cut -d, -f4)"
            tm="$(echo "$last" | cut -d, -f5)"
            vf="$(echo "$last" | cut -d, -f6)"
            printf "  size=%-6s %4s ms  %s\n" "$sz" "$tm" "$vf"
        done
    done
done

for alg in $BIPARTITE_ALGOS; do
    for graph in $BIPARTITE_GRAPHS; do
        gname="$(basename "$graph" .txt)"
        for lang in $LANGS; do
            total=$((total + 1))
            printf "  [%3d] %-35s %-6s %-30s" "$total" "$alg" "$lang" "$gname"
            run_one "$alg" "$graph" "$lang"
            last="$(tail -1 "$CSV")"
            sz="$(echo "$last" | cut -d, -f4)"
            tm="$(echo "$last" | cut -d, -f5)"
            vf="$(echo "$last" | cut -d, -f6)"
            printf "  size=%-6s %4s ms  %s\n" "$sz" "$tm" "$vf"
        done
    done
done

echo ""
echo "Total runs: $total"

# ── Phase 4: Cross-validation ────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Phase 3: Cross-validation"
echo "============================================"
echo ""

cross_errors=0
cross_ok=0

# For each graph, all algorithms that ran on it must agree on matching size
for graph in $GENERAL_GRAPHS; do
    gname="$(basename "$graph" .txt)"
    sizes="$(grep ",$gname," "$CSV" | grep -v "^algo" | cut -d, -f4 | grep -v ERR | sort -u)"
    count="$(echo "$sizes" | wc -l | tr -d ' ')"
    if [ "$count" -eq 1 ] && [ -n "$sizes" ]; then
        printf "  %-50s  size=%-6s ✓\n" "$gname" "$sizes"
        cross_ok=$((cross_ok + 1))
    elif [ -z "$sizes" ]; then
        printf "  %-50s  NO DATA\n" "$gname"
    else
        printf "  %-50s  MISMATCH: %s ✗\n" "$gname" "$(echo $sizes | tr '\n' ' ')"
        cross_errors=$((cross_errors + 1))
    fi
done

for graph in $BIPARTITE_GRAPHS; do
    gname="$(basename "$graph" .txt)"
    sizes="$(grep ",$gname," "$CSV" | grep -v "^algo" | cut -d, -f4 | grep -v ERR | sort -u)"
    count="$(echo "$sizes" | wc -l | tr -d ' ')"
    if [ "$count" -eq 1 ] && [ -n "$sizes" ]; then
        printf "  %-50s  size=%-6s ✓\n" "$gname" "$sizes"
        cross_ok=$((cross_ok + 1))
    elif [ -z "$sizes" ]; then
        printf "  %-50s  NO DATA\n" "$gname"
    else
        printf "  %-50s  MISMATCH: %s ✗\n" "$gname" "$(echo $sizes | tr '\n' ' ')"
        cross_errors=$((cross_errors + 1))
    fi
done

echo ""
if [ "$cross_errors" -eq 0 ]; then
    echo "ALL CROSS-VALIDATION PASSED ✓ ($cross_ok graphs)"
else
    echo "CROSS-VALIDATION: $cross_ok passed, $cross_errors FAILED"
fi

# ── Phase 5: Generate markdown report ────────────────────────────────────
echo ""
echo "============================================"
echo "  Phase 4: Generating report"
echo "============================================"

REPORT="$RESULTS/report.md"

cat > "$REPORT" << 'HEADER'
# Combinatorial Suite — Full Test Report

HEADER

echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')" >> "$REPORT"
echo "" >> "$REPORT"

# ── Summary
val_pass=$(grep -c ",PASS" "$CSV" || true)
val_fail=$(grep -c ",FAIL" "$CSV" || true)
val_none=$(grep -c ",NONE" "$CSV" || true)
val_err=$(grep -c ",ERR," "$CSV" || true)

cat >> "$REPORT" << EOF
## Summary

| Metric | Count |
|--------|-------|
| Total runs | $total |
| Validation PASS | $val_pass |
| Validation FAIL | $val_fail |
| No validation | $val_none |
| Errors / timeouts | $val_err |
| Cross-validation OK | $cross_ok |
| Cross-validation FAIL | $cross_errors |

EOF

# ── Per-graph tables: general
echo "## General Matching Results" >> "$REPORT"
echo "" >> "$REPORT"

for graph in $GENERAL_GRAPHS; do
    gname="$(basename "$graph" .txt)"
    short="$(echo "$gname" | sed 's/general_unweighted_//')"

    echo "### $short" >> "$REPORT"
    echo "" >> "$REPORT"
    echo "| Algorithm | C++ size | C++ ms | Rust size | Rust ms | Python size | Python ms |" >> "$REPORT"
    echo "|-----------|----------|--------|-----------|---------|-------------|-----------|" >> "$REPORT"

    for alg in $GENERAL_ALGOS $MV_PURE; do
        aname="$(echo "$alg" | sed 's/micali-vazirani-pure/mv-pure/' | sed 's/micali-vazirani/mv-hybrid/' | sed 's/edmonds-blossom-/eb-/' | sed 's/gabow-/g-/')"
        row="| $aname"
        for lang in cpp rust python; do
            line="$(grep "^$alg,$gname,$lang," "$CSV")"
            if [ -n "$line" ]; then
                sz="$(echo "$line" | cut -d, -f4)"
                tm="$(echo "$line" | cut -d, -f5)"
                row="$row | $sz | $tm"
            else
                row="$row | — | —"
            fi
        done
        echo "$row |" >> "$REPORT"
    done
    echo "" >> "$REPORT"
done

# ── Per-graph tables: bipartite
echo "## Bipartite Matching Results" >> "$REPORT"
echo "" >> "$REPORT"

for graph in $BIPARTITE_GRAPHS; do
    gname="$(basename "$graph" .txt)"
    short="$(echo "$gname" | sed 's/bipartite_unweighted_//')"

    echo "### $short" >> "$REPORT"
    echo "" >> "$REPORT"
    echo "| Algorithm | C++ size | C++ ms | Rust size | Rust ms | Python size | Python ms |" >> "$REPORT"
    echo "|-----------|----------|--------|-----------|---------|-------------|-----------|" >> "$REPORT"

    for alg in $BIPARTITE_ALGOS; do
        row="| $alg"
        for lang in cpp rust python; do
            line="$(grep "^$alg,$gname,$lang," "$CSV")"
            if [ -n "$line" ]; then
                sz="$(echo "$line" | cut -d, -f4)"
                tm="$(echo "$line" | cut -d, -f5)"
                row="$row | $sz | $tm"
            else
                row="$row | — | —"
            fi
        done
        echo "$row |" >> "$REPORT"
    done
    echo "" >> "$REPORT"
done

# ── Validation failures detail
if [ "$cross_errors" -gt 0 ] || [ "$val_fail" -gt 0 ]; then
    echo "## Failures" >> "$REPORT"
    echo "" >> "$REPORT"
    echo '```' >> "$REPORT"
    grep ",FAIL" "$CSV" >> "$REPORT" || true
    echo '```' >> "$REPORT"
    echo "" >> "$REPORT"
fi

echo "" >> "$REPORT"
echo "---" >> "$REPORT"
echo "*Sequential execution on single core. Times are wall-clock milliseconds reported by each implementation.*" >> "$REPORT"

echo ""
echo "  Report written to: results/report.md"
echo "  Raw CSV data:      results/raw/all_results.csv"
echo "  Individual logs:   results/raw/"
echo ""

# ── final verdict ─────────────────────────────────────────────────────────
if [ "$cross_errors" -eq 0 ] && [ "$val_fail" -eq 0 ] && [ "$compile_errors" -eq 0 ]; then
    echo "============================================"
    echo "  ALL TESTS PASSED ✓"
    echo "============================================"
    exit 0
else
    echo "============================================"
    echo "  SOME ISSUES DETECTED — see report"
    echo "============================================"
    exit 1
fi
