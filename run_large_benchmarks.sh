#!/bin/sh
# run_large_benchmarks.sh – Run selected algorithms on large sparse graphs.
#
# Usage:
#   ./run_large_benchmarks.sh                                    # defaults
#   ./run_large_benchmarks.sh --sizes 100k 1m --langs cpp rust
#   ./run_large_benchmarks.sh --algos hk mv-pure gabow-opt
#   ./run_large_benchmarks.sh --mode plain greedy greedy-md
#   ./run_large_benchmarks.sh --runs 5 --timeout 600
#   ./run_large_benchmarks.sh --list
#
# Defaults:
#   sizes:    all found in data/large-benchmarks/
#   langs:    cpp rust python
#   algos:    all (filtered by feasibility)
#   mode:     plain (options: plain, greedy, greedy-md)
#   runs:     3 (reports median)
#   timeout:  300s per run
#   datadir:  data/large-benchmarks
#
# Produces:
#   results/large-benchmarks/<timestamp>/report.md
#   results/large-benchmarks/<timestamp>/results.csv
#   results/large-benchmarks/<timestamp>/raw/          (individual run logs)

set -e

# ── defaults ──────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$SCRIPT_DIR"
ALGO="$REPO/algorithms"
DATADIR="$REPO/data/large-benchmarks"
OUTDIR="$REPO/results/large-benchmarks"
RUNS=3
TIMEOUT=300
LIST_ONLY=0

# Filters (empty = all)
F_SIZES=""
F_LANGS=""
F_ALGOS=""
F_MODE=""

# ── algorithm registry ────────────────────────────────────────────────
# Short name → directory name, graph type, complexity class
# Graph type: general | bipartite
# Complexity: ve (O(VE)) | fast (O(E√V))

ALL_GENERAL="edmonds-simple edmonds-opt gabow-simple gabow-opt mv-pure"
ALL_BIPARTITE="hk"
ALL_ALGOS="$ALL_GENERAL $ALL_BIPARTITE"

alg_dir() {
    case "$1" in
        edmonds-simple) echo "edmonds-blossom-simple" ;;
        edmonds-opt)    echo "edmonds-blossom-optimized" ;;
        gabow-simple)   echo "gabow-simple" ;;
        gabow-opt)      echo "gabow-optimized" ;;
        mv-pure)        echo "micali-vazirani-pure" ;;
        hk)             echo "hopcroft-karp" ;;
    esac
}

alg_src() {
    # directory name → source file base (tr - _)
    alg_dir "$1" | tr '-' '_'
}

alg_complexity() {
    case "$1" in
        edmonds-simple|edmonds-opt|gabow-simple) echo "ve" ;;
        gabow-opt|mv-pure|hk) echo "fast" ;;
    esac
}

alg_type() {
    case "$1" in
        hk) echo "bipartite" ;;
        *)  echo "general" ;;
    esac
}

# ── parse arguments ───────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        --sizes)
            shift
            while [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; do
                F_SIZES="$F_SIZES $1"; shift
            done
            ;;
        --langs)
            shift
            while [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; do
                F_LANGS="$F_LANGS $1"; shift
            done
            ;;
        --algos)
            shift
            while [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; do
                F_ALGOS="$F_ALGOS $1"; shift
            done
            ;;
        --mode)
            shift
            while [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; do
                F_MODE="$F_MODE $1"; shift
            done
            ;;
        --runs)    shift; RUNS="$1"; shift ;;
        --timeout) shift; TIMEOUT="$1"; shift ;;
        --datadir) shift; DATADIR="$1"; shift ;;
        --outdir)  shift; OUTDIR="$1"; shift ;;
        --list)    LIST_ONLY=1; shift ;;
        --help|-h)
            sed -n '2,/^$/p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Apply defaults
[ -z "$F_LANGS" ] && F_LANGS="cpp rust python"

# Timestamp the output directory
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
OUTDIR="${OUTDIR}/${TIMESTAMP}"

# ── timeout command (macOS compat) ────────────────────────────────────
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

# ── discover available graph files ────────────────────────────────────
if [ ! -d "$DATADIR" ]; then
    echo "ERROR: Data directory not found: $DATADIR"
    echo "Run generate_large_benchmarks.py first."
    exit 1
fi

GENERAL_FILES=""
BIPARTITE_FILES=""

for f in "$DATADIR"/general_sparse_*.txt; do
    [ -f "$f" ] || continue
    # Extract size tag: general_sparse_100k_3.txt → 100k
    tag="$(basename "$f" .txt | sed 's/general_sparse_//' | sed 's/_[0-9]*$//')"
    if [ -n "$F_SIZES" ]; then
        echo "$F_SIZES" | grep -qw "$tag" || continue
    fi
    GENERAL_FILES="$GENERAL_FILES $f"
done

for f in "$DATADIR"/bipartite_sparse_*.txt; do
    [ -f "$f" ] || continue
    tag="$(basename "$f" .txt | sed 's/bipartite_sparse_//' | sed 's/_[0-9]*$//')"
    if [ -n "$F_SIZES" ]; then
        echo "$F_SIZES" | grep -qw "$tag" || continue
    fi
    BIPARTITE_FILES="$BIPARTITE_FILES $f"
done

# ── determine which algorithms to run ─────────────────────────────────
# Extract max vertex count from filenames to decide feasibility
max_v=0
for f in $GENERAL_FILES $BIPARTITE_FILES; do
    # Parse header line for vertex count
    v="$(head -1 "$f" | awk '{print $1}')"
    [ "$v" -gt "$max_v" ] 2>/dev/null && max_v="$v"
done

select_algos() {
    gtype="$1"  # general or bipartite
    result=""
    for alg in $ALL_ALGOS; do
        # Filter by user request
        if [ -n "$F_ALGOS" ]; then
            echo "$F_ALGOS" | grep -qw "$alg" || continue
        fi
        # Filter by graph type
        [ "$(alg_type "$alg")" = "$gtype" ] || continue
        result="$result $alg"
    done
    echo "$result"
}

GENERAL_ALGOS="$(select_algos general)"
BIPARTITE_ALGOS="$(select_algos bipartite)"

# Apply defaults
[ -z "$F_MODE" ] && F_MODE="plain"

# Validate mode values
for m in $F_MODE; do
    case "$m" in
        plain|greedy|greedy-md) ;;
        *) echo "ERROR: --mode values must be plain, greedy, or greedy-md (got: $m)"; exit 1 ;;
    esac
done

# ── build plan ────────────────────────────────────────────────────────
PLAN=""
plan_count=0

add_to_plan() {
    alg="$1"; graph="$2"; lang="$3"; mode="$4"
    gname="$(basename "$graph" .txt)"

    # Parse V from header
    v="$(head -1 "$graph" | awk '{print $1}')"

    # Skip O(VE) algorithms on large graphs (> 200K vertices) unless explicitly requested
    if [ -z "$F_ALGOS" ] && [ "$(alg_complexity "$alg")" = "ve" ] && [ "$v" -gt 200000 ]; then
        return
    fi

    plan_count=$((plan_count + 1))
    PLAN="$PLAN
$alg|$graph|$lang|$gname|$v|$mode"
}

add_modes() {
    alg="$1"; graph="$2"; lang="$3"
    for m in $F_MODE; do
        add_to_plan "$alg" "$graph" "$lang" "$m"
    done
}

for alg in $GENERAL_ALGOS; do
    for graph in $GENERAL_FILES; do
        for lang in $F_LANGS; do
            add_modes "$alg" "$graph" "$lang"
        done
    done
done

for alg in $BIPARTITE_ALGOS; do
    for graph in $BIPARTITE_FILES; do
        for lang in $F_LANGS; do
            add_modes "$alg" "$graph" "$lang"
        done
    done
done

# ── display plan ──────────────────────────────────────────────────────
echo "============================================="
echo "  Large-Scale Benchmark Plan"
echo "============================================="
echo ""
echo "  Run:            $TIMESTAMP"
echo "  Output:         $OUTDIR"
echo "  Data directory: $DATADIR"
echo "  Languages:     $F_LANGS"
echo "  Mode:           $F_MODE"
echo "  Runs per job:   $RUNS (report median)"
echo "  Timeout:        ${TIMEOUT}s per run"
echo "  Total jobs:     $plan_count"
echo ""

printf "  %-18s %-35s %-8s %-10s %s\n" "Algorithm" "Graph" "Lang" "Mode" "V"
printf "  %-18s %-35s %-8s %-10s %s\n" "---------" "-----" "----" "------" "-"
echo "$PLAN" | while IFS='|' read -r alg graph lang gname v greedy; do
    [ -z "$alg" ] && continue
    printf "  %-18s %-35s %-8s %-10s %s\n" "$alg" "$gname" "$lang" "$greedy" "$v"
done
echo ""

if [ "$LIST_ONLY" -eq 1 ]; then
    echo "(--list mode, not running.)"
    exit 0
fi

# ── compile ───────────────────────────────────────────────────────────
echo "============================================="
echo "  Compiling"
echo "============================================="
echo ""

compiled=""
compile_errors=0

needs_compile() {
    alg="$1"; lang="$2"
    echo "$compiled" | grep -q "${alg}:${lang}" && return 1
    return 0
}

mark_compiled() {
    compiled="$compiled ${1}:${2}"
}

echo "$PLAN" | while IFS='|' read -r alg graph lang gname v greedy; do
    [ -z "$alg" ] && continue
    [ "$lang" = "python" ] && continue
    needs_compile "$alg" "$lang" || continue

    dir="$(alg_dir "$alg")"
    base="$(alg_src "$alg")"

    case "$lang" in
        cpp)
            src="$ALGO/$dir/cpp/${base}.cpp"
            bin="$ALGO/$dir/cpp/${base}_cpp"
            [ -f "$src" ] || continue
            printf "  compile %-20s %-6s " "$alg" "cpp"
            if g++ -O3 -std=c++17 "$src" -o "$bin" 2>/dev/null; then
                echo "✓"
            else
                echo "✗"
                compile_errors=$((compile_errors + 1))
            fi
            ;;
        rust)
            src="$ALGO/$dir/rust/${base}.rs"
            bin="$ALGO/$dir/rust/${base}_rust"
            [ -f "$src" ] || continue
            printf "  compile %-20s %-6s " "$alg" "rust"
            if rustc -O "$src" -o "$bin" 2>/dev/null; then
                echo "✓"
            else
                echo "✗"
                compile_errors=$((compile_errors + 1))
            fi
            ;;
    esac
    mark_compiled "$alg" "$lang"
done

echo ""

# ── run benchmarks ────────────────────────────────────────────────────
echo "============================================="
echo "  Running ($RUNS runs each, reporting median)"
echo "============================================="
echo ""

mkdir -p "$OUTDIR/raw"
CSV="$OUTDIR/results.csv"
echo "algo,graph,lang,vertices,mode,matching_size,greedy_init_size,greedy_pct,median_ms,run1_ms,run2_ms,run3_ms,validation" > "$CSV"

job=0
echo "$PLAN" | while IFS='|' read -r alg graph lang gname v greedy; do
    [ -z "$alg" ] && continue
    job=$((job + 1))

    dir="$(alg_dir "$alg")"
    base="$(alg_src "$alg")"
    logbase="$OUTDIR/raw/${base}_${lang}_${gname}_${greedy}"

    printf "  [%3d/%d] %-18s %-6s %-10s %-30s " "$job" "$plan_count" "$alg" "$lang" "$greedy" "$gname"

    # Determine executable
    case "$lang" in
        cpp)
            bin="$ALGO/$dir/cpp/${base}_cpp"
            [ -x "$bin" ] || { echo "SKIP (not compiled)"; continue; }
            cmd="$bin"
            ;;
        rust)
            bin="$ALGO/$dir/rust/${base}_rust"
            [ -x "$bin" ] || { echo "SKIP (not compiled)"; continue; }
            cmd="$bin"
            ;;
        python)
            src="$ALGO/$dir/python/${base}.py"
            [ -f "$src" ] || { echo "SKIP (not found)"; continue; }
            cmd="uv run $src"
            ;;
    esac

    # Build command args
    extra_args=""
    [ "$greedy" = "greedy" ] && extra_args="--greedy"
    [ "$greedy" = "greedy-md" ] && extra_args="--greedy-md"

    # Run N times
    times=""
    size="ERR"
    greedy_init="NA"
    greedy_pct="NA"
    valid="NONE"
    status="OK"

    run_i=0
    while [ "$run_i" -lt "$RUNS" ]; do
        run_i=$((run_i + 1))
        logfile="${logbase}_run${run_i}.log"

        if run_with_timeout "$TIMEOUT" $cmd "$graph" $extra_args > "$logfile" 2>&1; then
            t="$(grep '^Time:' "$logfile" | awk '{print $2}')"
            s="$(grep '^Matching size:' "$logfile" | tail -1 | awk '{print $3}')"
            gi="$(grep '^Greedy init size:' "$logfile" | awk '{print $4}')"
            gp="$(grep '^Greedy/Final:' "$logfile" | awk '{print $2}')"
            vl="$(grep 'VALIDATION' "$logfile" | head -1)"

            [ -n "$t" ] && times="$times $t" || times="$times ERR"
            [ -n "$s" ] && size="$s"
            [ -n "$gi" ] && greedy_init="$gi"
            [ -n "$gp" ] && greedy_pct="$gp"

            case "$vl" in
                *PASSED*) valid="PASS" ;;
                *FAILED*) valid="FAIL"; status="FAIL" ;;
            esac
        else
            times="$times TIMEOUT"
            status="TIMEOUT"
        fi
    done

    # Compute median
    clean_times="$(echo "$times" | tr ' ' '\n' | grep -v ERR | grep -v TIMEOUT | sort -n)"
    n_good="$(echo "$clean_times" | grep -c . || true)"

    if [ "$n_good" -gt 0 ]; then
        med="$(echo "$clean_times" | awk -v n="$n_good" 'NR==int((n+1)/2){print;exit}')"
    else
        med="ERR"
    fi

    # Pad times to 3 fields for CSV
    t1="$(echo "$times" | awk '{print $1}')"
    t2="$(echo "$times" | awk '{print $2}')"
    t3="$(echo "$times" | awk '{print $3}')"
    [ -z "$t1" ] && t1="-"
    [ -z "$t2" ] && t2="-"
    [ -z "$t3" ] && t3="-"

    echo "$alg,$gname,$lang,$v,$greedy,$size,$greedy_init,$greedy_pct,$med,$t1,$t2,$t3,$valid" >> "$CSV"

    if [ "$greedy" = "greedy" ] || [ "$greedy" = "greedy-md" ]; then
        printf "size=%-8s median=%-8s %-6s greedy_init=%-8s (%s)\n" "$size" "${med}ms" "$valid" "$greedy_init" "$greedy_pct"
    else
        printf "size=%-8s median=%-8s %s\n" "$size" "${med}ms" "$valid"
    fi
done

# ── cross-validation ──────────────────────────────────────────────────
echo ""
echo "============================================="
echo "  Cross-Validation"
echo "============================================="
echo ""

cross_ok=0
cross_fail=0

# Collect unique graph names from CSV
graphs="$(tail -n +2 "$CSV" | cut -d, -f2 | sort -u)"

for gname in $graphs; do
    # Compare matching sizes across all algos/langs/modes for this graph
    sizes="$(grep ",$gname," "$CSV" | cut -d, -f6 | grep -v ERR | sort -u)"
    count="$(echo "$sizes" | grep -c . || true)"
    if [ "$count" -eq 1 ] && [ -n "$sizes" ]; then
        printf "  %-45s size=%-8s ✓\n" "$gname" "$sizes"
        cross_ok=$((cross_ok + 1))
    elif [ "$count" -eq 0 ]; then
        printf "  %-45s NO DATA\n" "$gname"
    else
        printf "  %-45s MISMATCH: %s ✗\n" "$gname" "$(echo $sizes | tr '\n' ' ')"
        cross_fail=$((cross_fail + 1))
    fi
done

echo ""

# ── generate report ───────────────────────────────────────────────────
echo "============================================="
echo "  Generating Report"
echo "============================================="

REPORT="$OUTDIR/report.md"

cat > "$REPORT" << 'EOF'
# Large-Scale Benchmark Report

EOF

echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')" >> "$REPORT"
echo "Runs per job: $RUNS (median reported)" >> "$REPORT"
echo "Timeout: ${TIMEOUT}s" >> "$REPORT"
echo "Mode: $F_MODE" >> "$REPORT"
echo "" >> "$REPORT"

# Summary
total_runs="$(tail -n +2 "$CSV" | wc -l | tr -d ' ')"
pass_count="$(grep -c ',PASS' "$CSV" || true)"
fail_count="$(grep -c ',FAIL' "$CSV" || true)"

cat >> "$REPORT" << EOF
## Summary

| Metric | Count |
|--------|------:|
| Total jobs | $total_runs |
| Validation PASS | $pass_count |
| Validation FAIL | $fail_count |
| Cross-validation OK | $cross_ok |
| Cross-validation FAIL | $cross_fail |

EOF

# Results table per graph – show plain vs greedy side by side
echo "## Results" >> "$REPORT"
echo "" >> "$REPORT"

for gname in $graphs; do
    echo "### $gname" >> "$REPORT"
    echo "" >> "$REPORT"

    # Determine which languages ran on this graph
    graph_langs="$(grep ",$gname," "$CSV" | cut -d, -f3 | sort -u)"

    # Build header: for each lang, show plain, greedy, greedy-md
    header="| Algorithm"
    separator="|----------"
    for lang in $graph_langs; do
        header="$header | ${lang} plain | ${lang} greedy | ${lang} grdy-md | ${lang} md-init | ${lang} md-%"
        separator="$separator|-------:|-------:|-------:|-------:|-------:"
    done
    header="$header | size |"
    separator="$separator|-------:|"
    echo "$header" >> "$REPORT"
    echo "$separator" >> "$REPORT"

    # Determine which algorithms ran on this graph
    graph_algos="$(grep ",$gname," "$CSV" | cut -d, -f1 | sort -u)"

    for alg in $graph_algos; do
        row="| $alg"
        alg_size=""
        for lang in $graph_langs; do
            plain_line="$(grep "^$alg,$gname,$lang,.*,plain," "$CSV")"
            greedy_line="$(grep "^$alg,$gname,$lang,.*,greedy," "$CSV" | grep -v greedy-md)"
            greedymd_line="$(grep "^$alg,$gname,$lang,.*,greedy-md," "$CSV")"

            if [ -n "$plain_line" ]; then
                plain_ms="$(echo "$plain_line" | cut -d, -f9)"
                alg_size="$(echo "$plain_line" | cut -d, -f6)"
            else
                plain_ms="–"
            fi

            if [ -n "$greedy_line" ]; then
                greedy_ms="$(echo "$greedy_line" | cut -d, -f9)"
            else
                greedy_ms="–"
            fi

            if [ -n "$greedymd_line" ]; then
                greedymd_ms="$(echo "$greedymd_line" | cut -d, -f9)"
                greedymd_init="$(echo "$greedymd_line" | cut -d, -f7)"
                greedymd_pct="$(echo "$greedymd_line" | cut -d, -f8)"
            else
                greedymd_ms="–"
                greedymd_init="–"
                greedymd_pct="–"
            fi

            row="$row | $plain_ms | $greedy_ms | $greedymd_ms | $greedymd_init | $greedymd_pct"
        done
        echo "$row | $alg_size |" >> "$REPORT"
    done
    echo "" >> "$REPORT"
done

# Greedy speedup table
echo "## Greedy Bootstrap Speedup" >> "$REPORT"
echo "" >> "$REPORT"
echo "| Algorithm | Graph | Lang | Plain ms | Greedy ms | Speedup | Grdy-MD ms | MD Speedup | MD Init | MD % |" >> "$REPORT"
echo "|-----------|-------|------|--------:|---------:|--------:|---------:|--------:|--------:|--------:|" >> "$REPORT"

for gname in $graphs; do
    graph_algos="$(grep ",$gname," "$CSV" | cut -d, -f1 | sort -u)"
    graph_langs="$(grep ",$gname," "$CSV" | cut -d, -f3 | sort -u)"
    for alg in $graph_algos; do
        for lang in $graph_langs; do
            plain_ms="$(grep "^$alg,$gname,$lang,.*,plain," "$CSV" | cut -d, -f9)"
            greedy_ms="$(grep "^$alg,$gname,$lang,.*,greedy," "$CSV" | grep -v greedy-md | cut -d, -f9)"
            greedymd_ms="$(grep "^$alg,$gname,$lang,.*,greedy-md," "$CSV" | cut -d, -f9)"
            greedymd_init="$(grep "^$alg,$gname,$lang,.*,greedy-md," "$CSV" | cut -d, -f7)"
            greedymd_pct="$(grep "^$alg,$gname,$lang,.*,greedy-md," "$CSV" | cut -d, -f8)"

            [ -z "$plain_ms" ] && plain_ms="–"
            [ -z "$greedy_ms" ] && greedy_ms="–"
            [ -z "$greedymd_ms" ] && greedymd_ms="–"
            [ -z "$greedymd_init" ] && greedymd_init="–"
            [ -z "$greedymd_pct" ] && greedymd_pct="–"

            # Compute greedy speedup
            if [ "$plain_ms" != "–" ] && [ "$plain_ms" != "ERR" ] && \
               [ "$greedy_ms" != "–" ] && [ "$greedy_ms" != "ERR" ] && \
               [ "$greedy_ms" -gt 0 ] 2>/dev/null; then
                speedup="$(awk "BEGIN{printf \"%.2f\", $plain_ms / $greedy_ms}")"
            else
                speedup="–"
            fi

            # Compute greedy-md speedup
            if [ "$plain_ms" != "–" ] && [ "$plain_ms" != "ERR" ] && \
               [ "$greedymd_ms" != "–" ] && [ "$greedymd_ms" != "ERR" ] && \
               [ "$greedymd_ms" -gt 0 ] 2>/dev/null; then
                md_speedup="$(awk "BEGIN{printf \"%.2f\", $plain_ms / $greedymd_ms}")"
            else
                md_speedup="–"
            fi

            short_gname="$(echo "$gname" | sed 's/general_sparse_/g_/' | sed 's/bipartite_sparse_/b_/')"
            echo "| $alg | $short_gname | $lang | $plain_ms | $greedy_ms | ${speedup}× | $greedymd_ms | ${md_speedup}× | $greedymd_init | $greedymd_pct |" >> "$REPORT"
        done
    done
done

echo "" >> "$REPORT"

# Language speedup table (C++ = 1.0×) – using plain mode for fair comparison
echo "## Language Speedups (C++ = 1.0×, plain mode)" >> "$REPORT"
echo "" >> "$REPORT"
echo "| Algorithm | Graph | C++ ms | Rust ms | Rust/C++ | Python ms | Python/C++ |" >> "$REPORT"
echo "|-----------|-------|-------:|--------:|---------:|----------:|-----------:|" >> "$REPORT"

for gname in $graphs; do
    graph_algos="$(grep ",$gname," "$CSV" | cut -d, -f1 | sort -u)"
    for alg in $graph_algos; do
        cpp_ms="$(grep "^$alg,$gname,cpp,.*,plain," "$CSV" 2>/dev/null | cut -d, -f9)"
        rust_ms="$(grep "^$alg,$gname,rust,.*,plain," "$CSV" 2>/dev/null | cut -d, -f9)"
        py_ms="$(grep "^$alg,$gname,python,.*,plain," "$CSV" 2>/dev/null | cut -d, -f9)"

        [ -z "$cpp_ms" ] && cpp_ms="–"
        [ -z "$rust_ms" ] && rust_ms="–"
        [ -z "$py_ms" ] && py_ms="–"

        # Compute ratios
        if [ "$cpp_ms" != "–" ] && [ "$cpp_ms" != "ERR" ] && [ "$cpp_ms" -gt 0 ] 2>/dev/null; then
            if [ "$rust_ms" != "–" ] && [ "$rust_ms" != "ERR" ] 2>/dev/null; then
                rust_ratio="$(awk "BEGIN{printf \"%.1f\", $rust_ms / $cpp_ms}")"
            else
                rust_ratio="–"
            fi
            if [ "$py_ms" != "–" ] && [ "$py_ms" != "ERR" ] 2>/dev/null; then
                py_ratio="$(awk "BEGIN{printf \"%.1f\", $py_ms / $cpp_ms}")"
            else
                py_ratio="–"
            fi
        else
            rust_ratio="–"
            py_ratio="–"
        fi

        short_gname="$(echo "$gname" | sed 's/general_sparse_/g_/' | sed 's/bipartite_sparse_/b_/')"
        # Format ratios: append × only if not a dash
        if [ "$rust_ratio" = "–" ]; then rust_fmt="–"; else rust_fmt="${rust_ratio}×"; fi
        if [ "$py_ratio" = "–" ]; then py_fmt="–"; else py_fmt="${py_ratio}×"; fi
        echo "| $alg | $short_gname | $cpp_ms | $rust_ms | $rust_fmt | $py_ms | $py_fmt |" >> "$REPORT"
    done
done

echo "" >> "$REPORT"
echo "---" >> "$REPORT"
echo "*Median of $RUNS runs. Wall-clock ms reported by each implementation. Timeout: ${TIMEOUT}s.*" >> "$REPORT"

echo ""
echo "  Run:     $TIMESTAMP"
echo "  Report:  $OUTDIR/report.md"
echo "  CSV:     $OUTDIR/results.csv"
echo "  Logs:    $OUTDIR/raw/"

# Create/update 'latest' symlink
BASEDIR="$(dirname "$OUTDIR")"
ln -sfn "$TIMESTAMP" "$BASEDIR/latest"
echo "  Latest:  $BASEDIR/latest -> $TIMESTAMP"
echo ""

# ── verdict ───────────────────────────────────────────────────────────
if [ "$cross_fail" -eq 0 ] && [ "$fail_count" -eq 0 ]; then
    echo "============================================="
    echo "  ALL VALIDATIONS PASSED ✓"
    echo "============================================="
    exit 0
else
    echo "============================================="
    echo "  ISSUES DETECTED – see report"
    echo "============================================="
    exit 1
fi
