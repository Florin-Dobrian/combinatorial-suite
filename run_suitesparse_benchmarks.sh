#!/bin/sh
# run_suitesparse_benchmarks.sh – Run matching algorithms on SuiteSparse graphs.
#
# Usage:
#   ./run_suitesparse_benchmarks.sh                                    # defaults
#   ./run_suitesparse_benchmarks.sh --langs cpp rust
#   ./run_suitesparse_benchmarks.sh --algos g2-csr mv-csr
#   ./run_suitesparse_benchmarks.sh --mode plain greedy greedy-md
#   ./run_suitesparse_benchmarks.sh --graphs fe_body auto ecology1
#   ./run_suitesparse_benchmarks.sh --runs 5 --timeout 600
#   ./run_suitesparse_benchmarks.sh --list
#
# Defaults:
#   graphs:   all .txt files in data/general-unweighted/suitesparse/
#   langs:    cpp rust python
#   algos:    auto (all 5 for V ≤ 200k, GO+MV for larger)
#   mode:     plain greedy-md
#   runs:     3 (reports median)
#   timeout:  600s per run
#
# Produces:
#   results/suitesparse/<timestamp>/report.md
#   results/suitesparse/<timestamp>/results.csv
#   results/suitesparse/<timestamp>/logs/

set -e

# ── defaults ──────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$SCRIPT_DIR"
ALGO="$REPO/algorithms"
DATADIR="$REPO/data/general-unweighted/suitesparse"
OUTDIR="$REPO/results/suitesparse"
RUNS=3
TIMEOUT=600
LIST_ONLY=0

# Compiler flags (single source of truth; used at compile time and in environment.txt)
CPP_FLAGS="-O3 -std=c++17"
RUST_FLAGS="-O"

# Capture original command-line args for environment.txt (before parsing consumes them)
ORIG_ARGS="$*"

# Filters (empty = all)
F_GRAPHS=""
F_LANGS=""
F_ALGOS=""
F_MODE=""

# ── algorithm registry ────────────────────────────────────────────────
ALL_ALGOS="e1-vv e1-csr e2-vv e2-csr g1-vv g1-csr g2-vv g2-csr mv-vv mv-csr"

alg_dir() {
    case "$1" in
        e1-vv)  echo "e1-vv" ;;
        e1-csr) echo "e1-csr" ;;
        e2-vv)  echo "e2-vv" ;;
        e2-csr) echo "e2-csr" ;;
        g1-vv)  echo "g1-vv" ;;
        g1-csr) echo "g1-csr" ;;
        g2-vv)  echo "g2-vv" ;;
        g2-csr) echo "g2-csr" ;;
        mv-vv)  echo "mv-vv" ;;
        mv-csr) echo "mv-csr" ;;
    esac
}

alg_src() {
    alg_dir "$1" | tr '-' '_'
}

alg_complexity() {
    case "$1" in
        e1-vv|e1-csr|e2-vv|e2-csr|g1-vv|g1-csr) echo "ve" ;;
        g2-vv|g2-csr|mv-vv|mv-csr) echo "fast" ;;
    esac
}

# ── parse arguments ───────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        --graphs)
            shift
            while [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; do
                F_GRAPHS="$F_GRAPHS $1"; shift
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
[ -z "$F_MODE" ] && F_MODE="plain greedy-md"

# Validate mode values
for m in $F_MODE; do
    case "$m" in
        plain|greedy|greedy-md) ;;
        *) echo "ERROR: --mode values must be plain, greedy, or greedy-md (got: $m)"; exit 1 ;;
    esac
done

# Hostname-prefixed run ID: groups runs by machine when sorted
HOST_SHORT="$(hostname -s)"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
RUN_ID="${HOST_SHORT}_${TIMESTAMP}"
OUTDIR="${OUTDIR}/${RUN_ID}"

# ── Environment capture ───────────────────────────────────────────────
ENV_FILE=""
RUN_START_EPOCH=""

write_env_header() {
    mkdir -p "$OUTDIR"
    ENV_FILE="$OUTDIR/environment.txt"
    RUN_START_EPOCH="$(date +%s)"
    {
        echo "Run: $RUN_ID"
        echo "Date: $(date '+%Y-%m-%d %H:%M:%S %z (%Z)')"
        echo "Host: $HOST_SHORT"
        echo "Hostname: $(hostname)"
        echo ""
        echo "== Command =="
        echo "${0##*/} $ORIG_ARGS"
        echo ""
        echo "== Run parameters =="
        echo "Algorithms: $(echo ${F_ALGOS:-(all)})"
        echo "Languages: $(echo ${F_LANGS:-(all)})"
        echo "Modes: $(echo ${F_MODE:-(all)})"
        echo "Graphs: $(echo ${F_GRAPHS:-(all)})"
        echo "Runs per job: $RUNS"
        echo "Timeout: ${TIMEOUT}s"
        echo ""
        echo "== System =="
        if [ "$(uname)" = "Darwin" ]; then
            echo "OS: macOS $(sw_vers -productVersion 2>/dev/null)"
            echo "Kernel: $(uname -sr)"
            echo "Model: $(sysctl -n hw.model 2>/dev/null)"
            echo "CPU: $(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
            P="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)"
            E="$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null)"
            TOTAL="$(sysctl -n hw.ncpu 2>/dev/null)"
            if [ -n "$P" ] && [ -n "$E" ]; then
                echo "Cores: $TOTAL (${P} performance + ${E} efficiency)"
            else
                echo "Cores: $TOTAL"
            fi
            MEM_BYTES="$(sysctl -n hw.memsize 2>/dev/null)"
            if [ -n "$MEM_BYTES" ]; then
                echo "RAM: $(( MEM_BYTES / 1024 / 1024 / 1024 )) GB"
            fi
        else
            echo "OS: $(uname -sr)"
            echo "CPU: $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | sed 's/.*: //')"
            echo "Cores: $(nproc 2>/dev/null)"
            MEM_KB="$(grep MemTotal /proc/meminfo 2>/dev/null | awk '{print $2}')"
            if [ -n "$MEM_KB" ]; then
                echo "RAM: $(( MEM_KB / 1024 / 1024 )) GB"
            fi
        fi
        echo ""
        echo "== Compilers =="
        if command -v g++ >/dev/null 2>&1; then
            echo "C++ (g++): $(g++ --version | head -1)"
            echo "  Path: $(command -v g++)"
            echo "  Build flags: $CPP_FLAGS"
        fi
        if command -v rustc >/dev/null 2>&1; then
            echo "Rust (rustc): $(rustc --version)"
            echo "  Path: $(command -v rustc)"
            echo "  Build flags: $RUST_FLAGS"
        fi
        echo ""
        echo "== Repository =="
        if git -C "$REPO" rev-parse --git-dir >/dev/null 2>&1; then
            REMOTE="$(git -C "$REPO" config --get remote.origin.url 2>/dev/null)"
            BRANCH="$(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null)"
            COMMIT="$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null)"
            DIRTY="$(git -C "$REPO" status --porcelain 2>/dev/null)"
            [ -n "$REMOTE" ] && echo "Remote: $REMOTE"
            [ -n "$BRANCH" ] && echo "Branch: $BRANCH"
            [ -n "$COMMIT" ] && echo "Commit: $COMMIT"
            if [ -z "$DIRTY" ]; then
                echo "Working tree: clean"
            else
                echo "Working tree: dirty (uncommitted changes present)"
            fi
        fi
        echo ""
    } > "$ENV_FILE"
}

write_env_footer() {
    RUN_END_EPOCH="$(date +%s)"
    ELAPSED=$(( RUN_END_EPOCH - RUN_START_EPOCH ))
    H=$(( ELAPSED / 3600 ))
    M=$(( (ELAPSED % 3600) / 60 ))
    S=$(( ELAPSED % 60 ))
    {
        echo "== Timing =="
        echo "Started: $(date -r "$RUN_START_EPOCH" '+%Y-%m-%d %H:%M:%S %z (%Z)')"
        echo "Finished: $(date -r "$RUN_END_EPOCH" '+%Y-%m-%d %H:%M:%S %z (%Z)')"
        printf "Total: %dh %02dm %02ds\n" "$H" "$M" "$S"
    } >> "$ENV_FILE"
}

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

# ── discover graph files ──────────────────────────────────────────────
if [ ! -d "$DATADIR" ]; then
    echo "ERROR: Data directory not found: $DATADIR"
    echo "Download SuiteSparse matrices and convert with mtx_to_edgelist.py."
    exit 1
fi

GRAPH_FILES=""
for f in "$DATADIR"/*.txt; do
    [ -f "$f" ] || continue
    gname="$(basename "$f" .txt)"
    if [ -n "$F_GRAPHS" ]; then
        echo "$F_GRAPHS" | grep -qw "$gname" || continue
    fi
    GRAPH_FILES="$GRAPH_FILES $f"
done

if [ -z "$GRAPH_FILES" ]; then
    echo "ERROR: No graph files found in $DATADIR"
    exit 1
fi

# ── determine which algorithms to run ─────────────────────────────────
if [ -n "$F_ALGOS" ]; then
    RUN_ALGOS="$F_ALGOS"
else
    RUN_ALGOS="$ALL_ALGOS"
fi

# ── build plan ────────────────────────────────────────────────────────
PLAN=""
plan_count=0

add_to_plan() {
    alg="$1"; graph="$2"; lang="$3"; mode="$4"
    gname="$(basename "$graph" .txt)"
    v="$(head -1 "$graph" | awk '{print $1}')"

    # Auto-skip O(VE) algorithms on large graphs (> 200K) unless explicitly requested
    if [ -z "$F_ALGOS" ] && [ "$(alg_complexity "$alg")" = "ve" ] && [ "$v" -gt 200000 ]; then
        return
    fi

    # Auto-skip Python on very large graphs (> 2M) unless explicitly requested
    if [ "$lang" = "python" ] && [ "$v" -gt 2000000 ] && [ -z "$F_ALGOS" ]; then
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

for alg in $RUN_ALGOS; do
    for graph in $GRAPH_FILES; do
        for lang in $F_LANGS; do
            add_modes "$alg" "$graph" "$lang"
        done
    done
done

# ── display plan ──────────────────────────────────────────────────────
echo "============================================="
echo "  SuiteSparse Benchmark Plan"
echo "============================================="
echo ""
echo "  Run:            $RUN_ID"
echo "  Output:         $OUTDIR"
echo "  Data directory: $DATADIR"
echo "  Languages:     $F_LANGS"
echo "  Mode:           $F_MODE"
echo "  Runs per job:   $RUNS (report median)"
echo "  Timeout:        ${TIMEOUT}s per run"
echo "  Total jobs:     $plan_count"
echo ""

printf "  %-18s %-25s %-8s %-10s %s\n" "Algorithm" "Graph" "Lang" "Mode" "V"
printf "  %-18s %-25s %-8s %-10s %s\n" "---------" "-----" "----" "------" "-"
echo "$PLAN" | while IFS='|' read -r alg graph lang gname v greedy; do
    [ -z "$alg" ] && continue
    printf "  %-18s %-25s %-8s %-10s %s\n" "$alg" "$gname" "$lang" "$greedy" "$v"
done
echo ""

if [ "$LIST_ONLY" -eq 1 ]; then
    echo "(--list mode, not running.)"
    exit 0
fi

# Write environment.txt now that we know we're actually running
write_env_header

# ── compile ───────────────────────────────────────────────────────────
echo "============================================="
echo "  Compiling"
echo "============================================="
echo ""

compiled=""

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
            if g++ $CPP_FLAGS "$src" -o "$bin" 2>/dev/null; then
                echo "✓"
            else
                echo "✗"
            fi
            ;;
        rust)
            src="$ALGO/$dir/rs/${base}.rs"
            bin="$ALGO/$dir/rs/${base}_rs"
            [ -f "$src" ] || continue
            printf "  compile %-20s %-6s " "$alg" "rust"
            if rustc $RUST_FLAGS "$src" -o "$bin" 2>/dev/null; then
                echo "✓"
            else
                echo "✗"
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

mkdir -p "$OUTDIR/logs"
CSV="$OUTDIR/results.csv"
echo "algo,graph,lang,vertices,mode,matching_size,greedy_init_size,greedy_pct,phases,median_ms,run1_ms,run2_ms,run3_ms,validation,host" > "$CSV"

job=0
echo "$PLAN" | while IFS='|' read -r alg graph lang gname v greedy; do
    [ -z "$alg" ] && continue
    job=$((job + 1))

    dir="$(alg_dir "$alg")"
    base="$(alg_src "$alg")"
    logbase="$OUTDIR/logs/${base}_${lang}_${gname}_${greedy}"

    printf "  [%3d/%d] %-18s %-6s %-10s %-25s " "$job" "$plan_count" "$alg" "$lang" "$greedy" "$gname"

    # Determine executable
    case "$lang" in
        cpp)
            bin="$ALGO/$dir/cpp/${base}_cpp"
            [ -x "$bin" ] || { echo "SKIP (not compiled)"; continue; }
            cmd="$bin"
            ;;
        rust)
            bin="$ALGO/$dir/rs/${base}_rs"
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
    phases="NA"
    valid="NONE"

    run_i=0
    while [ "$run_i" -lt "$RUNS" ]; do
        run_i=$((run_i + 1))
        logfile="${logbase}_run${run_i}.log"

        if run_with_timeout "$TIMEOUT" $cmd "$graph" $extra_args > "$logfile" 2>&1; then
            t="$(grep '^Time:' "$logfile" | awk '{print $2}')"
            s="$(grep '^Matching size:' "$logfile" | tail -1 | awk '{print $3}')"
            gi="$(grep '^Greedy init size:' "$logfile" | awk '{print $4}')"
            gp="$(grep '^Greedy/Final:' "$logfile" | awk '{print $2}')"
            ph="$(grep '^Phases:' "$logfile" | awk '{print $2}')"
            vl="$(grep 'VALIDATION' "$logfile" | head -1)"

            [ -n "$t" ] && times="$times $t" || times="$times ERR"
            [ -n "$s" ] && size="$s"
            [ -n "$gi" ] && greedy_init="$gi"
            [ -n "$gp" ] && greedy_pct="$gp"
            [ -n "$ph" ] && phases="$ph"

            case "$vl" in
                *PASSED*) valid="PASS" ;;
                *FAILED*) valid="FAIL" ;;
            esac
        else
            times="$times TIMEOUT"
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

    echo "$alg,$gname,$lang,$v,$greedy,$size,$greedy_init,$greedy_pct,$phases,$med,$t1,$t2,$t3,$valid,$HOST_SHORT" >> "$CSV"

    if [ "$greedy" = "greedy" ] || [ "$greedy" = "greedy-md" ]; then
        printf "size=%-8s median=%-8s %-6s greedy_init=%-8s (%s)  phases=%s\n" "$size" "${med}ms" "$valid" "$greedy_init" "$greedy_pct" "$phases"
    else
        printf "size=%-8s median=%-8s %-6s phases=%s\n" "$size" "${med}ms" "$valid" "$phases"
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

graphs="$(tail -n +2 "$CSV" | cut -d, -f2 | sort -u)"

for gname in $graphs; do
    sizes="$(grep ",$gname," "$CSV" | cut -d, -f6 | grep -v ERR | sort -u)"
    count="$(echo "$sizes" | grep -c . || true)"
    if [ "$count" -eq 1 ] && [ -n "$sizes" ]; then
        printf "  %-35s size=%-8s ✓\n" "$gname" "$sizes"
        cross_ok=$((cross_ok + 1))
    elif [ "$count" -eq 0 ]; then
        printf "  %-35s NO DATA\n" "$gname"
    else
        printf "  %-35s MISMATCH: %s ✗\n" "$gname" "$(echo $sizes | tr '\n' ' ')"
        cross_fail=$((cross_fail + 1))
    fi
done

echo ""

# ── generate report ───────────────────────────────────────────────────
echo "============================================="
echo "  Generating Report"
echo "============================================="

# Disable set -e for this block: many grep calls legitimately return
# no matches (e.g., when only one mode was run), and that's not an error.
set +e

REPORT="$OUTDIR/report.md"

cat > "$REPORT" << 'EOF'
# SuiteSparse Benchmark Report

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

## Graph Properties

| Graph | V | E | Avg Degree | Source |
|-------|--:|--:|----------:|--------|
EOF

for f in $GRAPH_FILES; do
    gname="$(basename "$f" .txt)"
    header="$(head -1 "$f")"
    gv="$(echo "$header" | awk '{print $1}')"
    ge="$(echo "$header" | awk '{print $2}')"
    avg_deg="$(awk "BEGIN{printf \"%.1f\", 2*$ge/$gv}")"
    echo "| $gname | $gv | $ge | $avg_deg | SuiteSparse |" >> "$REPORT"
done

echo "" >> "$REPORT"

# Results table per graph – show plain vs greedy-md side by side
echo "## Results" >> "$REPORT"
echo "" >> "$REPORT"

for gname in $graphs; do
    echo "### $gname" >> "$REPORT"
    echo "" >> "$REPORT"

    graph_langs="$(grep ",$gname," "$CSV" | cut -d, -f3 | sort -u)"

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

            echo "| $alg | $gname | $lang | $plain_ms | $greedy_ms | ${speedup}× | $greedymd_ms | ${md_speedup}× | $greedymd_init | $greedymd_pct |" >> "$REPORT"
        done
    done
done

echo "" >> "$REPORT"

# Language speedup table (C++ = 1.0×)
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

        if [ "$rust_ratio" = "–" ]; then rust_fmt="–"; else rust_fmt="${rust_ratio}×"; fi
        if [ "$py_ratio" = "–" ]; then py_fmt="–"; else py_fmt="${py_ratio}×"; fi
        echo "| $alg | $gname | $cpp_ms | $rust_ms | $rust_fmt | $py_ms | $py_fmt |" >> "$REPORT"
    done
done

echo "" >> "$REPORT"
echo "---" >> "$REPORT"
echo "*Median of $RUNS runs. Wall-clock ms reported by each implementation. Timeout: ${TIMEOUT}s.*" >> "$REPORT"

echo ""
echo "  Run:         $RUN_ID"
echo "  Environment: $OUTDIR/environment.txt"
echo "  Report:      $OUTDIR/report.md"
echo "  Results:     $OUTDIR/results.csv"
echo "  Logs:        $OUTDIR/logs/"

# Create/update 'latest' symlink
BASEDIR="$(dirname "$OUTDIR")"
ln -sfn "$RUN_ID" "$BASEDIR/latest"
echo "  Latest:      $BASEDIR/latest -> $RUN_ID"
echo ""

# Append timing footer to environment.txt
write_env_footer

# Re-enable strict mode for the verdict
set -e

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
