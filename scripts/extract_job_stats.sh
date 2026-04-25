#!/usr/bin/env bash
#
# extract_job_stats.sh
#
# Print a pre-filled run_log.md entry for an SCC job. Pulls fields from:
#   - `qacct -j JOB_ID` (SGE accounting database — runtime, memory, exit code)
#   - `results/benchmarks/*.tsv` (Snakemake per-rule benchmarks)
#   - the run's pipeline log file
#
# Usage:
#   bash scripts/extract_job_stats.sh JOB_ID [LABEL]
#
# Example:
#   bash scripts/extract_job_stats.sh 4509000 "dedup-13M-first-attempt"
#
# Pipe the output into your run_log.md file with:
#   bash scripts/extract_job_stats.sh 4509000 >> run_log.md

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 JOB_ID [LABEL]" >&2
    exit 1
fi

JOB_ID="$1"
LABEL="${2:-unnamed}"
BENCHMARK_DIR="results/benchmarks"
LOG_FILE="logs/pipeline_${JOB_ID}.log"

# -----------------------------------------------------------------------------
# Pull SGE accounting fields. `qacct -j` gives one or more flat key/value
# blocks. We grep for the fields we want and strip whitespace.
# -----------------------------------------------------------------------------
qacct_field() {
    # $1 = field name (e.g. ru_wallclock)
    qacct -j "$JOB_ID" 2>/dev/null \
        | awk -v key="$1" '$1 == key { $1=""; sub(/^ +/,""); print; exit }'
}

EXIT_STATUS=$(qacct_field exit_status)
WALLCLOCK_SEC=$(qacct_field ru_wallclock)
CPU_SEC=$(qacct_field cpu)
MAXVMEM=$(qacct_field maxvmem)
HOSTNAME=$(qacct_field hostname)
QUEUE=$(qacct_field qname)
SLOTS=$(qacct_field slots)
QSUB_HOST=$(qacct_field qsub_host)
SUBMIT_TIME=$(qacct_field submission_time)
START_TIME=$(qacct_field start_time)
END_TIME=$(qacct_field end_time)

# Convert wallclock seconds to a friendly hh:mm:ss
if [[ -n "${WALLCLOCK_SEC:-}" ]]; then
    WALL_HMS=$(printf '%02d:%02d:%02d' \
        $((${WALLCLOCK_SEC%.*} / 3600)) \
        $(( (${WALLCLOCK_SEC%.*} % 3600) / 60 )) \
        $((${WALLCLOCK_SEC%.*} % 60)))
else
    WALL_HMS="(unknown — qacct may not have indexed this job yet)"
fi

# Determine pass/fail status
if [[ "${EXIT_STATUS:-}" == "0" ]]; then
    STATUS="Completed"
elif [[ -z "${EXIT_STATUS:-}" ]]; then
    STATUS="Unknown (qacct returned no record)"
else
    STATUS="Failed (exit $EXIT_STATUS)"
fi

# -----------------------------------------------------------------------------
# Emit the filled-in markdown block
# -----------------------------------------------------------------------------
cat <<EOF

### Run $(date +%F) — ${LABEL}

**SGE job ID:** \`${JOB_ID}\`
**Submitted by:** ${USER}
**Submission host:** ${QSUB_HOST:-?}
**Project:** ds596
**Submitted at:** ${SUBMIT_TIME:-?}
**Started at:** ${START_TIME:-?}
**Ended at:** ${END_TIME:-?}

#### Dataset

| Field | Value |
|---|---|
| Input file | _(fill in)_ |
| Source query | _(fill in)_ |
| Sequence count (input) | _(fill in)_ |
| Compressed size on disk | _(fill in)_ |
| Lectin anchors merged | _(fill in)_ |
| Deduplication step? | _(fill in)_ |
| Tagging scheme | _(fill in)_ |

#### Processing notes

_(fill in — what was different about this run, parameters, etc.)_

#### Resource request

| Field | Value |
|---|---|
| Cores | ${SLOTS:-?} |
| Compute node | ${HOSTNAME:-?} |
| Queue | ${QUEUE:-?} |

#### Job outcome

| Field | Value |
|---|---|
| Status | ${STATUS} |
| Exit status | ${EXIT_STATUS:-?} |
| Wall time used | ${WALL_HMS} (${WALLCLOCK_SEC:-?} s) |
| CPU time used | ${CPU_SEC:-?} s |
| Max RSS | ${MAXVMEM:-?} |
| Log file | \`${LOG_FILE}\` |

#### Per-rule statistics

_(From Snakemake \`benchmark:\` directives — auto-extracted below.)_

EOF

# -----------------------------------------------------------------------------
# Walk Snakemake benchmark TSVs and emit a markdown table row per rule.
# Benchmark TSV format (Snakemake 8): tab-separated, with header row.
# Columns include: s, h:m:s, max_rss, max_vms, max_uss, max_pss,
#                  io_in, io_out, mean_load, cpu_time
# -----------------------------------------------------------------------------
echo "| Rule | Wall time | Max RSS (MB) | Mean CPU load |"
echo "|---|---|---|---|"

if [[ -d "${BENCHMARK_DIR}" ]]; then
    for tsv in "${BENCHMARK_DIR}"/*.tsv; do
        [[ -f "$tsv" ]] || continue
        rule_name=$(basename "$tsv" .tsv)

        # Read the second line (data row) and pluck the right fields.
        # Awk indexes match the standard Snakemake benchmark column order.
        read -r seconds hms max_rss max_vms max_uss max_pss io_in io_out mean_load cpu_time \
            < <(awk 'NR==2 { for (i=1;i<=NF;i++) printf "%s ", \$i; print "" }' "$tsv")

        printf "| %s | %s | %s | %s |\n" \
            "$rule_name" "${hms:-?}" "${max_rss:-?}" "${mean_load:-?}"
    done
else
    echo "| _(no benchmark directory found at ${BENCHMARK_DIR})_ | | | |"
fi

cat <<EOF

#### Outputs

_(fill in observed sizes after the run)_

| Path | Size | Notes |
|---|---|---|
| \`data/query_rep50_rep_seq.fasta\` | | |
| \`data/query_all.fasta\` | | |
| \`results/all_vs_all_alignment.tsv\` | | |
| \`results/map_nodes.csv\` | | |
| \`results/map_links.csv\` | | |

#### Analysis observations

_(fill in)_

#### Next actions

- [ ] _(fill in)_

---
EOF