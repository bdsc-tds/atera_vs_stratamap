# Shared config + SLURM dispatch helper, sourced by 02_preprocess.sh and 03_generate_figures.sh.
# Every compute step in this package - including "cheap" ones like plotting - goes through
# dispatch() so USE_SLURM=true means NOTHING runs on whatever node you launched the script from,
# not even a few-second ggplot call.

export REPRO_ROOT="${REPRO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export REPRO_DATA="${REPRO_DATA:-$REPRO_ROOT/data}"
export REPRO_RESULTS="${REPRO_RESULTS:-$REPRO_ROOT/results}"

USE_SLURM="${USE_SLURM:-false}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
SLURM_PARTITION="${SLURM_PARTITION:-cpu}"
if [ "$USE_SLURM" = "true" ] && [ -z "$SLURM_ACCOUNT" ]; then
  echo "USE_SLURM=true requires SLURM_ACCOUNT to be set" >&2
  exit 1
fi

R_CMD="${R_CMD:-Rscript --vanilla}"
PYTHON_CMD="${PYTHON_CMD:-python3}"

CODE_R="$REPRO_ROOT/code/R"
CODE_PY="$REPRO_ROOT/code/python"
LOGDIR="$REPRO_RESULTS/logs"
mkdir -p "$LOGDIR"

# dispatch NAME MEM TIME CPUS "command string" - every command in this package goes through this,
# never called directly - runs under sbatch (backgrounded, --wait) when USE_SLURM=true, or
# directly (sequentially, so memory-hungry steps don't pile up on one machine) otherwise.
dispatch() {
  local name="$1" mem="$2" time="$3" cpus="$4" cmd="$5"
  if [ "$USE_SLURM" = "true" ]; then
    sbatch --wait --parsable --job-name="$name" --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
      --time="$time" --mem="$mem" --cpus-per-task="$cpus" \
      --output="$LOGDIR/${name}_%j.out" --error="$LOGDIR/${name}_%j.err" --wrap="$cmd" &
  else
    echo "=== $name ==="
    eval "$cmd"
  fi
}
# Blocks until every dispatch() call issued since the last wait_all finishes; a no-op locally
# (dispatch already ran synchronously in that case).
wait_all() { [ "$USE_SLURM" = "true" ] && wait; return 0; }
