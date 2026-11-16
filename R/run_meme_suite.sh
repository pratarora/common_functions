#!/bin/bash
################################################################################
# MEME Suite runner (Mercader Lab) — one FASTA in, one output folder out.
# Called by run_meme_suite() in reusable_functions/meme_functions.R; can also be
# run by hand. Nothing species- or project-specific is hard-coded here.
#
# Environment variables
#   FASTA        (required) primary sequences (e.g. DAR regions)
#   OUT_DIR      (required) output folder; one sub-folder per tool, plus RUN_STATUS
#   MOTIF_DBS    (required) motif databases (.meme), separated by ":"
#   NEG_FASTA    (optional) control sequences; xstreme/sea get --n
#   TOOLS        comma list of: xstreme, sea, centrimo, meme-chip   (default xstreme)
#   MEME_THREADS OpenMP threads for xstreme / meme-chip MEME step    (default 1)
#   SEED         random seed passed to every tool that takes one     (default 0)
#   MEME_ENV     conda env that holds MEME 5.5.x                    (default /opt/conda/envs/archr_tools)
#   XSTREME_ARGS extra xstreme arguments (default below = lab settings)
#
# Writes OUT_DIR/RUN_STATUS = "OK" or "ERROR <tool>" (also the exit status), so a
# caller that launched this in the background can wait for the file.
################################################################################
set -u -o pipefail

MEME_ENV="${MEME_ENV:-/opt/conda/envs/archr_tools}"
export PATH="${MEME_ENV}/bin:${MEME_ENV}/libexec/meme-5.5.9:${PATH}"
export LD_LIBRARY_PATH="${MEME_ENV}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"   # fimo links libxslt from the env (source build has no RPATH)

FASTA="${FASTA:-}"; OUT_DIR="${OUT_DIR:-}"; MOTIF_DBS="${MOTIF_DBS:-}"
NEG_FASTA="${NEG_FASTA:-}"; TOOLS="${TOOLS:-xstreme}"
MEME_THREADS="${MEME_THREADS:-1}"; SEED="${SEED:-0}"
XSTREME_ARGS="${XSTREME_ARGS:---dna --evt 0.05 --minw 6 --maxw 15 --align center --meme-mod zoops --streme-totallength 4000000 --meme-searchsize 100000 --sea-noseqs}"

fail() { echo "ERROR: $*" >&2; mkdir -p "${OUT_DIR:-.}"; echo "ERROR $*" > "${OUT_DIR:-.}/RUN_STATUS"; exit 1; }

[[ -n "$FASTA" && -n "$OUT_DIR" && -n "$MOTIF_DBS" ]] || fail "FASTA, OUT_DIR and MOTIF_DBS must be set"
[[ -s "$FASTA" ]] || fail "FASTA missing or empty: $FASTA"
[[ -z "$NEG_FASTA" || -s "$NEG_FASTA" ]] || fail "NEG_FASTA missing or empty: $NEG_FASTA"
command -v xstreme >/dev/null || fail "MEME not found under ${MEME_ENV}"

DB_ARGS=(); DB_ARGS_CHIP=(); DB_LIST=()
IFS=':' read -ra DB_LIST <<< "$MOTIF_DBS"
for db in "${DB_LIST[@]}"; do
    [[ -f "$db" ]] || fail "motif database not found: $db"
    DB_ARGS+=(--m "$db"); DB_ARGS_CHIP+=(--db "$db")
done

mkdir -p "$OUT_DIR"; rm -f "$OUT_DIR/RUN_STATUS"
exec > >(tee -a "$OUT_DIR/meme_suite_run.log") 2>&1
echo "== MEME Suite $(meme -version) | $(date) =="
echo "FASTA: $FASTA ($(grep -c '>' "$FASTA") sequences) | NEG: ${NEG_FASTA:-none} | tools: $TOOLS | threads: $MEME_THREADS | seed: $SEED"
printf 'DB: %s\n' "${DB_LIST[@]}"

status=OK
IFS=',' read -ra TOOL_LIST <<< "$TOOLS"
for tool in "${TOOL_LIST[@]}"; do
    echo "--- $tool ---"
    case "$tool" in
        xstreme)
            neg=(); [[ -n "$NEG_FASTA" ]] && neg=(--n "$NEG_FASTA")
            # shellcheck disable=SC2086
            xstreme --oc "$OUT_DIR/xstreme" --seed "$SEED" --meme-p "$MEME_THREADS" $XSTREME_ARGS \
                "${DB_ARGS[@]}" --p "$FASTA" "${neg[@]}" || { status="ERROR xstreme"; break; } ;;
        sea)
            neg=(); [[ -n "$NEG_FASTA" ]] && neg=(--n "$NEG_FASTA")
            sea --oc "$OUT_DIR/sea" --seed "$SEED" --align center "${DB_ARGS[@]}" \
                --p "$FASTA" "${neg[@]}" || { status="ERROR sea"; break; } ;;
        centrimo)
            centrimo --oc "$OUT_DIR/centrimo" --ethresh 10.0 "$FASTA" "${DB_LIST[@]}" \
                || { status="ERROR centrimo"; break; } ;;
        meme-chip)
            meme-chip --oc "$OUT_DIR/meme-chip" -seed "$SEED" --ccut 100 --dna --order 2 --minw 6 --maxw 15 \
                "${DB_ARGS_CHIP[@]}" --meme-mod zoops --meme-nmotifs 3 -meme-p "$MEME_THREADS" \
                --meme-searchsize 100000 --streme-pvt 0.05 --streme-totallength 4000000 \
                --centrimo-score 5.0 --centrimo-ethresh 10.0 "$FASTA" \
                || { status="ERROR meme-chip"; break; } ;;
        *) status="ERROR unknown tool $tool"; break ;;
    esac
done

echo "$status" > "$OUT_DIR/RUN_STATUS"
echo "== $status | $(date) =="
[[ "$status" == "OK" ]]
