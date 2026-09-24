#!/bin/bash
# run_facebook_delta033_fix.sh
#
# Chạy lại TOÀN BỘ 12 dòng (6 B x 2 algo) cho delta=0.33, thống nhất
# dùng fb_cost4860.bin (thay vì lẫn lộn fb_cost5060.bin cho B=60-100
# và fb_cost4860.bin cho B=50 như trước). Ghi vào file CSV riêng để
# không lẫn với facebook_rkcl_atkc.csv cũ — sau khi chạy xong, xoá các
# dòng delta=0.33 cũ (graph=fb_cost5060.bin hoặc fb_cost50100.bin)
# trong facebook_rkcl_atkc.csv rồi nối file này vào.

set -u

BIN=./ic
GRAPH=fb_cost4860.bin
CSV=facebook_delta033_v2.csv

B_VALUES="50 60 70 80 90 100"
DELTA=0.33

export KIC_SEED=42
export KIC_LAMBDA=1.0

RKCL_MC=50
RKCL_SCAN_CAP=50

ATKC_MC=10
ATKC_CHASE_CAP=10
ATKC_MAX_SCAN=30
ATKC_AUG_SCAN_CAP=30

run_count=0
total=12

for B in $B_VALUES; do
    run_count=$((run_count+1))
    echo ""
    echo "=========================================="
    echo "[$run_count/$total] RKCL  graph=$GRAPH B=$B delta=$DELTA"
    echo "=========================================="
    START=$(date +%s)
    $BIN --graph "$GRAPH" --B_abs "$B" --alg rkcl --delta "$DELTA" \
         --mc $RKCL_MC --scan_cap $RKCL_SCAN_CAP \
         --consistency 1 --csv "$CSV"
    END=$(date +%s)
    echo "[RKCL] B=$B delta=$DELTA elapsed=$((END-START))s"

    run_count=$((run_count+1))
    echo ""
    echo "=========================================="
    echo "[$run_count/$total] ATKC  graph=$GRAPH B=$B delta=$DELTA"
    echo "=========================================="
    START=$(date +%s)
    $BIN --graph "$GRAPH" --B_abs "$B" --alg atkc --delta "$DELTA" --atkc_eps 0.5 \
         --mc $ATKC_MC --chase_cap $ATKC_CHASE_CAP \
         --max_scan $ATKC_MAX_SCAN --aug_scan_cap $ATKC_AUG_SCAN_CAP \
         --consistency 1 --csv "$CSV"
    END=$(date +%s)
    echo "[ATKC] B=$B delta=$DELTA elapsed=$((END-START))s"
done

echo ""
echo "Done. Results in $CSV"
echo "Buoc tiep theo: xoa cac dong delta=0.33 cu (graph chua 'fb_cost5060' hoac 'fb_cost50100')"
echo "trong facebook_rkcl_atkc.csv, roi noi $CSV vao."
