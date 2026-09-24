#!/bin/bash
# run_facebook_batch.sh
#
# Chạy toàn bộ ma trận thực nghiệm RKCL + ATKC trên facebook, theo yêu cầu
# của cô Dung cho camera-ready CSONET 2026:
#   - dataset: facebook (fb.bin)
#   - ứng dụng: IM (Influence Maximization, oracle IC)
#   - Budget B: 50, 60, 70, 80, 90, 100 (tuyệt đối)
#   - cost node: random trong [10,20]  (dùng cho delta=0.1)
#                random trong [50,100] (dùng cho delta=0.5, xem ghi chú)
#   - delta: 0.1 và 0.5
#
# GHI CHÚ QUAN TRỌNG VỀ DELTA:
# Ban đầu dự kiến delta=0.1 và 0.5, nhưng cost[10,20]+B[50,100] khiến
# delta=0.5 vi phạm giả định granularity (cần c(e)>=25..50, vượt cost
# tối đa 20). Đã thử đổi cost lên [50,100] để delta=0.5 khả thi, nhưng
# khi đó q=floor(1/delta)=2 với MỌI delta trong (1/3, 0.5] (không chỉ
# delta=0.5) — nghĩa là "giảm delta xuống chút" quanh 0.45-0.49 không
# hề làm solution lớn hơn. Cô Dung đã xác nhận: dùng delta=0.33 (q=3)
# thay vì 0.5, với dataset cost[50,60] (hẹp hơn [50,100] theo chỉ đạo).
#
# Vậy cấu hình CUỐI CÙNG đã thống nhất với cô Dung:
#   delta=0.1  -> fb_cost1020.bin (cost~Uniform[10,20]), q=10
#   delta=0.33 -> fb_cost5060.bin (cost~Uniform[50,60]), q=3
#
# Tham số mc/scan_cap/chase_cap/max_scan/aug_scan_cap dựa trên đo thời
# gian thật trên máy 8-core của Nam Sơn:
#   RKCL full facebook: mc=50, scan_cap=50           -> ~2m40s/lần
#   ATKC full facebook: mc=10, chase_cap=10,
#                        max_scan=30, aug_scan_cap=30 -> ~3m8s/lần
# Các cap này là chệch thực dụng khỏi lý thuyết (đã ghi rõ trong code),
# không phải hành vi lý thuyết chuẩn của Algorithm 2/4.
#
# Tổng: 6 B x 2 delta x 2 algo = 24 lần chạy, ước tính ~70-90 phút.

set -u

BIN=./ic
COST_LOW=fb_cost1020.bin    # delta = 0.1
COST_HIGH=fb_cost4860.bin   # delta = 0.33 (cost floor 48 < B_min 50, tránh crash "no feasible element")
CSV=facebook_rkcl_atkc.csv

B_VALUES="50 60 70 80 90 100"
DELTAS="0.1 0.33"

export KIC_SEED=42
export KIC_LAMBDA=1.0

RKCL_MC=50
RKCL_SCAN_CAP=50

ATKC_MC=10
ATKC_CHASE_CAP=10
ATKC_MAX_SCAN=30
ATKC_AUG_SCAN_CAP=30

run_count=0
total=$(( $(echo $B_VALUES | wc -w) * $(echo $DELTAS | wc -w) * 2 ))

for delta in $DELTAS; do
    if [ "$delta" = "0.1" ]; then
        graph=$COST_LOW
    else
        graph=$COST_HIGH
    fi

    for B in $B_VALUES; do
        run_count=$((run_count+1))
        echo ""
        echo "=========================================="
        echo "[$run_count/$total] RKCL  graph=$graph B=$B delta=$delta"
        echo "=========================================="
        START=$(date +%s)
        $BIN --graph "$graph" --B_abs "$B" --alg rkcl --delta "$delta" \
             --mc $RKCL_MC --scan_cap $RKCL_SCAN_CAP \
             --consistency 1 --csv "$CSV"
        END=$(date +%s)
        echo "[RKCL] B=$B delta=$delta elapsed=$((END-START))s"

        run_count=$((run_count+1))
        echo ""
        echo "=========================================="
        echo "[$run_count/$total] ATKC  graph=$graph B=$B delta=$delta"
        echo "=========================================="
        START=$(date +%s)
        $BIN --graph "$graph" --B_abs "$B" --alg atkc --delta "$delta" --atkc_eps 0.5 \
             --mc $ATKC_MC --chase_cap $ATKC_CHASE_CAP \
             --max_scan $ATKC_MAX_SCAN --aug_scan_cap $ATKC_AUG_SCAN_CAP \
             --consistency 1 --csv "$CSV"
        END=$(date +%s)
        echo "[ATKC] B=$B delta=$delta elapsed=$((END-START))s"
    done
done

echo ""
echo "Done. Results appended to $CSV"
