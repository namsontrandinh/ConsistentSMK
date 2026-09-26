#!/bin/bash
# run_facebook_multistream.sh
#
# Chay MultiStream (baseline, Cui et al., trich dan truc tiep trong
# Related Work cua paper) tren facebook, chi de lam moc so sanh CHAT
# LUONG voi RKCL/ATKC — khong tinh consistency (MultiStream khong co
# dam bao consistency, va tinh full consistency se mat ~200+ gio do
# cach cai dat prefix-rerun khong co cap, xem giai thich trong chat).
#
# Chay tren CA HAI dataset cost da dung cho RKCL/ATKC de bieu do
# "value vs B" so sanh cong bang (cung phan bo cost):
#   fb_cost1020.bin  (cost~[10,20], dung cho delta=0.1 cua RKCL/ATKC)
#   fb_cost4860.bin  (cost~[48,60], dung cho delta=0.33 cua RKCL/ATKC)
#
# Tong: 2 dataset x 6 B = 12 lan chay, ~6 phut/lan (uoc luong tu do
# thuc te tren may 8-core) => ~70-75 phut.

set -u

BIN=./ic
CSV=facebook_multistream.csv

DATASETS="fb_cost1020.bin fb_cost4860.bin"
B_VALUES="50 60 70 80 90 100"

export KIC_SEED=42
export KIC_LAMBDA=1.0

MS_EPS=0.1
MS_MC=50

run_count=0
total=$(( $(echo $DATASETS | wc -w) * $(echo $B_VALUES | wc -w) ))

for graph in $DATASETS; do
    for B in $B_VALUES; do
        run_count=$((run_count+1))
        echo ""
        echo "=========================================="
        echo "[$run_count/$total] MultiStream  graph=$graph B=$B"
        echo "=========================================="
        START=$(date +%s)
        $BIN --graph "$graph" --B_abs "$B" --alg multistream --ms_eps $MS_EPS \
             --mc $MS_MC --consistency 0 --csv "$CSV"
        END=$(date +%s)
        echo "[MultiStream] graph=$graph B=$B elapsed=$((END-START))s"
    done
done

echo ""
echo "Done. Results in $CSV"
