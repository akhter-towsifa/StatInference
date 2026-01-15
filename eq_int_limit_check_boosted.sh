# !/bin/bash

bin_edges=( 0.0,1.0
            0.0,0.9528,1.0
)

counter=1

for bin_edge in "${bin_edges[@]}"; do
    echo ${counter} >> /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/muTau_2022/boosted/boosted.txt

    python3 bin_opt/rebinAndRunLimits.py --input /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/hh_boosted_muTau_2022_13p6TeV.txt --output /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/muTau_2022/boosted/num_bins_${counter} --bin-edges ${bin_edge} --poi r # 2>&1 | grep "Expected 95% CL limit:"  >> /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/eTau_2022/boosted/boosted.txt

    counter=$((counter + 1))
done 
