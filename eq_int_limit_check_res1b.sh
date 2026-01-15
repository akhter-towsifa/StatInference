# !/bin/bash

bin_edges=( 0.0,1.0
            0.0,0.8168,1.0
            0.0,0.582,0.9368,1.0
            0.0,0.4322,0.8168,0.9666,1.0
            0.0,0.336,0.6892,0.9012,0.978,1.0
            0.0,0.281,0.5824,0.817,0.9368,0.984,1.0
            0.0,0.2448,0.5014,0.733,0.88,0.9554,0.9872,1.0
            0.0,0.215,0.4326,0.6502,0.8172,0.9166,0.9666,0.9894,1.0
            0.0,0.1922,0.3778,0.583,0.7524,0.868,0.9372,0.9736,0.991,1.0
            0.0,0.1736,0.336,0.5262,0.6894,0.8172,0.9014,0.9506,0.9784,0.9922,1.0
)

counter=1

for bin_edge in "${bin_edges[@]}"; do
    echo ${counter} >> /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/eTau_2022/res1b/res1b.txt

    python3 bin_opt/rebinAndRunLimits.py --input /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/hh_res1b_eTau_2022_13p6TeV.txt --output /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/eTau_2022/res1b/num_bins_${counter} --bin-edges ${bin_edge} --poi r # 2>&1 | grep "Expected 95% CL limit:"  >> /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/eTau_2022/res1b/res1b.txt

    counter=$((counter + 1))
done 
