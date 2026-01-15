# !/bin/bash

bin_edges=( 0.0,1.0
            0.0,0.975,1.0
            0.0,0.876,0.9954,1.0
            0.0,0.7452,0.975,0.9982,1.0
            0.0,0.6386,0.9338,0.9912,0.9988,1.0
            0.0,0.5558,0.876,0.975,0.9954,0.9992,1.0
            0.0,0.4932,0.8094,0.95,0.9882,0.9972,0.9994,1.0
            0.0,0.4392,0.7454,0.9164,0.9752,0.9934,0.9982,0.9996,1.0
            0.0,0.3984,0.6832,0.876,0.9572,0.986,0.9956,0.9986,0.9996,1.0
            0.0,0.3598,0.6386,0.832,0.9338,0.9752,0.9912,0.9968,0.9988,0.9996,1.0 )

counter=1

for bin_edge in "${bin_edges[@]}"; do
    echo ${counter} >> /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/eTau_2022/res2b/res2b.txt

    python3 bin_opt/rebinAndRunLimits.py --input /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/hh_res2b_eTau_2022_13p6TeV.txt --output /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/eTau_2022/res2b/num_bins_${counter} --bin-edges ${bin_edge} --poi r #2>&1 | grep "Expected 95% CL limit:"  >> /eos/user/t/toakhter/bin_opt_tests/bin_opt_test_oct_hamburg_v2/equal_integral/res2b/res2b.txt

    counter=$((counter + 1))
done 
