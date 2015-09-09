#!/bin/bash
#R --vanilla --slave -f rp.R --args loi /data_ls/galaxy/data/cri/survival_analysis/LUMINAL.RData genelist.txt 2 results.txt results.html results rfs er tamoxifen 0 0 0 0 0 0 0 0

sh r_wrapper.sh rp.R loi /data_ls/galaxy/data/cri/survival_analysis/LUMINAL.RData genelist.txt 2 test_results.txt test_results.html test_results rfs all tamoxifen none none none none none none none none
