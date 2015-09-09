#!/bin/bash
R --vanilla --slave -f rp_loi.R --args ../../../tool-data/cri/survival_analysis/LUMINAL.RData genelist.txt 2 results.txt results.html results rfs er tamoxifen
