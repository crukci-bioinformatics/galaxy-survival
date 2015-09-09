#!/bin/bash
R --vanilla --slave -f coxph_loi.R --args ../../../tool-data/cri/survival_analysis/LUMINAL.RData genelist.txt 1 results.txt rfs er tamoxifen
