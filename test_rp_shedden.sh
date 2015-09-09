#!/bin/bash
R --vanilla --slave -f rp_shedden.R --args ../../../tool-data/cri/survival_analysis/Shedden3.RData genelist.txt 2 results.txt results.html results rfs kra.pos
