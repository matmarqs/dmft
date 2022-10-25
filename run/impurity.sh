#!/bin/sh

[ ! -d "./out" ] && mkdir -p "./out"    # create ./out, if does not exists
"./nca-latest" "out=out" "Sig=in/Sigmapp.inp" "Ac=in/Ac-last-0.2.inp" "cix=cix-latest.dat" \
    "U=4.0" "T=0.1" "Ed=-2.0"
#"./nca-2004" "out=out" "Sig=in/Sigmapp2004.inp" "Ac=in/Ac.inp" "cix=cix-2004.dat" \
#    "U=4.0" "T=0.1" "Ed={-2.0}"
#"./nca-2004" "out=out" "Sig=in/sig-SIAM.inp" "Ac=in/Ac-SIAM-0.02.inp" "cix=cix-2004.dat" \
#    "U=4" "T=0.1" "Ed=-2.0"
#"./nca-2004" "out=out" "Sig=in/ex-Sigma.000" "Ac=in/ex-Ac.12" "cix=in/ex-cix4567.cix" \
#    "U=4" "T=0.1" "Ed={-20.67, -20.58, -20.64}"

[ -f "./cores.dat" ] && rm "./cores.dat"
