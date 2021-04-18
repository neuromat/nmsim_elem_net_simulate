#! /bin/bash
# Last edited on 2019-05-28 16:29:01 by jstolfi

# Reads the trace of a neuron in an elem_level simulation.
# Plots the potential the neuron as a function of time.

tmin=$1; shift;    # First time to plot.
tmax=$1; shift;    # Last time to plot.
fname="$1"; shift  # Trace file name 

# END COMMAND LINE PARSING
# ----------------------------------------------------------------------

# Prefix for temporary file names
tmp="/tmp/$$"

tfile="${tmp}_I.png"
pfile="${fname%.*}_I.png"

# Extract the neuron number from the file name:
nnum="${fname%.txt}"
echo "nnum = ${nnum}"
nnum="${nnum#*_n}"
echo "nnum = ${nnum}"
nnum=$(( 0 + 10#${nnum} ))
echo "nnum = ${nnum}"

export GDFONTPATH="."

gnuplot <<EOF
  hpx = 1600
  vpx = 200
  tmin = ${tmin}
  tmax = ${tmax}
  tfile = "${tfile}"
  nnum = "${nnum}"
  vars = "external input {I}"
  load "plot_common.gpl"
  set yrange [-10:+10]
  set tmargin 0.5; set bmargin 0.5
  set ytics 5
  set mytics 5
  
  plot \
    "< egrep -e '^[ ]*[0-9]' ${fname}" using 1:7 title "I" with lines lt 1 lc rgb '#0022cc'
EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile}
  display ${pfile}
  rm ${tfile}
fi
