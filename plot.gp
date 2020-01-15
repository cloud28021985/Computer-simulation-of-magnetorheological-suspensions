# 2019 Computer simulations of anisotropic structures in magnetorheological elastomers


set term pdfcairo \
font 'Liberation Serif, 12' \
solid \
size 8cm, 6cm


set o 'figs/sigma_vs_gamma.pdf'
set grid
unset key
set xlabel '{/:Italic γ}' font ',18'
set bmargin 4.0
set ylabel '{/:Italic σ}' font ',18'
plot \
'data/sigma_vs_gamma_experiment_1.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'red', \
'data/sigma_vs_gamma_experiment_2.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'blue', \
'data/sigma_vs_gamma_experiment_3.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'green', \
'data/sigma_vs_gamma_experiment_4.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'black'


set o 'figs/sigma_vs_gamma_sat.pdf'
set grid
unset key
set xlabel '{/:Italic γ}' font ',18'
set bmargin 4.0
set ylabel '{/:Italic σ_s}' font ',18'
plot \
'data/sigma_vs_gamma_experiment_sat_1.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'red', \
'data/sigma_vs_gamma_experiment_sat_2.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'blue', \
'data/sigma_vs_gamma_experiment_sat_3.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'green', \
'data/sigma_vs_gamma_experiment_sat_4.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'black'


set o 'figs/sigma_mean_vs_gamma.pdf'
set xlabel '{/:Italic γ}' font ',18'
set bmargin 4.0
set ylabel '{/:Italic σ}, kPa' font ',18'
plot [] [0.0:]\
'data/sigma_mean_vs_gamma_h_1.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'blue', \
'data/sigma_mean_vs_gamma_h_2.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'red'


set o 'figs/sigma_mean_vs_gamma_paper.pdf'
set xlabel '{/:Italic γ}' font ',18'
set bmargin 4.0
set ylabel '{/:Italic σ}, kPa' font ',18'
set label 1 '1' at 0.23, 1.25 font ',18'
set label 2 '2' at 0.23, 2.45 font ',18'
plot [] [0.0:2.6]\
'data/sigma_mean_vs_gamma_h_1.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'black', \
'data/sigma_mean_vs_gamma_h_2.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'black'


set o 'figs/sigma_mean_vs_gamma_sat.pdf'
unset label
set xlabel '{/:Italic γ}' font ',18'
set bmargin 4.0
set ylabel '{/:Italic σ_s}, kPa' font ',18'
plot [] [0.0:]\
'data/sigma_mean_vs_gamma_sat.txt' u 1:2 w lp pt 5 ps 0.2 lw 0.5 lt rgb 'black'




















