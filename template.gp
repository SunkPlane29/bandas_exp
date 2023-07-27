datafile = "EoS-insensitive_posterior_samples.dat"
set terminal png size 400,300 enhanced font "Helvetica,10"
set output "bandas_exp_MR1.png"
plot [7:15][0.5:3] datafile u 5:1 with points

set output "bandas_exp_MR2.png"
plot [7:15][0.5:3] datafile u 6:2 with points

set output "bandas_exp_MΛ1.png"
plot datafile u 1:3 with points

set output "bandas_exp_MΛ2.png"
plot datafile u 2:4 with points
