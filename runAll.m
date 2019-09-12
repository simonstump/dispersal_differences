%This code runs all of the simulations for you, and then makes all of the figures.

run_holdAbund(1)
run_holdAbund(2)
run_holdAbund(3)
run_holdAbund(11)
run_invasion(1,0)
run_invasion(2,0)
run_invasion(3,0)
run_invasion(4,0)
run_invasion(6,0)
run_invasion(7,0)
run_approx_check(1)
run_approx_check(2)
run_approx_check(3)
run_invasion(1,1)
run_invasion(2,1)
run_invasion(3,1)
run_invasion(4,1)
run_invasion(6,1)
run_invasion(7,1)

clc, clear
run_for_distribution

clc, clear
for_emergent_tradeoff


clc, cleardraw_fig2bcclc, clear
draw_fig3
clc, clear
draw_figA1
clc, clear
draw_figA2clc, clear
draw_figA3A4

clc, clear
draw_figA5