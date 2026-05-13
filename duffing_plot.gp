set terminal pngcairo size 1200, 800 enhanced font "Times,12"
set output "duffing_plots.png"

set multiplot layout 2,2

set title "Смещение x(t)"
set xlabel "Время t"
set ylabel "Смещение x"
plot "duffing_data.txt" using 1:2 with lines linewidth 2 linecolor rgb "blue" notitle

set title "Фазовый портрет"
set xlabel "Смещение x"
set ylabel "Скорость v"
set size square
plot "duffing_data.txt" using 2:3 with lines linewidth 1 linecolor rgb "red" notitle

set title "Относительная ошибка энергии (доказательство точности)"
set xlabel "Время t"
set ylabel "Относительная ошибка"
set logscale y
plot "duffing_data.txt" using 1:5 with lines linewidth 1 linecolor rgb "green" notitle

set title "Полная механическая энергия"
set xlabel "Время t"
set ylabel "Энергия"
plot "duffing_data.txt" using 1:4 with lines linewidth 1 linecolor rgb "purple" notitle

unset multiplot
set output
