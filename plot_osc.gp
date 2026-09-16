set terminal qt size 900,1000 enhanced font "Arial,11"

E(x,v) = 0.5*(v**2 + x**2)
E0     = 0.5

set grid
set xlabel "t"
set multiplot layout 3,1 title "Гармонический осциллятор"

set ylabel "x(t)"
set title "Положение x(t)"
plot 'oscdata' using 1:2 with lines lw 2 lc rgb "#00aa00" title "RK4", \
     cos(x)     with lines lw 2 lc rgb "black" dashtype 2 title "cos(t)"

set ylabel "E(t)"
set title "Полная энергия E(t)"
plot 'oscdata' using 1:(E($2,$3)) with lines lw 2 lc rgb "#00aa00" title "E(t) числ.", \
     E0 with lines lw 2 lc rgb "black" dashtype 2 title "E точная = 0.5"

set ylabel "x - cos(t)"
set title "Ошибка x(t) - cos(t)"
plot 'oscdata' using 1:($2 - cos($1)) with lines lw 2 lc rgb "red" title "ошибка", \
     0 with lines lw 1 lc rgb "black" dashtype 2 notitle

unset multiplot

pause -1 "Нажмите Enter, чтобы закрыть окно"
