set terminal qt size 900,1000 enhanced font "Arial,11"

MASS  = 1.0
K     = 1.0
OMEGA = sqrt(K / MASS)
X0    = 1.0
V0    = 0.0

E(x,v) = 0.5 * MASS * v**2 + 0.5 * K * x**2
E0     = E(X0, V0)

X_exact(t) = X0 * cos(OMEGA * t) + (V0 / OMEGA) * sin(OMEGA * t)
V_exact(t) = -X0 * OMEGA * sin(OMEGA * t) + V0 * cos(OMEGA * t)

set grid
set xlabel "t"
set multiplot layout 3,1 title "Гармонический осциллятор"

set ylabel "x(t)"
set title "Положение x(t)"
plot 'oscdata' using 1:2 with lines lw 2 lc rgb "#00aa00" title "RK4", \
     X_exact(x) with lines lw 2 lc rgb "black" dashtype 2 title "точное x(t)"

set ylabel "E(t)"
set title "Полная энергия E(t)"
plot 'oscdata' using 1:(E($2,$3)) with lines lw 2 lc rgb "#00aa00" title "E(t) числ.", \
     E0 with lines lw 2 lc rgb "black" dashtype 2 title sprintf("E точная = %.6g", E0)

set ylabel "x - x_exact(t)"
set title "Ошибка x(t) - x_exact(t)"
plot 'oscdata' using 1:($2 - X_exact($1)) with lines lw 2 lc rgb "red" title "ошибка", \
     0 with lines lw 1 lc rgb "black" dashtype 2 notitle

unset multiplot

pause -1 "Нажмите Enter, чтобы закрыть окно"
