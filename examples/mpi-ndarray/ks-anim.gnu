
set terminal {terminal}
set output '{output}'

set style line 1 lw 2 lc rgbcolor "blue"

unset key; unset xtics

set yrange [-3.5:3.5]
set xrange [0:255]

plot '-' with lines ls 1
{data}e
