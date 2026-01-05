set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set grid
set key outside right top
set style line 1 lc rgb '#0060ad' lt 1 lw 2
set style line 2 lc rgb '#dd181f' lt 1 lw 2
set style line 3 lc rgb '#00a000' lt 1 lw 2
set style line 4 lc rgb '#ffa000' lt 1 lw 2
set output 'resultados/grafico_acf.png'
set title 'Função de Autocorrelação (ACF)'
set xlabel 'Lag'
set ylabel 'Autocorrelação'
set style fill solid 0.5
set boxwidth 0.8
set yrange [-1:1]
plot 'resultados/grafico_acf.dat' using 1:2 with impulses lw 3 ls 1 title 'ACF', 0.590962 with lines ls 2 title 'IC 95%', -0.590962 with lines ls 2 notitle
