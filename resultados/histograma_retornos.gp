set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set grid
set key outside right top
set style line 1 lc rgb '#0060ad' lt 1 lw 2
set style line 2 lc rgb '#dd181f' lt 1 lw 2
set style line 3 lc rgb '#00a000' lt 1 lw 2
set style line 4 lc rgb '#ffa000' lt 1 lw 2
set output 'resultados/histograma_retornos.png'
set title 'Distribuição dos Retornos Diários'
set xlabel 'Valor'
set ylabel 'Frequência'
set style fill solid 0.5
set boxwidth 0.00110168
plot 'resultados/histograma_retornos.dat' using 1:2 with boxes ls 1 title 'Distribuição'
