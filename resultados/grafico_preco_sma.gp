set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set grid
set key outside right top
set style line 1 lc rgb '#0060ad' lt 1 lw 2
set style line 2 lc rgb '#dd181f' lt 1 lw 2
set style line 3 lc rgb '#00a000' lt 1 lw 2
set style line 4 lc rgb '#ffa000' lt 1 lw 2
set output 'resultados/grafico_preco_sma.png'
set title 'Preço do Ouro com Média Móvel'
set xlabel 'Tempo'
set ylabel 'Valor'
plot 'resultados/grafico_preco_sma_0.dat' using 1:2 with lines ls 1 title 'Preço', 'resultados/grafico_preco_sma_1.dat' using 1:2 with lines ls 2 title 'SMA(20)'
