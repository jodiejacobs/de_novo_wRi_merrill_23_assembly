set title 'read length'
set xlabel 'read length, bin width = 250'
set ylabel 'number of reads'

binwidth=250
set boxwidth binwidth
bin(x,width) = width*floor(x/width) + binwidth/2.0

set terminal png size 1024,1024
set output './canu.seqStore/readlengths-cor.png'
plot [] './canu.seqStore/readlengths-cor.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes title ''
