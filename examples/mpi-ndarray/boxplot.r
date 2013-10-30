library(ggplot2)

# read in send, sweep, recv level timings and generate plots
base = 'timings_heat_p64l3_recv'
t = read.csv(paste(base, '.csv', sep=''))
plt = ggplot(t, aes(factor(trial), log(delta))) + geom_boxplot(aes(fill=factor(level))) + labs(x='trial', y='log(elapsed time [s])')
ggsave(paste(base, '.pdf', sep=''), plt)

#d1 = read.csv("edison.r1_burgers_p04_l3.csv")
#d2 = read.csv("sweep_edison.r1_burgers_p04_l3.csv")
#p = ggplot(d1, aes(factor(trial), log(send))) + geom_boxplot(aes(fill=factor(level)))
#ggplot(d1) + geom_boxplot(aes(factor(trial), log(send))) + geom_boxplot(aes(fill=factor(level)))
#p + geom_boxplot(data=d2, aes(factor(trial), log(sweep))) + geom_boxplot(aes(fill=factor(level)))


ggplot(NULL, aes(factor(trial), log(delta))) + geom_boxplot(data=d1, aes(fill=factor(level))) + geom_boxplot(data=d2, aes(fill=factor(level)))



# read in speedups.csv and produce speedup/efficiency boxplots
# XXX: should add in line plot with theoretical speedup

s = read.csv('speedups.csv')
for (p in unique(s$prob)) {
  h = s[which(s$prob == p),]
  plt = ggplot(h, aes(factor(nproc), speedup)) + geom_boxplot() + labs(x='no. of processors')
  ggsave(paste(p, "_speedup.pdf", sep=""), plt)
  plt = ggplot(h, aes(factor(nproc), efficiency)) + geom_boxplot() + labs(x='no. of processors')
  ggsave(paste(p, "_efficiency.pdf", sep=""), plt)
}
