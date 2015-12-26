# load(file = 'data.Rdata')

rm(list = ls())

start_time = proc.time()

library("Biobase")
# library("gtools")
library("clues")

source('gen_data_by_lnn_cluster_wise.r')
source('lnn.r')
source('rough_analyze.r')

G = 10000
n = 100
times = 1

delta_1 = 1
k_1 = 2
lambda_1 = 3
nu_1 = 4

delta_2 = 0.1
k_2 = 0.2
lambda_2 = 5
nu_2 = 6

k_3 = 0.1
lambda_3 = 2
nu_3 = 3


psi = c(delta_1, k_1, lambda_1, nu_1,
        delta_2, k_2, lambda_2, nu_2,
        0,       k_3, lambda_3, nu_3)
b = c(2,2,2)
t_pi = c(0.15, 0.15, 0.70)


t_pi_prior = c(0.33, 0.33, 0.34)
k_prior = 0.1


file_to_save_psi = 'psi.txt'
file_to_save_t_pi = 't_pi.txt'
file_to_save_accuracy = 'accuracy.txt'
file_to_save_notes = 'notes.txt'
file_to_save_adjusted_rand = 'adjusted_rand.txt'

# write to file_to_save_psi, file_to_save_t_pi, file_to_save_accuracy
write(toString(Sys.time()), file = file_to_save_t_pi)
write(toString(Sys.time()), file = file_to_save_psi)
write(toString(Sys.time()), file = file_to_save_accuracy)
write(c("total","npositive","npositive_recognize","npositive_not_recognize","false_positive","nnegative","nnegative_recognize","nnegative_not_recognize","false_negative"), file = file_to_save_accuracy, ncolumns = 9, append = TRUE, sep = '\t')



# write to file_to_save_notes
write(toString(Sys.time()), file = file_to_save_notes)
write(c('k_prior=', k_prior), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('pi_prior=', t_pi_prior), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write('TRUE VALUES:', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('delta_1=', delta_1), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('k_1=', k_1), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('lambda_1=', lambda_1), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('nu_1=', nu_1), file = file_to_save_notes, append = TRUE, sep = '\t')

write(c('delta_2=', delta_2), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('k_2=', k_2), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('lambda_2=', lambda_2), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('nu_2=', nu_2), file = file_to_save_notes, append = TRUE, sep = '\t')

write(c('k_3=', k_3), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('lambda_3=', lambda_3), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('nu_3=', nu_3), file = file_to_save_notes, append = TRUE, sep = '\t')


write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('t_pi=', t_pi), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('G=', G), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('n=', n), file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('times=', times), file = file_to_save_notes, append = TRUE, sep = '\t')

write('', file = file_to_save_notes, append = TRUE, sep = '\t')
write(c('b=', b), file = file_to_save_notes, append = TRUE, sep = '\t')



# write to file_to_save_adjusted_rand
write(toString(Sys.time()), file = file_to_save_adjusted_rand)
write(c('Rand','HA','MA','FM','Jaccard'), file = file_to_save_adjusted_rand, ncolumns = 5, append = TRUE, sep = '\t')

for (i in 1:times)
{
  print(c("now running loop:", i))

  E_Set = gen_data(G, n, psi, t_pi)
  save.image(file = 'data.Rdata')

  result = lnn(E_Set, t_pi_prior, k_prior, b, verbose = 1)

  write(c(result[[2]]$convergence,result[[2]]$par), file = file_to_save_psi, ncolumns = 13, append = TRUE, sep = '\t')
  write(result[[3]], file = file_to_save_t_pi, ncolumns = 3, append = TRUE, sep = '\t')

  rough_analysis = rough_analyze(as.matrix(fData(result[[1]])))
  write(rough_analysis, file = file_to_save_accuracy, ncolumns = length(rough_analysis), append = TRUE, sep = '\t')

  write(unname(adjustedRand(as.matrix(fData(result[[1]]))[,1],as.matrix(fData(result[[1]]))[,2])), file = file_to_save_adjusted_rand, ncolumns = 5, append = TRUE, sep = '\t')


  save.image(file = paste(c('simulation_',i,'.Rdata'), collapse = ''))
  
  print(proc.time()-start_time)
  #print(fData(result[[1]])[3098,])
  #print(sum(fData(result[[1]])[,3]))
}

end_time = proc.time()
print(end_time-start_time)
