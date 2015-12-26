rm(list = ls())
load("es.rmNA_GSE37020.rda")

start_time = proc.time()

source('lnn.r')
source("gen_new_expression_set_from_es.rmNA_GSE37020.r")

E_Set = gen_new_expression_set(es.rmNA)

file_to_save_psi = 'psi.txt'
file_to_save_t_pi = 't_pi.txt'
file_to_save_cluster_result = 'cluster_result.csv'

write(toString(Sys.time()), file = file_to_save_t_pi)
write(toString(Sys.time()), file = file_to_save_psi)
write(toString(Sys.time()), file = file_to_save_cluster_result)

print("OK")
result.lnn = lnn(E_Set, plot = 1, verbose = 1)

write(c(result.lnn[[2]]$convergence,result.lnn[[2]]$par), file = file_to_save_psi, ncolumns = 6, append = TRUE, sep = '\t')
write(result.lnn[[3]], file = file_to_save_t_pi, ncolumns = 3, append = TRUE, sep = '\t')
write.csv(fData(result.lnn[[1]]), file = 'cluster_result.csv')
save.image(file = 'analysis_of_es.rmNA_GSE37020.Rdata')

print(proc.time()-start_time)
