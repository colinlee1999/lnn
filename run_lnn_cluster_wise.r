source("run_limma_es.rmNA_GSE37020.r")
source("lnn.r")
source('gen_E_set_from_cluster_data.r')

positive_E_Set = gen_E_set_from_cluster_data(positive_data, 1)
negative_E_Set = gen_E_set_from_cluster_data(negative_data, 2)
neutral_E_Set = gen_E_set_from_cluster_data(neutral_data, 3)

positive_result = lnn(positive_E_Set, is_sim = 1, verbose = 1)
negative_result = lnn(negative_E_Set, is_sim = 1, verbose = 1)
neutral_result = lnn(neutral_E_Set, is_sim = 1, verbose = 1)

