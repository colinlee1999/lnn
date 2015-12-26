source("run_lnn_es.rmNA_GSE37020.r")
source("run_limma_es.rmNA_GSE37020.r")

recognized.lnn = fData(result.lnn[[1]])

table(levels(recognized.lnn[recognized.lnn[,'Name'] == 1, 'ID'])[c(recognized.lnn[recognized.lnn[,'Name'] == 1, 'ID'])] %in% positive_probeIDs)
table(levels(recognized.lnn[recognized.lnn[,'Name'] == 2, 'ID'])[c(recognized.lnn[recognized.lnn[,'Name'] == 2, 'ID'])] %in% negative_probeIDs)
table(levels(recognized.lnn[recognized.lnn[,'Name'] == 3, 'ID'])[c(recognized.lnn[recognized.lnn[,'Name'] == 3, 'ID'])] %in% neutral_probeIDs)
