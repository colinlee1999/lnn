library(iCheck)
source("gen_new_expression_set_from_es.rmNA_GSE37020.r")
E_Set = gen_new_expression_set(es.rmNA)

res.limma = lmFitPaired(E_Set, probeID.var = "ID", gene.var = "Name", chr.var = "IlmnStrand")

result.limma = res.limma$frame.unsorted

significant = result.limma[result.limma$pval<0.05,]

positive_significant = significant[significant$stats>0,]
negative_significant = significant[significant$stats<0,]

write.csv(positive_significant$probeIDs, 'positive_probeID_by_limma.csv')
write.csv(negative_significant$probeIDs, 'negative_probeID_by_limma.csv')

positive_probeIDs = levels(positive_significant$probeIDs)[c(positive_significant$probeIDs)]
negative_probeIDs = levels(negative_significant$probeIDs)[c(negative_significant$probeIDs)]

data = assayData(E_Set)$exprs

positive_data = data[positive_probeIDs,]
negative_data = data[negative_probeIDs,]
neutral_data = data[rownames(data)[(!rownames(data) %in% c(positive_probeIDs,negative_probeIDs))],]


# hist of tau_g by cluster
positive_data_tau_g = 1 / (apply(positive_data, 1, mad)^2)
negative_data_tau_g = 1 / (apply(negative_data, 1, mad)^2)
neutral_data_tau_g = 1 / (apply(neutral_data, 1, mad)^2)

# par(mfrow=c(3,1))
# hist(log10(positive_data_tau_g), xlim=c(0,5))
# hist(log10(negative_data_tau_g), xlim=c(0,5))
# hist(log10(neutral_data_tau_g), xlim=c(0,5))

# hist of k by cluster
positive_data_mu_g = apply(positive_data, 1, median)
negative_data_mu_g = apply(negative_data, 1, median)
neutral_data_mu_g = apply(neutral_data, 1, median)

positive_data_mu_g_var = mad(positive_data_mu_g)^2
negative_data_mu_g_var = mad(negative_data_mu_g)^2
neutral_data_mu_g_var = mad(negative_data_mu_g)^2

positive_data_k = positive_data_mu_g_var * positive_data_tau_g
negative_data_k = negative_data_mu_g_var * negative_data_tau_g
neutral_data_k = neutral_data_mu_g_var * neutral_data_tau_g

# hist(log10(positive_data_k), xlim=c(0,3))
# hist(log10(negative_data_k), xlim=c(0,3))
# hist(log10(neutral_data_k), xlim=c(0,3))

num_positive = nrow(positive_data)
num_negative = nrow(negative_data)
num_neutral = nrow(neutral_data)

par(mfrow=c(1,1))

ttframe = data.frame(mu = c(positive_data_mu_g, negative_data_mu_g, neutral_data_mu_g),
	mem = c(rep(1, num_positive), rep(2, num_negative), rep(3, num_neutral)))

boxplot(mu~mem, data = ttframe)

ttframe = data.frame(mu = log10(c(1/positive_data_tau_g, 1/negative_data_tau_g, 1/neutral_data_tau_g)),
	mem = c(rep(1, num_positive), rep(2, num_negative), rep(3, num_neutral)))

boxplot(mu~mem, data = ttframe)


















# fDat=fData(E_Set)
# truemem=fDat$X1
# pos1=which(truemem==1)[1]
# pos2=which(truemem==2)[1]
# pos3=which(truemem==3)[1]

# dat=exprs(E_Set)

# par(mfrow=c(1,1))
# hist(dat[pos1,],main="overexprs")
# hist(dat[pos2,],main="underexprs")
# hist(dat[pos3,],main="non")

# print(mean(dat[pos1,]))
# print(mean(dat[pos2,]))
# print(mean(dat[pos3,]))