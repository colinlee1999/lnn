library(iCheck)
load("simulation_100.Rdata")

# generate simulated data set from conditional normal distribution
set.seed(1234567)
es.sim = genSimData.BayesNormal(nCpGs = 100,
                                nCases = 20, nControls = 20,
                                mu.n = -2, mu.c = 2,
                                d0 = 20, s02 = 0.64, s02.c = 1.5, testPara = "var",
                                outlierFlag = FALSE,
                                eps = 1.0e-3, applier = lapply)
print(es.sim)
# although the generated data is not from
# paired design, we use it to illusrate the
# usage of the function lmFitPaired
a = result[[1]]
res.limma = lmFitPaired(
  es = a,
  formula = ~1,
  pos.var.interest = 0,
  pvalAdjMethod = "fdr",
  alpha = 0.05,
  probeID.var = "X1",
  gene.var = "X2",
  chr.var = "X3",
  verbose = TRUE)

qwl=res.limma$frame.unsorted
hist(-log10(qwl$pval))
obsmem=qwl$pval<1.0e-180

fDat=fData(a)
truemem=fDat$X1
pos1=which(truemem==1)[1]
pos2=which(truemem==2)[1]
pos3=which(truemem==3)[1]

dat=exprs(a)

par(mfrow=c(3,1))
hist(dat[pos1,],main="overexprs")
hist(dat[pos2,],main="underexprs")
hist(dat[pos3,],main="non")

print(mean(dat[pos1,]))
print(mean(dat[pos2,]))
print(mean(dat[pos3,]))