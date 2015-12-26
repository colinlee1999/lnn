library("Biobase")
library("gtools")

source('gen_data_by_lnn.r')
source('lnn.r')

N = 1000
n = 100

delta_1 = 1
delta_2 = 1
k = 3
lambda = 2
nu = 1

psi = c(delta_1, delta_2, k, lambda, nu)

theta = c(0.1, 0.1, 0.8)

gen = 1
if (gen)
{
  E_Set = gen_data(N, n, psi, theta)
  save.image(file = 'data.Rdata')
} else
{
  load(file = 'data.Rdata')
}

b = c(2,2,2)
result = lnn(E_Set, b, 1)
print(result[3])
# print(E_Set)
# print(mleinfo)
#print (pData(featureData(result)))
