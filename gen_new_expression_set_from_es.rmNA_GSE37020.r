gen_new_expression_set <- function(es.rmNA)
{
	load("es.rmNA_GSE37020.rda")
	require(data.table)
	pDat = as.matrix(data.table(pData(es.rmNA), key = "age"))

	str1 = "Cyto+ve"
	str2 = "Cyto-ve"

	group1_names = pDat[pDat[,39] == str1, 2]
	group2_names = pDat[pDat[,39] == str2, 2]

	group1_data = assayData(es.rmNA)[['exprs']][,group1_names]
	group2_data = assayData(es.rmNA)[['exprs']][,group2_names]

	dat = group1_data - group2_data
	# dat = log(group1_data) - log(group2_data)
	rownames(dat) = rownames(assayData(es.rmNA)[['exprs']])

	# dat = - dat

	return(ExpressionSet(assayData = dat, featureData = featureData(es.rmNA)))
}