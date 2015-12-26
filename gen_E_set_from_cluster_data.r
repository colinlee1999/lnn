gen_E_set_from_cluster_data <- function(data, cluster)
{
	num_gene = nrow(data)
	
	feature_data = matrix(rep(0,3*num_gene),num_gene,3)
	rownames(feature_data) = rownames(data)

	feature_data[,1] = cluster

	return (ExpressionSet(assayData = data, featureData = new("AnnotatedDataFrame", data = data.frame(feature_data))))
}