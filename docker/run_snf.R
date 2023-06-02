library(SNFtool)

args<-commandArgs(TRUE)

argnum <- length(args)

# need the hyperparams AND at least two input files
if(argnum < 6){
    message('SNF runs on a minimum of two input data types. Please add an additional files.')
    quit(status=1)
}

K<-as.integer(args[1])        # number of neighbors, usually (10~30)
alpha<-as.numeric(args[2])    # hyperparameter, usually (0.3~0.8)
T <- as.numeric(args[3])      # Number of Iterations, usually (10~20)
num_clusters <- as.integer(args[4])
input_files <- args[5:argnum]

# change the working directory to co-locate with the input files:
working_dir <- dirname(input_files[[1]])
setwd(working_dir)

mtx_list <- lapply(input_files, read.table)

# In the normalization and distance calculations, it is expected that the matrix is (n_samples, n_features) so we need
# to first transpose the canonical orientation of our input matrices, which are (n_genes, n_samples)
mtx_list <- lapply(mtx_list, t)

# It's recommended to normalize. This function performs standardization column-wise. Since we have
# transposed the inputs, this normalizes the expressions for each probe/gene/feature.
mtx_list <- lapply(mtx_list, standardNormalization)

## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; if the data is discrete, we recommend the users to use ""
dist_list <- lapply(mtx_list, function(x) dist2(as.matrix(x),as.matrix(x)))

## next, construct similarity graphs
affinity_list <- lapply(dist_list, function(x) affinityMatrix(x, K, alpha))

W = SNF(affinity_list, K, T)

# perform clustering on the fused network.
clustering = spectralClustering(W, num_clusters);

similarity_mtx_output_file = paste(working_dir,'snf_similarities.tsv', sep='/')
clustering_output_file = paste(working_dir,'snf_clusters.tsv', sep='/')
write.table(W, similarity_mtx_output_file, sep='\t', row.names=F, col.names=F)
write.table(clustering, clustering_output_file, sep='\t', row.names=F, col.names=F)

# create the expected outputs file:
json_str = paste0(
       '{"snf_similarity":"', similarity_mtx_output_file, '",',
       '"snf_clustering":"', clustering_output_file, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)