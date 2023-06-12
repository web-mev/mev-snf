library(SNFtool)

args<-commandArgs(TRUE)

argnum <- length(args)

# need the hyperparams AND at least two input files
if(argnum < 6){
    message('SNF runs on a minimum of two input data types. Please add an additional files.')
    quit(status=1)
}

K<-as.integer(args[1])               # number of neighbors, usually (10~30)
alpha<-as.numeric(args[2])           # hyperparameter, usually (0.3~0.8)
T <- as.numeric(args[3])             # Number of Iterations, usually (10~20)
num_clusters <- as.integer(args[4])
primary_input_file <- args[5]        # the matrix which dictates the final sample/aliquot IDs
input_files <- args[6:argnum]        # other matrices to use

# change the working directory to co-locate with the input files:
working_dir <- dirname(primary_input_file)
setwd(working_dir)

# read the primary file and extract the column names which will dictate the final sample IDs
primary_mtx = read.table(primary_input_file, sep='\t', row.names=1, header=T, check.names=F)
sample_ids = colnames(primary_mtx)

# check that we don't have duplicated data and read the other matrix files
if (primary_input_file %in% input_files) {
    message('You have selected the same input file twice. If you choose a file to be the primary input data, you cannot choose the same file again for the remaining inputs.')
    quit(status=1)
}
mtx_list <- lapply(input_files, read.table, sep='\t', row.names=1, header=T)

# check that the other matrices have the same number of columns (using the 'primary'
# matrix as the one which determines the correct dimension)
n <- length(sample_ids)
df_sizes = lapply(mtx_list, function(x) dim(x)[2])
all_same_size <- all(unlist(lapply(df_sizes, function(x) x == n)))
if (!all_same_size){
    message(sprintf('Not all of the input matrices had the same number of columns. We were expecting all inputs to have %d columns, but we found the following numbers of columns: %s. Please check your inputs.', n, paste(df_sizes, collapse=', ')))
    quit(status=1)
}

# we can't check by name, so we have to check the contents of each matrix to ensure the user did not duplicate
# the data. Since WebMeV assigns a unique path/UUID to each file input, we can't check for duplicated data by file path.
check_for_duplicated_data <- function(df, primary_df) {
    out <- tryCatch(
        {
            is_equal = df == primary_df
            if(all(is_equal)){
                return(TRUE)
            }
            return(FALSE)
        },
        error=function(x){
            return(FALSE)
        }
    )
    return(out)
}
duplicated <- lapply(mtx_list, check_for_duplicated_data, primary_mtx)
if (any(unlist(duplicated))){
    message('You have selected the same input file twice. If you choose a file to be the primary input data, you cannot choose the same file again for the remaining inputs.')
    quit(status=1)
}


# append the primary_mtx dataframe to that list:
mtx_list <- c(list(primary_mtx), mtx_list)

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

# Add some row/col names such that WebMeV can ingest the output. 
# In general, the columns of each matrix could be named
# by an identifier OTHER than the patient ID (e.g. an aliquot ID like in TCGA). In 
# such a case, there is no obvious way to name the output sample-by-sample
# similarity matrix except to have users define the 'primary' matrix as we have above.
# That matrix will set the row/colnames
rownames(W) = sample_ids
colnames(W) = sample_ids

# perform clustering on the fused network.
clustering = data.frame(cluster=spectralClustering(W, num_clusters), row.names=sample_ids)

similarity_mtx_output_file = paste(working_dir,'snf_similarities.tsv', sep='/')
clustering_output_file = paste(working_dir,'snf_clusters.tsv', sep='/')
write.table(W, similarity_mtx_output_file, sep='\t', quote=F)
write.table(clustering, clustering_output_file, sep='\t', quote=F)

# create the expected outputs file:
json_str = paste0(
       '{"snf_similarity":"', similarity_mtx_output_file, '",',
       '"snf_clustering":"', clustering_output_file, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)