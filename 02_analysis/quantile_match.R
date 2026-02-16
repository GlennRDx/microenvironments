load(file='KI_ISS_scRNA_seq_new_1.RData')

#------------------------------------------------------------------------------#

# Subset both datasets to include only "Endothelial" cells for initial testing
test_ISS = ISS_data_1[,names(cell_type_dict_ISS)[cell_type_dict_ISS == "Endothelial"]]
test_RNA = scData_3p_logRPM[,names(cell_type_dict_scRNA)[cell_type_dict_scRNA == "Endothelial"]]

# Exploratory Plot: Compare quantiles of ABCC11 expression between ISS and RNA-seq
plot(quantile(test_ISS['ABCC11',], seq(0,1,0.001)), quantile(test_RNA['ABCC11',], seq(0,1,0.001)))

# Visualise the cumulative distribution functions (CDFs) for ABCC11
plot(seq(0,1,0.001), quantile(test_ISS['ABCC11',], seq(0,1,0.001)), type='o')
points(seq(0,1,0.001), quantile(test_RNA['ABCC11',], seq(0,1,0.001)), type='o')

# Create a quantile dictionary from ISS data based on the number of cells in the RNA-seq data
quant_dict = quantile(test_ISS['ABCC11',], seq(0,1,1/(dim(test_RNA)[2]-1)))

# Map RNA-seq expression values to ISS quantiles based on their rank
test_match = quant_dict[rank(test_RNA['ABCC11',], ties.method='first')]

# Verify the match: The Q-Q plot should now be a straight diagonal line
plot(quantile(test_ISS['ABCC11',], seq(0,1,0.001)), quantile(test_match, seq(0,1,0.001)))

# Check the number of overlapping genes between datasets
length(intersect(rownames(test_ISS), rownames(test_RNA)))

# Function to automate quantile matching for a specific gene across datasets
match_quantiles_per_genes <- function(ISS_data, scRNAseq_data, gene)
{
  # Generate quantiles from ISS reference
  quant_dict = quantile(ISS_data[gene,], seq(0,1,1/(dim(scRNAseq_data)[2]-1)))
  # Reassign RNA-seq values based on their rank using the ISS distribution
  return(quant_dict[rank(scRNAseq_data[gene,], ties.method='first')])
}

# Test the function on a single gene
match_quantiles_per_genes(test_ISS, test_RNA, 'ABCC11')

# Apply the matching function to all common genes for the current cell type
result = sapply(intersect(rownames(test_ISS), rownames(test_RNA)), match_quantiles_per_genes, ISS_data=test_ISS, scRNAseq_data=test_RNA)
rownames(result) = colnames(test_RNA)
result = t(result)

# Validation plot for another gene (ACTA2) after transformation
plot(quantile(test_ISS['ACTA2',], seq(0,1,0.001)), quantile(result['ACTA2',], seq(0,1,0.001)))

# Perform PCA to visualise how well the distributions align globally
PCA = prcomp(t(cbind(log2(test_ISS[rownames(result),]+1), log2(result+1))))
plot(PCA$x, col=c(rep(1,dim(test_ISS)[2]),rep(2,dim(result)[2])), pch=16, cex=0.5)

# seems to work better

# now run this for all cell types and bind them together

# Initialise an empty dataframe to store combined results
all_results = data.frame(row.names=intersect(rownames(test_ISS), rownames(test_RNA)))

# Loop through every unique cell type defined in the scRNA-seq dictionary
for (cell_type in unique(cell_type_dict_scRNA))
{
  # Subset ISS and RNA-seq data for the specific cell type
  test_ISS = ISS_data_1[,names(cell_type_dict_ISS)[cell_type_dict_ISS == cell_type]]
  test_RNA = scData_3p_logRPM[,names(cell_type_dict_scRNA)[cell_type_dict_scRNA == cell_type]]
  
  # Perform quantile matching for all overlapping genes
  result = sapply(intersect(rownames(test_ISS), rownames(test_RNA)), match_quantiles_per_genes, ISS_data=test_ISS, scRNAseq_data=test_RNA)
  
  # Clean up row names and merge into the master dataframe
  rownames(result) = colnames(test_RNA)
  all_results = cbind(all_results, t(result))
}

# Check final dimensions and preview the head of the matched dataset
dim(all_results)
top(all_results)

# Export the final matched dataset to a CSV file
write.csv(all_results, file='scRNAseq_data_matched.txt')