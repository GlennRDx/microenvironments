# Paths and directories
root_dir <- "/Users/glero527/Documents/microenvironments"
data_dir <- file.path(root_dir, "01_data")
figs_dir <- file.path(root_dir, "03_figures")

# General parameters
classification_method <- "surface_class_nn" # Options: "surface_class_nn", "surface_class_clusters", "surface_class_radius"

# 0.1.1
wts_params <- list(
  n_subsamples = 1000, # Number of cells to keep per cell type (to reduce computational load)
  n_top_genes = 4000,   # Number of highly deviant genes for selection
  dims    = 1:15     # Number of dimensions to use for dimensionality reduction techniques (PCA, UMAP, t-SNE). Number is based on elbow plot.
)

# 0.2.1
iss_params <- list(
  radius      = 25,    # Radius for neighbourhood definition
  knn         = 20,    # Number of neighbours for KNN method
  centre_x    = 3249.27,  # Center X coordinate for zooming/plotting
  centre_y    = 2027.26,  # Center Y coordinate for zooming/plotting
  threshold   = 0.9  # Threshold for Core vs Surface classification
)

# 0.3.1
lda_params <- list(
  seed        = 123,
  train_prop  = 2/3,   # Proportion of shared genes for training
  split_prop  = 0.9    # Train/Test split ratio for cells
)

# 0.3.2
# svr_params <- list(
#   seed        = 123,
#   train_split = 0.7,
#   kernel      = "linear",
#   type        = "eps-regression"
# )

# 0.3.3
# rf_params <- list(
#   seed        = 42,
#   n_x_bins    = 5,    # Spatial bins for cross-validation
#   n_y_bins    = 5,
#   train_split = 0.8,
#   ntree       = 500,
#   nodesize    = 5
# )

# 0.4.1
prediction_params <- list(
  prob_threshold = 0.85 # Probability threshold for surface/core predictions
)

# 0.5.1
validation_params <- list(
  n_pseudo_reps = 2,      # Number of pseudo-replicates
  marker_method = "negbinom"  # Options: "MAST", "wilcox", "DESeq2", "negbinom"
)
