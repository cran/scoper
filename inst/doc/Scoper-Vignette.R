## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Load scoper
library("scoper")

## ---- eval=FALSE, warning=FALSE, message=FALSE---------------------------
#  # Load scoper
#  library("scoper")
#  defineClonesScoper(db,
#                     model = c("identical", "hierarchical", "spectral"),
#                     method = c("nt", "aa", "single", "average", "complete",
#                                "novj", "vj"),
#                     germline_col = "GERMLINE_IMGT", sequence_col = "SEQUENCE_IMGT",
#                     junction_col = "JUNCTION",
#                     v_call_col = "V_CALL", j_call_col = "J_CALL",
#                     clone_col = c("CLONE", "clone_id"),
#                     targeting_model = NULL, len_limit = NULL, first = FALSE,
#                     cdr3 = FALSE, mod3 = FALSE, max_n = NULL, threshold = NULL,
#                     base_sim = 0.95, iter_max = 1000, nstart = 1000, nproc = 1,
#                     verbose = FALSE, log_verbose = FALSE, out_dir = ".",
#                     summerize_clones = FALSE)

## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Clonal assignment using hierarchical model
results <- defineClonesScoper(db = ExampleDb, clone_col = "CLONE",
                              model = "hierarchical", method = "single", 
                              threshold = 0.15, summerize_clones = TRUE)
# cloned data (a data.frame)
cloned_db <- results$db
# print effective threshold (a numeric)
results$eff_threshold
# get inter and intra conal distances (a data.frame)
df <- results$inter_intra
# histogram of inter versus intra clonal distances  (a ggplot).
results$plot_inter_intra


# Clonal assignment using spectral model
# IMGT_V object from shazam package to identify sequence limit length
library("shazam")
results <- defineClonesScoper(db = ExampleDb, clone_col = "clone_id",
                              model = "spectral", method = "vj", 
                              len_limit = shazam::IMGT_V,
                              targeting_model = shazam::HH_S5F,
                              sequence_col = "SEQUENCE_IMGT",
                              germline_col = "GERMLINE_IMGT_D_MASK",
                              threshold = 0.15, 
                              summerize_clones = TRUE)
# cloned data (a data.frame)
cloned_db <- results$db
# print effective threshold (a numeric)
results$eff_threshold
# get inter and intra conal distances (a data.frame)
df <- results$inter_intra
# histogram of inter versus intra clonal distances  (a ggplot).
results$plot_inter_intra


