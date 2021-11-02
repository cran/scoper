# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "ExampleDb.rda"), envir=e1)
#load(file.path("tests", "data-tests", "Example10x.rda"), envir=e1)
load(file.path("..", "data-tests", "ExampleDb.rda"), envir=e1)
load(file.path("..", "data-tests", "Example10x.rda"), envir=e1)
db <- get("ExampleDb", envir=e1)
db_sc <- get("Example10x", envir=e1)
rm(e1)

#ensure older version of sample() used
R_v <- paste(version$major, version$minor,sep=".")
if ( numeric_version(R_v) >= numeric_version("3.6.0") ) {
    expect_warning(RNGkind(sample.kind="Round"))
}

# Check for pipelines environment
# pipeline_env <- Sys.getenv("CI") == "true"
# cat("Bitbucket Pipelines:", pipeline_env, "\n")
pipeline_env <- TRUE

#### clone - identicalClones ####

test_that("Test identicalClones", {
    # Truth
    expects <- as.integer(c(20, 21, 23, 26, 27, 28, 30, 44, 50, 100))
    
    # Reproduce example
    db <- identicalClones(ExampleDb, method ="nt", 
                          junction = "junction", v_call = "v_call", 
                          j_call = "j_call", summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        db <- identicalClones(ExampleDb, method ="nt",
                              junction = "junction", v_call = "v_call",
                              j_call = "j_call", summarize_clones = FALSE,
                              nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

#### clone - hierarchicalClones ####

test_that("Test hierarchicalClones", {
    # Truth
    expects <- as.integer(c(7, 8, 8, 8, 8, 9, 10, 11, 12, 683))
    
    # Reproduce example
    db <- hierarchicalClones(ExampleDb, threshold = 0.15,
                             method = "nt", linkage = "single",
                             junction = "junction", 
                             v_call = "v_call", j_call = "j_call", 
                             summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        db <- hierarchicalClones(ExampleDb, threshold = 0.15,
                                 method = "nt", linkage = "single",
                                 junction = "junction",
                                 v_call = "v_call", j_call = "j_call",
                                 summarize_clones = FALSE,
                                 nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

#### clone - spectralClones - novj method ####

test_that("Test spectralClones - novj", {
    # Truth
    expects <- c(6, 7, 7, 7, 7, 8, 10, 11, 12, 679)
    
    # Reproduce example
    set.seed(12345)
    db <- spectralClones(ExampleDb, method = "novj", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.20,
                         summarize_clones = FALSE)
    clones  <- as.numeric(tail(sort(table(db$clone_id)), 10))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        set.seed(12345)
        db <- spectralClones(ExampleDb, method = "novj",
                             junction = "junction", v_call = "v_call",
                             j_call = "j_call", threshold=0.20,
                             summarize_clones = FALSE,
                             nproc=2)
        clones  <- as.numeric(tail(sort(table(db$clone_id)), 10))
        expect_identical(clones, expects)
    }
})

#### clone - spectralClones - vj method ####

test_that("Test spectralClones - vj", {
    # Truth
    expects <- as.integer(c(11, 12, 12, 13, 14, 15, 16, 29, 35, 683))
    
    # Reproduce example
    set.seed(12345)
    db <- spectralClones(ExampleDb, method = "vj", 
                         germline = "germline_alignment_d_mask",
                         sequence = "sequence_alignment", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.15,
                         summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        set.seed(12345)
        db <- spectralClones(ExampleDb, method = "vj",
                             germline = "germline_alignment_d_mask",
                             sequence = "sequence_alignment",
                             junction = "junction", v_call = "v_call",
                             j_call = "j_call", threshold=0.15,
                             summarize_clones = FALSE,
                             nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

#### Single cell 

test_that("Test assigning clones works for heavy-only sc data", {
    # Issue https://bitbucket.org/kleinstein/scoper/issues/23
    db_sc$chain <- "light"
    db_sc$chain[grepl("IGH",db_sc[['v_call']])] <- "heavy"
    db_sc_heavy <- db_sc %>%
        filter(chain == "heavy")
    expect_warning(cloned <- identicalClones(db_sc_heavy, method="aa",
                               cell_id = "cell_id",
                               locus = "locus", nproc=1),
                   "Single cell mode requested, but")
})
