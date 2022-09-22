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


#### Single cell 

test_that("Test hierarchicalClones light chain split works", {
    
    db <- data.frame(
        sequence_id=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8"),
        cell_id=c("cell1","cell1","cell2","cell2","cell3","cell3","cell4","cell4"),
        v_call=c("IGHV1*01","IGLV1*01","IGHV1*01","IGLV1*01","IGHV1*01","IGLV2*01", "IGHV3*01,IGHV1*01","IGLV1*01"),
        d_call=c("IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA),
        j_call=c("IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01*01","IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01"),
        junction=c("TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTTTTC","TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTC")
    )
    db$chain <- "light"
    db$chain[grepl("IGH",db[['v_call']])] <- "heavy"
    db$locus <- alakazam::getLocus(db$v_call)
    db$junction_len <- stringi::stri_length(db[['junction']])  

    # prepare_db uses groupGenes to find vj gene groups, but 
    # uses junc_len = NULL
    # It adds the junction length groups outside groupGenes
    # requires both cell_id and locus
    # scoper:::prepare_db
    db_only_heavy_T_locus <- scoper:::prepare_db(db,only_heavy = T, 
                                                      cell_id="cell_id",
                                                      locus="locus", first=F)$db
    expect_true(all(db_only_heavy_T_locus$vj_group =="G1"))
    
    # scoper:::prepare_db
    db_only_heavy_F_locus_first_T <- scoper:::prepare_db(db,only_heavy = F, 
                                           cell_id="cell_id",
                                           locus="locus",
                                           first=T)$db
    expect_equal(db_only_heavy_F_locus_first_T$vj_group, c("G1","G1","G1","G1","G2","G2","G3","G3"))
    
    db_only_heavy_F_locus_first_F <- scoper:::prepare_db(db,only_heavy = F, 
                                                       cell_id="cell_id",
                                                       locus="locus",
                                                       first=F)$db
    expect_equal(db_only_heavy_F_locus_first_F$vj_group, c("G1","G1","G1","G1","G2","G2","G1","G1"))
    
    
    ## Test hierachicalClones
    ## case only_heavy   split_light first
    ## 1    T            F           F
    ## 2    T            F           T
    ## 3    T            T           F
    ## 4    T            T           T
    ## 5    F            F           F
    ## 6    F            F           T
    ## 7    F            T           F
    ## 8    F            T           T
    
    # Case 1
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=TRUE,
        split_light=FALSE,
        first=F, # default is first=F
        nproc=1)
    expect_true(all(clones@db[['clone_id']] == "1"))
    
    # Case 2
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=TRUE,
        split_light=FALSE,
        first=T, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']],c("1","1","1","1","1","1","2","2"))
    
    # Case 3
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=TRUE,
        split_light=TRUE,
        first=F, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']], c("1","1","2","2","2","2","3","3"))
    
    # Case 4
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=TRUE,
        split_light=TRUE,
        first=T, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']], c("1","1","2","2","3","3","4","4"))
    
    # Case 5
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=FALSE,
        split_light=FALSE,
        first=F, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']],c("1","1","1","1","1","1","2","2"))    
    
    # Case 6
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=FALSE,
        split_light=FALSE,
        first=T, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']],c("1","1","1","1","2","2","3","3"))    
    
    # Case 7
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=FALSE,
        split_light=TRUE,
        first=F, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']],c("1","1","2","2","2","2","3","3"))    
    
    # Case 8
    clones <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=FALSE,
        split_light=TRUE,
        first=T, # default is first=F
        nproc=1)
    expect_equal(clones@db[['clone_id']],c("1","1","2","2","3","3","4","4"))      
    
    ## Wwhat happens with the second groupGenes (inside the light chain split) when
    # first=false, but the "linker" ambiguous call was left out of the same cluster id
    # because of the distance threshold? This groupGenes could be splitting again by vj calls...
    # seq2 is the linker, and the juction is one nt different. Everything else is the same.
    db <- data.frame(
        sequence_id=c("seq1","seq2","seq3","seq4","seq5","seq6"),
        cell_id=c("1","2","3","1","2","3"),
        v_call=c("IGHV1", "IGHV1,IGHV2","IGHV2","IGLV1","IGLV1","IGLV1"),
        j_call=c("IGHJ1","IGHJ1","IGHJ1","IGHJ1","IGHJ1","IGHJ1"),
        junction=c("TCGAAATTC","TCGAACTTC","TCGAAATTC","TCGAAATTC","TCGAAATTC","TCGAAATTC")
    )
    db$locus <- alakazam::getLocus(db$v_call)
    db
    clones_split_F <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=TRUE,
        split_light=FALSE,
        first=F, # default is first=F
        nproc=1)
    
    clones_split_T <- hierarchicalClones(
        db,
        threshold=0,
        cell_id='cell_id',
        locus='locus',
        only_heavy=TRUE,
        split_light=TRUE,
        first=F, # default is first=F
        nproc=1)
    # expecting same results because the light chains are the same.
    expect_equal(clones_split_F@db[['clone_id']], clones_split_T@db[['clone_id']])
    
    ## TODO: multiple light chains?
    
    ## TODO: summaries etc are using clone ids created before the light chain split?
})
