######################### this section contains source code from the Seurat package: https://github.com/satijalab/seurat/blob/13b615c27eeeac85e5c928aa752197ac224339b9/R/clustering.R #############
######################### slight modifications have been made to this particular section of source code better accomdate the purpose of this benchmark project ######################################

CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

RunLeiden <- function(
  object,
  method = c("matrix", "igraph"),
  partition.type = c(
    'RBConfigurationVertexPartition',
    'ModularityVertexPartition',
    'RBERVertexPartition',
    'CPMVertexPartition',
    'MutableVertexPartition',
    'SignificanceVertexPartition',
    'SurpriseVertexPartition'
  ),
  initial.membership = NULL,
  node.sizes = NULL,
  resolution.parameter = 1,
  random.seed = 0,
  n.iter = 10
) {
  if (!py_module_available(module = 'leidenalg')) {
    stop(
      "Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).",
      call. = FALSE
    )
  }
  switch(
    EXPR = method,
    "matrix" = {
      input <- as(object = object, Class = "matrix")
    },
    "igraph" = {
      input <- if (inherits(x = object, what = 'list')) {
        graph_from_adj_list(adjlist = object)
      } else if (inherits(x = object, what = c('dgCMatrix', 'matrix', 'Matrix'))) {
        if (inherits(x = object, what = 'Graph')) {
          object <- as(object = object, Class = "dgCMatrix")
        }
        graph_from_adjacency_matrix(adjmatrix = object, weighted = TRUE)
      } else if (inherits(x = object, what = 'igraph')) {
        object
      } else {
        stop(
          "Method for Leiden not found for class", class(x = object),
          call. = FALSE
        )
      }
    },
    stop("Method for Leiden must be either 'matrix' or igraph'")
  )
  #run leiden from CRAN package (calls python with reticulate)
  partition <- leiden(
    object = input,
    partition_type = partition.type,
    initial_membership = initial.membership,
    weights = NULL,
    node_sizes = node.sizes,
    resolution_parameter = resolution.parameter,
    seed = random.seed,
    n_iterations = n.iter
  )
  return(partition)
}

FindClusters.default.alt <- function(
  object,
  modularity.fxn = 1,
  initial.membership = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  method = "matrix",
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  if (is.null(x = object)) {
    stop("Please provide an SNN graph")
  }
  if (tolower(x = algorithm) == "louvain") {
    algorithm <- 1
  }
  if (tolower(x = algorithm) == "leiden") {
    algorithm <- 4
  }
  if (nbrOfWorkers() > 1) {
    clustering.results <- future_lapply(
      X = resolution,
      FUN = function(r) {
        if (algorithm %in% c(1:3)) {
          ids <- RunModularityClustering(
            SNN = object,
            modularity = modularity.fxn,
            resolution = r,
            algorithm = algorithm,
            n.start = n.start,
            n.iter = n.iter,
            random.seed = random.seed,
            print.output = verbose,
            temp.file.location = temp.file.location,
            edge.file.name = edge.file.name
          )
        } else if (algorithm == 4) {
          ids <- RunLeiden(
            object = object,
            method = method,
            partition.type = "RBConfigurationVertexPartition",
            initial.membership = initial.membership,
            node.sizes = node.sizes,
            resolution.parameter = r,
            random.seed = random.seed,
            n.iter = n.iter
          )
        } else {
          stop("algorithm not recognised, please specify as an integer or string")
        }
        names(x = ids) <- colnames(x = object)
        ids <- GroupSingletons(ids = ids, SNN = object, verbose = verbose)
        results <- list(factor(x = ids))
        names(x = results) <- paste0('res.', r)
        return(results)
      }
    )
    clustering.results <- as.data.frame(x = clustering.results)
  } else {
    clustering.results <- data.frame(row.names = colnames(x = object))
    for (r in resolution) {
      if (algorithm %in% c(1:3)) {
        ids <- RunModularityClustering(
          SNN = object,
          modularity = modularity.fxn,
          resolution = r,
          algorithm = algorithm,
          n.start = n.start,
          n.iter = n.iter,
          random.seed = random.seed,
          print.output = verbose,
          temp.file.location = temp.file.location,
          edge.file.name = edge.file.name)
      } else if (algorithm == 4) {
        ids <- RunLeiden(
          object = object,
          method = method,
          partition.type = "RBConfigurationVertexPartition",
          initial.membership = initial.membership,
          node.sizes = node.sizes,
          resolution.parameter = r,
          random.seed = random.seed,
          n.iter = n.iter
        )
      } else {
        stop("algorithm not recognised, please specify as an integer or string")
      }
      names(x = ids) <- colnames(x = object)
      #ids <- GroupSingletons(ids = ids, SNN = object, group.singletons = group.singletons, verbose = verbose)
      clustering.results[, paste0("res.", r)] <- factor(x = ids)
    }
  }
  return(clustering.results)
}


FindClusters.alt <- function(
  object,
  graph.name = NULL,
  modularity.fxn = 1,
  initial.membership = NULL,
  node.sizes = NULL,
  resolution = 0.8,
  method = "matrix",
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  temp.file.location = NULL,
  edge.file.name = NULL,
  verbose = TRUE,
  ...
) {
  CheckDots(...)
  graph.name <- graph.name %||% paste0(DefaultAssay(object = object), "_snn")
  if (!graph.name %in% names(x = object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if (!is(object = object[[graph.name]], class2 = "Graph")) {
    stop("Provided graph.name does not correspond to a graph object.")
  }
  clustering.results <- FindClusters.default.alt(
    object = object[[graph.name]],
    modularity.fxn = modularity.fxn,
    initial.membership = initial.membership,
    node.sizes = node.sizes,
    resolution = resolution,
    method = method,
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    group.singletons = group.singletons,
    temp.file.location = temp.file.location,
    edge.file.name = edge.file.name,
    verbose = verbose,
    ...
  )
  colnames(x = clustering.results) <- paste0(graph.name, "_", colnames(x = clustering.results))
  object <- AddMetaData(object = object, metadata = clustering.results)
  Idents(object = object) <- colnames(x = clustering.results)[ncol(x = clustering.results)]
  levels <- levels(x = object)
  levels <- tryCatch(
    expr = as.numeric(x = levels),
    warning = function(...) {
      return(levels)
    },
    error = function(...) {
      return(levels)
    }
  )
  Idents(object = object) <- factor(x = Idents(object = object), levels = sort(x = levels))
  object[['seurat_clusters']] <- Idents(object = object)
  cmd <- LogSeuratCommand(object = object, return.command = TRUE)
  slot(object = cmd, name = 'assay.used') <- DefaultAssay(object = object[[graph.name]])
  object[[slot(object = cmd, name = 'name')]] <- cmd
  return(object)
}


#Seuratv4: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
#data1: gene x cell raw
#  (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python

combine_Seuratv4 = function(data1, data2, data1_name, data2_name,#input
                            sf = 6000, center = F, scale = F, #normalization and scale
                            data1_npcs = 100, data2_npcs=100, data1_reduce = T, data2_reduce = T, #max dimension reduction
                            data1_dims, data2_dims, n_neighbors,smooth=F,#integrate data
                            clustering = T, res = 0.8, res_step = 0.4, #clustering 
                            alg=4, ncl_truth, max_iter = 50, seed = 1234, n.iter=1000
                            ){
  library(Seurat)
  library(dplyr)
  library(future)
  library(reticulate)
  library(leiden)
  #create combined Seurat object
  combine_SO = CreateSeuratObject(counts = data1, assay=data1_name)
  data2_assay = CreateAssayObject(counts = data2)
  combine_SO[[data2_name]] = data2_assay
  #Assays(combine_SO)
  
  #preprocess and dimension reduction
  DefaultAssay(combine_SO) = data1_name
  VariableFeatures(combine_SO) = rownames(data1)
  combine_SO = NormalizeData(combine_SO, scale.factor = sf)%>%
    ScaleData(do.scale=scale, do.center=center) %>% 
    RunPCA(npcs = data1_npcs, reduction.name=paste0(data1_name,"_pca")) 
  if(!data1_reduce){
    combine_SO@reductions[[paste0(data1_name,"_pca")]]@cell.embeddings = t(combine_SO@assays[[data1_name]]@scale.data)
    colnames(combine_SO@reductions[[paste0(data1_name,"_pca")]]@cell.embeddings) = paste("PC", 1:nrow(data1), sep="_")
  }

  DefaultAssay(combine_SO) = data2_name
  VariableFeatures(combine_SO) = rownames(data2)
  combine_SO = NormalizeData(combine_SO, scale.factor = sf)%>%
    ScaleData(do.scale = scale, do.center = center) %>% 
    RunPCA(npcs = data2_npcs, reduction.name=paste0(data2_name,"_pca")) 
  if(!data2_reduce){
    combine_SO@reductions[[paste0(data2_name,"_pca")]]@cell.embeddings = t(combine_SO@assays[[data2_name]]@scale.data)
    colnames(combine_SO@reductions[[paste0(data2_name,"_pca")]]@cell.embeddings) = paste("PC", 1:nrow(data2), sep="_")
  }
  
  #integrate data
  #combine_SO = FindMultiModalNeighbors(combine_SO, 
  #                                     reduction.list = list(paste0(data1_name,"_pca"),
  #                                                           paste0(data2_name,"_pca")),
  #                                     dims.list = list(1:data1_dims, 1:data2_dims),
  #                                     k.nn = n_neighbors, smooth = smooth)
  
  #clustering
  if(clustering == T){
    res_init = res
    num_cl = 0
    iter0 = 1
    while(num_cl != ncl_truth & n_neighbors > 0 & iter0 < 15){
      #combine_SO = FindMultiModalNeighbors(combine_SO,
      #                                 reduction.list = list(paste0(data1_name,"_pca"),
      #                                                       paste0(data2_name,"_pca")),
      #                                 dims.list = list(1:data1_dims, 1:data2_dims),
      #                                 k.nn = n_neighbors, smooth = smooth)
       
      if(num_cl != 0){
        if(num_cl > ncl_truth){
          if(n_neighbors <= 10){ 
            n_neighbors = n_neighbors + 1
          }else{
            n_neighbors = n_neighbors + 5
          }
        }else if(num_cl < ncl_truth){
          if(n_neighbors <= 10){
            n_neighbors = n_neighbors - 1
          }else{
            n_neighbors = n_neighbors - 5
          }
        }
      }

      set.seed(seed)
      combine_SO = FindMultiModalNeighbors(combine_SO,
                                       reduction.list = list(paste0(data1_name,"_pca"),
                                                             paste0(data2_name,"_pca")),
                                       dims.list = list(1:data1_dims, 1:data2_dims),
                                       k.nn = n_neighbors, smooth = smooth)
      iter0 = iter0 + 1
      iter = 1
      res = res_init
      while(num_cl != ncl_truth & iter < max_iter & res <= 1 & res > 9.1e-4){
        set.seed(seed)
        reticulate::py_set_seed(seed)
        combine_SO = FindClusters.alt(combine_SO, graph.name = "wsnn", algorithm = alg,
                                   resolution = res, random.seed = seed, n.iter=n.iter)
        num_cl = length(table(combine_SO@meta.data$seurat_clusters))
        if(num_cl > ncl_truth){
          #while(res - res_step/iter < 1e-3){
          #  res_step = res_step/2
          #}
          res = res - res_step/iter
        }else if(num_cl < ncl_truth){
          res = res + res_step/iter
        }
        print(paste("number of clusters eis", num_cl))
        print(paste("resolution is", res))
        print(iter)
        iter = iter + 1
      }
    print(n_neighbors)
    }
  }
  
  return(list(SO = combine_SO, num_nn = n_neighbors))
}
