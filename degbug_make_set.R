test_set = function (enadata, dimensions = 2, norm.by = fun_sphere_norm, 
          rotation.by = ena.svd, rotation.params = NULL, rotation.set = NULL, 
          endpoints.only = TRUE, center.align.to.origin = TRUE, node.position.method = lws.positions.sq, 
          as.list = TRUE, ...) 
{ browser()
  if (as.list == F) {
    warning(paste0("Usage of ENAdata and ENAset objects will be deprecated ", 
                   "and potentially removed altogether in future versions."))
    if (!is(enadata, "ENAdata")) {
      stop(paste0("Use of ena.make.set with as.list=FALSE requires `enadata` ", 
                  "be an ENAdata object. Re-run the accumulation with as.list=FALSE"))
    }
    set <- ENAset$new(enadata = enadata, dimensions = dimensions, 
                      rotation.by = ifelse(!is.null(rotation.by) && identical(rotation.by, 
                                                                              ena.svd), ena.svd.R6, rotation.by), rotation.params = rotation.params, 
                      rotation.set = rotation.set, norm.by = norm.by, node.position.method = ifelse(identical(node.position.method, 
                                                                                                              lws.positions.sq), lws.positions.sq.R6, node.position.method), 
                      endpoints.only = endpoints.only, center.align.to.origin = center.align.to.origin, 
                      ...)
    return(set$process())
  }
  else {
    if ("ENAdata" %in% class(enadata)) {
      warning(paste0("Usage of ENAdata objects will be deprecated and ", 
                     "potentially removed altogether in future versions. See ", 
                     "ena.accumulate.data() or ena.set()."))
      enadata <- ena.set(enadata)
    }
    enadata$`_function.params`$center.align.to.origin <- center.align.to.origin
    enadata$`_function.params`$rotation.by <- rotation.by
    enadata$`_function.params`$rotation.params <- rotation.params
    code_columns <- rENA:::svector_to_ut(enadata$rotation$codes)
    line.weights <- norm.by(as.matrix(enadata$connection.counts))
    colnames(line.weights) <- code_columns
    line.weights.dt <- data.table::as.data.table(line.weights)
    for (i in seq(ncol(line.weights.dt))) rENA:::set(line.weights.dt, 
                                              j = i, value = as.ena.co.occurrence(line.weights.dt[[i]]))
    enadata$line.weights <- cbind(enadata$meta.data, line.weights.dt)
    class(enadata$line.weights) <- c("ena.line.weights", 
                                     "ena.matrix", class(enadata$line.weights))
    if (!is.null(rotation.set)) {
      if (inherits(rotation.set, "ena.rotation.set")) {
        if (center.align.to.origin) {
          points.for.projection <- line.weights
          points.for.projection[rowSums(as.matrix(line.weights)) != 
                                  0, ] <- center.projection(lws = line.weights[rowSums(as.matrix(line.weights)) != 
                                                                                 0, ], rotation = rotation.set)
        }
        else {
          points.for.projection <- center.projection(lws = line.weights, 
                                                     rotation = rotation.set)
        }
      }
      else {
        stop("Supplied rotation.set is not an instance of ENARotationSet")
      }
    }
    else {
      if (center.align.to.origin) {
        points.for.projection <- line.weights
        points.for.projection[rowSums(as.matrix(line.weights)) != 
                                0, ] <- center_data_c(line.weights[rowSums(as.matrix(line.weights)) != 
                                                                     0, ])
      }
      else {
        points.for.projection <- center_data_c(line.weights)
      }
    }
    colnames(points.for.projection) <- code_columns
    enadata$model$points.for.projection = data.table::as.data.table(points.for.projection)
    for (i in seq(ncol(enadata$model$points.for.projection))) {
      set(enadata$model$points.for.projection, j = i, value = as.ena.co.occurrence(enadata$model$points.for.projection[[i]]))
    }
    enadata$model$points.for.projection <- as.ena.matrix(cbind(enadata$meta.data, 
                                                               enadata$model$points.for.projection), "ena.points")
    if (!is.null(rotation.by) && is.null(rotation.set)) {
      rotation <- do.call(rotation.by, list(enadata, rotation.params))
      enadata$rotation.matrix <- as.data.table(rotation$rotation, 
                                               keep.rownames = "codes")
      for (i in seq(ncol(enadata$rotation.matrix))) {
        if (i == 1) {
          set(enadata$rotation.matrix, j = i, value = as.ena.metadata(enadata$rotation.matrix[[i]]))
        }
        else {
          set(enadata$rotation.matrix, j = i, value = as.ena.dimension(enadata$rotation.matrix[[i]]))
        }
      }
      class(enadata$rotation.matrix) <- c("ena.rotation.matrix", 
                                          class(enadata$rotation.matrix))
      enadata$rotation$rotation.matrix <- enadata$rotation.matrix
      enadata$rotation$eigenvalues <- rotation$eigenvalues
      if (center.align.to.origin) {
        enadata$rotation$center.vec = colMeans(line.weights[rowSums(as.matrix(line.weights)) != 
                                                              0, ])
      }
      else {
        enadata$rotation$center.vec = colMeans(line.weights)
      }
    }
    else if (!is.null(rotation.set)) {
      if (is(rotation.set, "ena.rotation.set")) {
        enadata$rotation.matrix <- rotation.set$rotation.matrix
        enadata$rotation$rotation.matrix <- rotation.set$rotation.matrix
        enadata$rotation$nodes <- rotation.set$nodes
        enadata$rotation$eigenvalues <- rotation.set$eigenvalues
        enadata$rotation$center.vec = rotation.set$center.vec
      }
      else {
        stop("Supplied rotation.set is not an instance of ENARotationSet")
      }
    }
    else {
      stop("Unable to find or create a rotation set")
    }
    if (!is.null(enadata$rotation.matrix)) {
      points <- points.for.projection %*% as.matrix(enadata$rotation.matrix)
      points.dt <- as.data.table(points)
      for (i in seq(ncol(points.dt))) {
        set(points.dt, j = i, value = as.ena.dimension(points.dt[[i]]))
      }
      if (grepl(x = enadata$model$model.type, pattern = "Trajectory")) {
        enadata$points <- cbind(enadata$trajectories, 
                                points.dt)
      }
      else {
        enadata$points <- cbind(enadata$meta.data, points.dt)
      }
      enadata$points <- as.ena.matrix(enadata$points, "ena.points")
    }
    else {
      stop(paste0("There is no rotation matrix, if you supplied a custom ", 
                  "rotation.set, be sure it contains a rotation.matrix"))
    }
    if (exists("rotation") && !is.null(rotation) && is.null(rotation.set)) {
      positions <- node.position.method(enadata)
      if (all(names(positions) %in% c("node.positions", 
                                      "centroids"))) {
        enadata$rotation$nodes <- as.data.table(positions$node.positions)
        colnames(enadata$rotation$nodes) <- colnames(points)
        rownames(enadata$rotation$nodes) <- enadata$rotation$codes
        for (i in seq(ncol(enadata$rotation$nodes))) {
          set(enadata$rotation$nodes, j = i, value = as.ena.dimension(enadata$rotation$nodes[[i]]))
        }
        enadata$rotation$nodes <- data.table(code = structure(enadata$rotation$codes, 
                                                              class = c("code", class(enadata$rotation$codes))), 
                                             enadata$rotation$nodes)
        class(enadata$rotation$nodes) = c("ena.nodes", 
                                          class(enadata$rotation$nodes))
        enadata$model$centroids <- as.data.table(positions$centroids)
        for (i in seq(ncol(enadata$model$centroids))) {
          set(enadata$model$centroids, j = i, value = as.ena.dimension(enadata$model$centroids[[i]]))
        }
        colnames(enadata$model$centroids) <- colnames(as.matrix(enadata$rotation.matrix))
        enadata$model$centroids = cbind(data.table(unit = enadata$model$unit.labels), 
                                        enadata$model$centroids)
        set(enadata$model$centroids, j = 1L, value = as.ena.metadata(enadata$model$centroids[[1L]]))
        enadata$model$centroids <- as.ena.matrix(enadata$model$centroids)
      }
      else {
        stop(paste0("The node position method didn't return back the ", 
                    "expected objects:\n", "\tExpected: c('node.positions','centroids')\n", 
                    "\tReceived: ", names(positions), sep = ""))
      }
    }
    else if (!is.null(rotation.set)) {
      enadata$rotation$nodes <- rotation.set$nodes
    }
    if (is.null(enadata$rotation$nodes)) {
      stop("Unable to determine the node positions either by calculating\n                    them using `node.position.method` or using a supplied\n                    `rotation.set`")
    }
    var_rot_data <- var(points)
    diagonal_variance <- as.vector(diag(var_rot_data))
    enadata$model$variance <- diagonal_variance/sum(diagonal_variance)
    names(enadata$model$variance) <- colnames(enadata$rotation$rotation.matrix)[-1]
    enadata$plots <- list()
    enadata$`_function.params`$norm.by <- norm.by
    return(enadata)
  }
}
