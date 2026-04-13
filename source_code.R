# source_code.R

# This file contains all custom R functions used for the 
# paper "Experimental Designs for Electric Bus Routes."

library(igraph)

# get_all_paths_with_L_edges
# Returns a list of all paths in G that contain exactly L edges.
# G: an undirected graph with m vertices
# L: number of edges in the path
get_all_paths_with_L_edges <- function(G,L){
  
  # Initialize an empty list to store paths of length 5
  paths <- list()
  
  # Get all vertex IDs
  vertex_ids <- V(G)
  
  # Iterate through all possible source vertices
  for (start_node in vertex_ids) {
    # Iterate through all possible target vertices (excluding the start node)
    for (end_node in vertex_ids[vertex_ids != start_node]) {
      # Find all simple paths between start_node and end_node with a cutoff of L
      current_paths <- all_simple_paths(G, from = start_node, to = end_node, cutoff = L)
      
      # Filter for paths that have exactly L edges (length of the path is number of vertices - 1)
      for (path in current_paths) {
        if (length(path) - 1 == L) {
          vertex_names <- V(G)[path]$name
          
          # Convert vertex sequence to edge pairs
          edge_pairs <- cbind(head(vertex_names, -1), tail(vertex_names, -1))
          
          
          # Create custom edge names with "_"
          edge_names <- apply(edge_pairs, 1, function(x) paste(sort(x), collapse = "_"))
          
          paths <- c(paths, list(edge_names))
        }
      }
    }
  }
  
  
  
  return(paths)
}

# This function generates a random connected graph with m vertices and n edges.
generate_random_connected_graph <- function(m, n) {
  # Ensure valid parameters for a connected graph
  if (n < (m - 1)) {
    stop("A connected graph with 'm' vertices requires at least 'm-1' edges.")
  }
  
  repeat {
    # Generate a random graph with m vertices and n edges
    g <- sample_gnm(n = m, m = n, directed = FALSE, loops = FALSE)
    
    # Check if the graph is connected
    if (is_connected(g)) {
      return(g)
    }
  }
}

# This is a helper function used to generate a 2^p full factorial design.
generate_2p_factorial_design <- function(p) {
  
  
  levels <- c(-1, 1) 
  factors_list <- rep(list(levels), p)
  names(factors_list) <- paste0("X", 1:p)
  design <- expand.grid(factors_list)
  
  return(design)
}


# This is a helper function that is called by SA_design. It is used to generate random neighbors of current design points.
get_random_neighbor <- function(x_current, p, X_paths){
  num_edges = length(x_current)-p
  x_edges = x_current[1:num_edges]
  
  new_path = X_paths[sample(1:nrow(X_paths),1),]
  x_neighbor = c(new_path, sample(x = c(-1, 1), size = p, replace = TRUE) )
  
  return(x_neighbor)
  
}

# get_bayes_D
# Returns the Bayesian D-optimality criterion for a design with model matrix X and matrix R (the inverse of the prior variance-covariance matrix)
get_bayes_D <- function(X, R){
  return( (det( t(X)%*%X/nrow(X) + R))^{1/ncol(X)}   )
}


# SA_design: This function implements Algorithm 1
# X_init: an n times |E|+p initial design matrix. Each set of L nonzero edges in X_init must form a valid path in G.
# p: number of factors
# L: an integer in {1,2,...,m} which is the maximum number of stops a bus can visit in a single trip
# R: a diagonal (|E| + p) times (|E| + p) prior covariance matrix
# X_paths: a matrix corresponding to all paths of length L in G (coded using x_{jk}(a) = 0 or 1)
# n_iter: number of iterations (default is 10000)
# This function uses fast exchange (rank 2) updates.
SA_design <- function(X_init, p, L, R, X_paths, 
                      n_iter = 10000){
  
  q = ncol(X_init)
  num_edges = q - p
  if(any(rowSums(X_init[,1:num_edges])  > L)){
    stop("At least one row has more than L edges.")
  }
  
  
  
  X_current = X_init
  A_current = t(X_current)%*%X_current/n + R
  det_current = det(A_current)
  phi_current = det_current^(1/q)
  A_current_inv = solve(A_current)
  
  
  # phi_current = get_bayes_D(X=X_current, R=R)
  
  phi_best = phi_current
  X_best = X_current
  
  
  
  # phi_vec = numeric(n_iter)
  for(t in 1:n_iter){
    
    # print(paste("Beginning Iteration", t, sep = " "))
    
    # randomly select row of X
    rand_row_index = sample(1:n, 1)
    rand_row = X_current[rand_row_index,]
    
    # generate a random neighbor of x_current
    x_neighbor = get_random_neighbor(x_current = rand_row, p = p,
                                     X_paths = X_paths)    
    
    # update X_next
    X_next = X_current
    X_next[rand_row_index,] = x_neighbor
    
    # update phi
    # phi_next = get_bayes_D(X=X_next, R=R)
    fx = matrix(rand_row,ncol = 1)/sqrt(n)
    fy = matrix(x_neighbor,ncol = 1)/sqrt(n)
    vx = t(fx)%*%A_current_inv%*%fx
    vy = t(fy)%*%A_current_inv%*%fy
    vxy = t(fx)%*%A_current_inv%*%fy
    det_next = det_current*(1+vy - vx + vxy^2 -vy*vx)
    phi_next = det_next^(1/q)
    
    # decide whether or not to accept
    U = runif(1, min = 0, max = 1)
    accept_prob = exp(-(phi_current - phi_next)/(1/(t+1)))
    # print(paste("Acceptance Probability: ", accept_prob, sep = ""))
    if(U <= accept_prob){
      # perform the exchange
      X_current = X_next
      F1 = cbind(fy,-fx)
      F2 = cbind(fy, fx)
      tmp1 = solve( diag(2) + t(F2)%*%A_current_inv%*%F1)
      A_current_inv = A_current_inv - A_current_inv%*%F1%*%tmp1%*%t(F2)%*%A_current_inv
      det_current = det_next
      phi_current = phi_next
    }
    
    if(phi_current >= phi_best){
      
      phi_best = phi_current
      X_best = X_current
      
    }
    
    # else do nothing
    # phi_vec[t] = phi_current
    
  }
  
  # plot(1:n_iter, phi_vec, type = "l", xlab = "Iteration", ylab = "Bayesian D-efficiency")
  # abline(h = phi_best, col = "red")
  
  
  
  return(list(model.matrix = X_best, phi = phi_best))
  
}

# This function implements the Fedorov exchange algorithm with the Bayesian D-optimality criterion.
# It uses fast exchange (rank 2) updates.
fedorov_bayes <- function(X_init, X_all, R, n_repeat = 5){
  
  X_best = X_init
  q = ncol(X_init)
  # phi_best = get_bayes_D(X=X_best, R=R)
  
  n = nrow(X_best)
  A_best = t(X_best)%*%X_best/n + R
  det_best = det(A_best)
  phi_best = det_best^(1/q)
  A_best_inv = solve(A_best)
  
  n = nrow(X_init)
  N = nrow(X_all)
  
  for(rep in 1:n_repeat){
    
    
    
    for(i in 1:n){
      
      # X_temp = X_best
      phi_temps = numeric(N)
      det_temps = numeric(N)
      for(j in 1:N){
        # X_temp[i,] = X_all[j,]
        # phi_temps[j] = get_bayes_D(X=X_temp, R=R) 
        # Using exchange formula to speed this up
        fx = matrix(X_best[i,],ncol = 1)/sqrt(n)
        fy = matrix(X_all[j,],ncol = 1)/sqrt(n)
        vx = t(fx)%*%A_best_inv%*%fx
        vy = t(fy)%*%A_best_inv%*%fy
        vxy = t(fx)%*%A_best_inv%*%fy
        det_temps[j] = det_best*(1+vy - vx + vxy^2 -vy*vx)
        phi_temps[j] = (det_temps[j])^(1/q)
        if(is.na(phi_temps[j])){ phi_temps[j] = 0}
      }
      
      max_ind = which.max(phi_temps)
      if(phi_temps[max_ind] > phi_best){
        
        # update A_best_inv
        fx = matrix(X_best[i,],ncol = 1)/sqrt(n)
        fy = matrix(X_all[max_ind,],ncol = 1)/sqrt(n)
        F1 = cbind(fy,-fx)
        F2 = cbind(fy, fx)
        tmp1 = solve( diag(2) + t(F2)%*%A_best_inv%*%F1)
        A_best_inv = A_best_inv - A_best_inv%*%F1%*%tmp1%*%t(F2)%*%A_best_inv
        
        X_best[i,] = X_all[max_ind,]
        det_best = det_temps[max_ind]
        phi_best = phi_temps[max_ind]
        
      }
      
    }
    
  }
  
  
  return(list(model.matrix = X_best, phi = phi_best))
  
  
}

CRPS_analytic_normal <- function(Y, means, variances){
  
  sds = sqrt(variances)
  Zs = (Y-means)/sds
  return(  sds*Zs*(2*pnorm(Zs) - 1) + 2*dnorm(Zs) - 1/sqrt(pi))
  
}


fit_model <- function(X, Y, mu, R){
  XtX_inv = solve(t(X)%*%X + R) 
  theta_hat = XtX_inv%*%(t(X)%*%Y + t(R)%*%mu)
  posterior_variances = diag(XtX_inv)
  # We can find credible intervals for each effect (very similar to confidence interval)
  z_crit = qnorm(0.025, lower.tail = FALSE)
  lower_credible_interval = theta_hat - z_crit*sqrt(posterior_variances)
  upper_credible_interval = theta_hat + z_crit*sqrt(posterior_variances)
  credible_intervals = cbind(lower_credible_interval, upper_credible_interval)
  
  summary_table = round(cbind(theta_hat, posterior_variances, credible_intervals),3)
  colnames(summary_table) = c("Posterior Mean", "Posterior Variance", "Lower 2.5%","Upper 97.5%")
  # print(summary_table)
  # need to estimate sigma2
  yhat = X%*%theta_hat
  sigma2hat = mean((yhat-Y)^2)
  alpha_star_update_factor = 0.5*nrow(X)
  beta_star_update_factor = 0.5*(sum(Y^2) + t(mu)%*%solve(R)%*%mu + t(theta_hat)%*%XtX_inv%*%theta_hat)
  return(list(mean = theta_hat, var = posterior_variances, sigma2hat = sigma2hat, var_cov = XtX_inv,
              summary = summary_table, 
              alpha_star_update_factor = alpha_star_update_factor,
              beta_star_update_factor = as.numeric(beta_star_update_factor)))
}

# obtains model predictions using the "plug-in" estimator for the variance sigma^2.
predict_model <- function(model, X){
  theta_hat = model$mean
  posterior_variances = model$var
  posterior_var_cov = model$var_cov
  sigma2hat = model$sigma2hat
  yhat = X%*%theta_hat
  # pred_var_cov = X%*%posterior_var_cov%*%t(X) + sigma2hat*diag(nrow(X))
  pred_var_cov = sigma2hat*(X%*%posterior_var_cov%*%t(X) + diag(nrow(X)))
  pred_var = diag(pred_var_cov)
  pred_SE = sqrt(pred_var)
  pred_uppers = yhat + 1.96*pred_SE
  pred_lowers = yhat - 1.96*pred_SE
  pred_mat = cbind(yhat, pred_var, pred_lowers, pred_uppers)
  colnames(pred_mat) = c("Predicted Mean", "Variance of Prediction", "95% Lower", "95% Upper")
  return(pred_mat)
}

# obtains model predictions using an IG(alpha0, beta0) prior on sigma^2
predict_model_t <- function(model, X, alpha0, beta0, mu, R){
  alpha_star = alpha0 + model$alpha_star_update_factor
  beta_star = beta0 + model$beta_star_update_factor
  theta_hat = model$mean
  posterior_variances = model$var
  posterior_var_cov = model$var_cov
  yhat = X%*%theta_hat
  # pred_var_cov = X%*%posterior_var_cov%*%t(X) + (beta_star/alpha_star)*diag(nrow(X))
  pred_var_cov = (beta_star/alpha_star)*(X%*%posterior_var_cov%*%t(X) + diag(nrow(X)))
  pred_var = diag(pred_var_cov)
  pred_SE = sqrt(pred_var)
  tcrit = qt(0.025, df = 2*alpha_star, lower.tail = F)
  pred_uppers = yhat + tcrit*pred_SE
  pred_lowers = yhat - tcrit*pred_SE
  pred_mat = cbind(yhat, pred_var, pred_lowers, pred_uppers)
  colnames(pred_mat) = c("Predicted Mean", "Variance of Prediction", "95% Lower", "95% Upper")
  return(pred_mat)
}


# This function is used to help quickly produce the design tables in Appendix C.
# edges: character vector like c("u_v", "v_w", ...)
get_path_from_strings <- function(edges, sep = "_") {
  stopifnot(is.character(edges), length(edges) > 0)
  
  # Parse endpoints for each edge
  pairs <- strsplit(edges, sep, fixed = TRUE)
  if (any(vapply(pairs, length, 1L) != 2)) {
    stop("Each edge string must contain exactly one separator and define two vertices.")
  }
  L <- vapply(pairs, `[`, character(1), 1)
  R <- vapply(pairs, `[`, character(1), 2)
  
  # Build degree counts
  deg <- table(c(L, R))
  vertices <- names(deg)
  # Endpoints are degree-1 vertices
  endpoints <- names(deg[deg == 1])
  if (length(endpoints) != 2) {
    stop("Edges do not form a single path (need exactly two endpoints).")
  }
  
  # Build adjacency: for each vertex, list incident edge indices
  # and the opposite vertex for that edge
  # Use integer indices for speed
  nE <- length(edges)
  # map vertex name -> list of (edge index, other vertex)
  adj <- setNames(vector("list", length(vertices)), vertices)
  for (i in seq_len(nE)) {
    u <- L[i]; v <- R[i]
    adj[[u]] <- c(adj[[u]], list(list(e = i, to = v)))
    adj[[v]] <- c(adj[[v]], list(list(e = i, to = u)))
  }
  
  used <- rep(FALSE, nE)
  v_order <- character(0)
  e_order <- integer(0)
  
  cur <- endpoints[1]
  v_order <- c(v_order, cur)
  # Walk the path
  for (step in seq_len(nE)) {
    # find the only unused edge from current vertex
    nexts <- adj[[cur]]
    # among incident edges, pick the one not yet used
    idx <- which(!used[vapply(nexts, function(x) x$e, integer(1))])
    if (length(idx) != 1) {
      stop("Traversal failed: expected exactly one next unused edge from current vertex.")
    }
    choice <- nexts[[idx]]
    e <- choice$e
    nxt <- choice$to
    used[e] <- TRUE
    e_order <- c(e_order, e)
    v_order <- c(v_order, nxt)
    cur <- nxt
  }
  
  list(
    vertex_order = v_order,
    edge_order = e_order,
    # oriented edge labels in the walked order (left-right as you traverse)
    edge_labels_in_order = sprintf("%s%s%s",
                                   v_order[-length(v_order)], sep, v_order[-1])
  )
}


# Coordinate-Exchange algorithm (Bayesian D-optimality) with fast rank-2 updates
# - X_init: initial design matrix (n x q)
# - X_all : full candidate set (N x q)
# - R     : prior precision (q x q)
# - coord_cols: integer vector of column indices on which to perform coordinate updates.
# - n_repeat: number of full sweeps through all rows/coordinates
coordinate_exchange_bayes <- function(X_init, X_all, R, coord_cols = NULL, n_repeat = 5) {
  
  X_best <- X_init
  n <- nrow(X_best)
  q <- ncol(X_best)
  N <- nrow(X_all)
  
  if (is.null(coord_cols)) coord_cols <- seq_len(q)  
  coord_cols <- sort(unique(coord_cols))
  
  # Initialize A, its inverse, and Bayesian D-criterion (per-parameter geometric mean)
  A_best <- crossprod(X_best) / n + R
  det_best <- det(A_best)
  phi_best <- det_best^(1 / q)
  A_best_inv <- solve(A_best)
  
  # Helper to compute det update via rank-2 replacement of row i: x -> y
  det_after_swap <- function(Ainv, detA, x_old, x_new) {
    fx  <- matrix(x_old, ncol = 1) / sqrt(n)
    fy  <- matrix(x_new, ncol = 1) / sqrt(n)
    vx  <- as.numeric(t(fx) %*% Ainv %*% fx)
    vy  <- as.numeric(t(fy) %*% Ainv %*% fy)
    vxy <- as.numeric(t(fx) %*% Ainv %*% fy)
    detA * (1 + vy - vx + vxy^2 - vy * vx)
  }
  
  # Helper to update A_inv after accepting a swap (rank-2 Woodbury)
  update_Ainv_after_swap <- function(Ainv, x_old, x_new) {
    fx  <- matrix(x_old, ncol = 1) / sqrt(n)
    fy  <- matrix(x_new, ncol = 1) / sqrt(n)
    F1  <- cbind(fy, -fx)
    F2  <- cbind(fy,  fx)
    tmp <- solve(diag(2) + t(F2) %*% Ainv %*% F1)  
    Ainv - Ainv %*% F1 %*% tmp %*% t(F2) %*% Ainv
  }
  
  for (rep in seq_len(n_repeat)) {
    improved_any <- FALSE
    
    # Loop over runs (rows) and coordinates (columns)
    for (i in seq_len(n)) {
      x_i <- X_best[i, ]
      
      for (k in coord_cols) {
        # Find candidates in X_all that differ from x_i in exactly the kth coordinate
        # (i.e., all other q-1 columns equal). This enforces coordinate exchange.
        same_except_k <- rowSums(X_all[, -k, drop = FALSE] ==
                                   matrix(x_i[-k], nrow = N, ncol = q - 1, byrow = TRUE)) == (q - 1)
        
        cand_idx <- which(same_except_k)
        
        # If nothing differs only on k, skip this coordinate
        if (length(cand_idx) == 0) next
        
        # Evaluate all candidates with fast determinant update
        det_cands <- numeric(length(cand_idx))
        for (t in seq_along(cand_idx)) {
          j      <- cand_idx[t]
          det_j  <- det_after_swap(A_best_inv, det_best, x_old = x_i, x_new = X_all[j, ])
          det_cands[t] <- if (is.finite(det_j) && det_j > 0) det_j else 0
        }
        
        best_local_idx <- which.max(det_cands)
        if (length(best_local_idx) == 0) next
        
        det_best_local <- det_cands[best_local_idx]
        phi_best_local <- det_best_local^(1 / q)
        
        # Accept only if strictly better
        if (phi_best_local > phi_best) {
          j_star <- cand_idx[best_local_idx]
          
          # Update 
          A_best_inv <- update_Ainv_after_swap(A_best_inv, x_old = x_i, x_new = X_all[j_star, ])
          X_best[i, ] <- X_all[j_star, ]
          det_best    <- det_best_local
          phi_best    <- phi_best_local
          x_i         <- X_best[i, ]     
          improved_any <- TRUE
        }
      }
    }
    
    # Optional early stop if a full sweep made no improvement
    if (!improved_any) {
      print("stopped early")
      break
    }
  }
  
  list(model.matrix = X_best, phi = phi_best)
}

