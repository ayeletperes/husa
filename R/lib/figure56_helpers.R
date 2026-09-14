# Heatmap helpers shared by figure 4, supp 4 and supp 5: the blockiness-maximising
# dendrogram flip optimiser, the phylo clade collapse, the allele-group / usage builders,
# and the ComplexHeatmap annotation pieces (consensus-difference row labels, subgroup
# glyphs, usage boxplots, count colour ramp). Lifted from the manuscript tree unchanged.

suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(grid); library(circlize); library(Biostrings)
  library(gridtext); library(GetoptLong); library(digest); library(mallinfo); library(dendextend); library(ape); library(jsonlite)
})

########################################

# Optimized Matrix Ordering Code - Complete Implementation

## Helper Functions


#### Helper functions ############

`%||%` <- function(x, y) if (is.null(x)) y else x
address <- function(x) paste0(attr(x, "class"), "_", digest::digest(x))


get_family <- function(vec) {
  family_vec <- sapply(strsplit(vec, "-"), "[[", 1)
  # remove D
  family_vec <- gsub("D", "", family_vec)
  names(family_vec) <- vec
  return(family_vec)
}

#### Optimized blockiness score function ##########

#' Generate block kernel with exponentially decreasing weights
#' @param size size of the kernel
#' @return block kernel
generate_block_kernel <- function(size) {
  if (size %% 2 != 1) stop("Kernel size must be odd.")
  center <- ceiling(size / 2)
  outer <- abs(row(matrix(1, size, size)) - center) + abs(col(matrix(1, size, size)) - center)
  kernel <- 2 / 2^outer
  kernel
}

#' Calculate blockiness score (optimized version)
#' @param mat_ matrix to calculate blockiness for
#' @param kernel_size size of the kernel
#' @param custom_kernel custom kernel to use
#' @return blockiness score
calc_blockiness_optimized <- function(mat_, kernel_size = 5, custom_kernel = NULL) {
  kernel <- if (is.null(custom_kernel)) generate_block_kernel(kernel_size) else custom_kernel
  kcenter <- ceiling(kernel_size / 2)
  mat_binary <- (mat_ > 0) * 1
  
  # Early exit for empty matrix
  if (sum(mat_binary) == 0) return(0)
  
  # Pad the matrix to handle boundaries efficiently
  pad_size <- kcenter - 1
  padded_mat <- matrix(0, 
                       nrow = nrow(mat_binary) + 2 * pad_size,
                       ncol = ncol(mat_binary) + 2 * pad_size)
  padded_mat[(pad_size + 1):(nrow(padded_mat) - pad_size),
             (pad_size + 1):(ncol(padded_mat) - pad_size)] <- mat_binary
  
  # Use 2D convolution for blockiness calculation
  blockiness_score <- 0
  non_zero <- which(mat_binary > 0, arr.ind = TRUE)
  
  # Vectorized computation for all non-zero cells
  for (i in 1:nrow(non_zero)) {
    row_idx <- non_zero[i, 1] + pad_size
    col_idx <- non_zero[i, 2] + pad_size
    
    # Extract neighborhood and compute weighted sum
    neighborhood <- padded_mat[(row_idx - pad_size):(row_idx + pad_size),
                               (col_idx - pad_size):(col_idx + pad_size)]
    blockiness_score <- blockiness_score + sum(neighborhood * kernel)
  }
  
  return(blockiness_score)
}

#' Create a cache for blockiness calculations
#' @return cache object with calc and clear functions
create_blockiness_cache <- function() {
  cache <- new.env(hash = TRUE, parent = emptyenv())
  
  get_cache_key <- function(row_order, col_order, subgroup = NULL) {
    jsonlite::toJSON(
      list(rows = unname(row_order), cols = unname(col_order), subgroup = subgroup),
      auto_unbox = TRUE,
      null = "null"
    )
  }
  
  calc_blockiness_cached <- function(mat_subset, row_order, col_order, kernel_size = 5, subgroup = NULL) {
    row_order <- row_order[row_order %in% rownames(mat_subset)]
    col_order <- col_order[col_order %in% colnames(mat_subset)]
    
    key <- get_cache_key(row_order, col_order, subgroup)
    
    if (exists(key, envir = cache)) {
      return(get(key, envir = cache))
    }
    
    # Reorder matrix and calculate
    mat_ordered <- mat_subset[row_order, col_order, drop = FALSE]
    score <- calc_blockiness_optimized(mat_ordered, kernel_size)
    
    # Store in cache
    assign(key, score, envir = cache)
    return(score)
  }
  
  clear_cache <- function() {
    # Remove all objects
    obj_names <- ls(envir = cache, all.names = T)
    rm(list = obj_names, envir = cache)
    # Force garbage collection
    mallinfo::malloc.trim(0L)
    invisible(NULL)
  }
  
  finalize_cache <- function() {
    clear_cache()
    if (exists("cache") && !is.null(cache)) {
      rm(list = ls(envir = cache, all.names = TRUE), envir = cache)
      cache <<- NULL  # Use superassignment to ensure it's cleared
    }
    invisible(NULL)
  }
  
  list(
    calc = calc_blockiness_cached,
    clear = clear_cache,
    finalize = finalize_cache,
    size = function() length(ls(envir = cache))
  )
}

apply_flip_sequence <- function(dend, flip_sequence) {
  out <- dend
  if (length(flip_sequence) == 0L) {
    return(out)
  }
  for (flip in flip_sequence) {
    out <- flip_node_fast(out, flip$node_id)
  }
  out
}

score_by_order <- function(mat_subset, row_order, col_order, kernel_size, cache = NULL, subgroup = NULL) {
  row_order <- row_order[row_order %in% rownames(mat_subset)]
  col_order <- col_order[col_order %in% colnames(mat_subset)]
  if (!length(row_order) || !length(col_order)) {
    return(0)
  }
  if (!is.null(cache)) {
    return(cache$calc(mat_subset, row_order, col_order, kernel_size, subgroup))
  }
  calc_blockiness_optimized(mat_subset[row_order, col_order, drop = FALSE], kernel_size)
}

#' Add node IDs to a dendrogram
#' @param dend dendrogram
#' @param start_id starting ID
#' @return dendrogram with node IDs
add_node_ids <- function(dend, start_id = 1) {
  id_counter <- start_id
  
  add_ids_recursive <- function(d) {
    if (is.leaf(d)) {
      return(d)
    }
    
    attr(d, "id") <- id_counter
    id_counter <<- id_counter + 1
    
    d[[1]] <- add_ids_recursive(d[[1]])
    d[[2]] <- add_ids_recursive(d[[2]])
    
    return(d)
  }
  
  add_ids_recursive(dend)
}

#' Flip a node in a dendrogram
#' @param dend dendrogram
#' @param node_id ID of the node to flip
#' @return the flipped dendrogram
# flip_node_fast <- function(dend, node_id) {
#   if (is.leaf(dend)) {
#     return(dend)
#   }
  
#   current_id <- attr(dend, "id")
#   if (!is.null(current_id) && current_id == node_id) {
#     temp <- dend[[1]]
#     dend[[1]] <- dend[[2]]
#     dend[[2]] <- temp
#     return(dend)
#   }
  
#   if (!is.leaf(dend[[1]])) {
#     dend[[1]] <- flip_node_fast(dend[[1]], node_id)
#   }
#   if (!is.leaf(dend[[2]])) {
#     dend[[2]] <- flip_node_fast(dend[[2]], node_id)
#   }
  
#   return(dend)
# }


flip_node_fast <- function(dend, node_id, label_cache = NULL) {
  if (is.leaf(dend)) return(dend)
  cur_id <- attr(dend,"id")
  if (!is.null(cur_id) && cur_id == node_id) {
    tmp <- dend[[1]]; dend[[1]] <- dend[[2]]; dend[[2]] <- tmp
    return(dend)
  }
  # descend only if needed
  if (!is.leaf(dend[[1]])) dend[[1]] <- flip_node_fast(dend[[1]], node_id, label_cache)
  if (!is.leaf(dend[[2]])) dend[[2]] <- flip_node_fast(dend[[2]], node_id, label_cache)
  dend
}


#' Get all internal node IDs from a dendrogram
#' @param dend dendrogram
#' @return vector of internal node IDs
get_internal_node_ids_fast <- function(dend) {
  ids <- integer()
  stack <- list(dend)
  
  while (length(stack) > 0) {
    current <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    
    if (!is.leaf(current)) {
      node_id <- attr(current, "id")
      if (!is.null(node_id)) {
        ids <- c(ids, node_id)
      }
      
      if (!is.leaf(current[[2]])) stack[[length(stack) + 1]] <- current[[2]]
      if (!is.leaf(current[[1]])) stack[[length(stack) + 1]] <- current[[1]]
    }
  }
  
  return(ids)
}

#' Build a node index for fast lookups
#' @param dend dendrogram
#' @return named list where names are node IDs and values are nodes
build_node_index <- function(dend) {
  node_index <- list()
  stack <- list(dend)
  
  while (length(stack) > 0) {
    current <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    
    if (!is.leaf(current)) {
      node_id <- attr(current, "id")
      if (!is.null(node_id)) {
        node_index[[as.character(node_id)]] <- current
      }
      
      # Add children to stack
      if (!is.leaf(current[[2]])) {
        stack[[length(stack) + 1]] <- current[[2]]
      }
      if (!is.leaf(current[[1]])) {
        stack[[length(stack) + 1]] <- current[[1]]
      }
    }
  }
  
  return(node_index)
}

#' Fast node lookup using pre-built index
#' @param node_index pre-built index from build_node_index()
#' @param target_id ID to find
#' @return node or NULL
get_node_by_id_fast <- function(node_index, target_id) {
  return(node_index[[as.character(target_id)]])
}

#' Check if a dendrogram node contains any of the specified labels
#' @param node dendrogram node
#' @param target_labels labels to check for
#' @return TRUE if node contains any target labels
contains_any_labels <- function(node, target_labels) {
  node_labels <- labels(node)
  any(node_labels %in% target_labels)
}

#' Get the smallest subtree containing all specified labels
#' This finds the minimal subtree that includes all family labels
#' @param dend full dendrogram
#' @param target_labels labels that must be in the subtree
#' @return minimal subtree containing all target labels
# find_minimal_subtree <- function(dend, target_labels) {
#   if (is.leaf(dend)) {
#     if (attr(dend, "label") %in% target_labels) {
#       return(dend)
#     } else {
#       return(NULL)
#     }
#   }
  
#   # Check if current node contains all target labels
#   current_labels <- labels(dend)
#   contains_targets <- any(current_labels %in% target_labels)
  
#   if (!contains_targets) {
#     return(NULL)
#   }
  
#   # Check children
#   left_contains <- contains_any_labels(dend[[1]], target_labels)
#   right_contains <- contains_any_labels(dend[[2]], target_labels)
  
#   # If all targets are in one child, recurse into that child
#   left_labels <- if (left_contains) labels(dend[[1]]) else character(0)
#   right_labels <- if (right_contains) labels(dend[[2]]) else character(0)
  
#   targets_in_left <- sum(target_labels %in% left_labels)
#   targets_in_right <- sum(target_labels %in% right_labels)
#   total_targets <- length(target_labels)
  
#   if (targets_in_left == total_targets) {
#     return(find_minimal_subtree(dend[[1]], target_labels))
#   } else if (targets_in_right == total_targets) {
#     return(find_minimal_subtree(dend[[2]], target_labels))
#   } else {
#     # Targets are split between children, so this is the minimal subtree
#     return(dend)
#   }
# }

find_minimal_subtree <- function(dend, target_labels) {
  # mark presence bottom-up
  postorder <- list()
  stack <- list(list(node=dend, visited=FALSE))
  has_target <- new.env(parent = emptyenv())

  while (length(stack)) {
    cur <- stack[[length(stack)]]; stack <- stack[-length(stack)]
    node <- cur$node
    if (is.leaf(node)) {
      lab <- attr(node,"label")
      has_target[[address(node)]] <- lab %in% target_labels
      postorder[[length(postorder)+1]] <- node
    } else if (!cur$visited) {
      stack[[length(stack)+1]] <- list(node=node, visited=TRUE)
      stack[[length(stack)+1]] <- list(node=node[[2]], visited=FALSE)
      stack[[length(stack)+1]] <- list(node=node[[1]], visited=FALSE)
    } else {
      a <- has_target[[address(node[[1]])]] %||% FALSE
      b <- has_target[[address(node[[2]])]] %||% FALSE
      has_target[[address(node)]] <- a || b
      postorder[[length(postorder)+1]] <- node
    }
  }

  # climb down to the smallest node that still contains any target
  cur <- dend
  repeat {
    if (is.leaf(cur)) return(cur)
    a <- has_target[[address(cur[[1]])]] %||% FALSE
    b <- has_target[[address(cur[[2]])]] %||% FALSE
    if (a && !b) { cur <- cur[[1]]; next }
    if (b && !a) { cur <- cur[[2]]; next }
    # split across both children -> current is minimal
    return(cur)
  }
}


build_node_index_with_labels <- function(dend) {
  node_index <- list()
  label_cache <- list()
  stack <- list(dend)
  while (length(stack)) {
    cur <- stack[[length(stack)]]; stack <- stack[-length(stack)]
    if (!is.leaf(cur)) {
      id <- attr(cur, "id")
      if (!is.null(id)) node_index[[as.character(id)]] <- cur

      # collect labels for the two children without calling labels() on the parent
      left  <- cur[[1]]; right <- cur[[2]]
      if (!is.leaf(left))  stack[[length(stack)+1]] <- left
      if (!is.leaf(right)) stack[[length(stack)+1]] <- right

      left_labels  <- if (is.leaf(left))  attr(left, "label")  else NULL
      right_labels <- if (is.leaf(right)) attr(right, "label") else NULL

      # climb children iteratively to get all labels, but cache as we go
      get_all_labels <- function(node) {
        if (is.leaf(node)) return(attr(node, "label"))
        nid <- attr(node, "id")
        if (!is.null(nid) && !is.null(label_cache[[as.character(nid)]])) {
          return(label_cache[[as.character(nid)]])
        }
        q <- list(node)
        out <- character()
        while (length(q)) {
          n <- q[[length(q)]]; q <- q[-length(q)]
          if (is.leaf(n)) {
            out <- c(out, attr(n,"label"))
          } else {
            # prefer cache if present
            nid2 <- attr(n,"id")
            if (!is.null(nid2) && !is.null(label_cache[[as.character(nid2)]])) {
              out <- c(out, label_cache[[as.character(nid2)]])
            } else {
              q[[length(q)+1]] <- n[[1]]
              q[[length(q)+1]] <- n[[2]]
            }
          }
        }
        # store on this node if it has an id
        nid <- attr(node,"id")
        if (!is.null(nid)) label_cache[[as.character(nid)]] <<- out
        out
      }

      # cache labels for this internal node
      if (!is.null(id)) {
        label_cache[[as.character(id)]] <- get_all_labels(cur)
      }
    }
  }
  list(node_index = node_index, label_cache = label_cache)
}


#' Get potential flips that affect family labels
#' Works with full dendrogram but only returns flips affecting target labels
#' @param dend full dendrogram (not pruned)
#' @param labels_to_keep family labels we care about
#' @param consider_mixed_flips if TRUE, include flips that mix family and non-family labels
#' @param ... additional arguments (not used)
#' @return list of possible flips
# get_potential_flips <- function(dend, labels_to_keep, consider_mixed_flips = TRUE) {
#   if (length(labels_to_keep) <= 1) {
#     return(NULL)
#   }
  
#   # Work with the full dendrogram
#   node_index <- build_node_index(dend)
#   internal_nodes <- get_internal_node_ids_fast(dend)
#   possible_flips <- list()
  
#   for (node_id in internal_nodes) {
#     node <- get_node_by_id_fast(node_index, node_id)
#     if (!is.null(node) && length(node) == 2) {
#       # Check if this flip would affect any family labels
#       left_labels <- labels(node[[1]])
#       right_labels <- labels(node[[2]])
      
#       left_has_family <- any(left_labels %in% labels_to_keep)
#       right_has_family <- any(right_labels %in% labels_to_keep)
      
#       # Skip if neither side has family labels
#       if (!left_has_family && !right_has_family) {
#         next
#       }
      
#       # Option to skip mixed flips (family with non-family)
#       if (!consider_mixed_flips) {
#         left_all_family <- all(left_labels %in% labels_to_keep)
#         right_all_family <- all(right_labels %in% labels_to_keep)
        
#         # Skip if one side is pure family and other is pure non-family
#         if ((left_all_family && !right_has_family) || 
#             (right_all_family && !left_has_family)) {
#           next
#         }
#       }
      
#       # Create flip description
#       if (is.leaf(node[[1]]) && is.leaf(node[[2]])) {
#         left_label <- attr(node[[1]], "label")
#         right_label <- attr(node[[2]], "label")
        
#         # Mark family labels for clarity
#         left_mark <- if (left_label %in% labels_to_keep) paste0("[F]", left_label) else left_label
#         right_mark <- if (right_label %in% labels_to_keep) paste0("[F]", right_label) else right_label
        
#         flip_desc <- paste(left_mark, "↔", right_mark)
        
#       } else if (is.leaf(node[[1]]) && !is.leaf(node[[2]])) {
#         left_label <- attr(node[[1]], "label")
#         left_mark <- if (left_label %in% labels_to_keep) paste0("[F]", left_label) else left_label
        
#         # Summarize right subtree
#         right_family <- sum(right_labels %in% labels_to_keep)
#         right_other <- length(right_labels) - right_family
#         right_desc <- if (right_other > 0) {
#           paste0("(", right_family, "F+", right_other, "O)")
#         } else {
#           paste0("(", right_family, "F)")
#         }
        
#         flip_desc <- paste(left_mark, "↔", right_desc)
        
#       } else if (!is.leaf(node[[1]]) && is.leaf(node[[2]])) {
#         right_label <- attr(node[[2]], "label")
#         right_mark <- if (right_label %in% labels_to_keep) paste0("[F]", right_label) else right_label
        
#         # Summarize left subtree
#         left_family <- sum(left_labels %in% labels_to_keep)
#         left_other <- length(left_labels) - left_family
#         left_desc <- if (left_other > 0) {
#           paste0("(", left_family, "F+", left_other, "O)")
#         } else {
#           paste0("(", left_family, "F)")
#         }
        
#         flip_desc <- paste(left_desc, "↔", right_mark)
        
#       } else {
#         # Both are subtrees - summarize each
#         left_family <- sum(left_labels %in% labels_to_keep)
#         left_other <- length(left_labels) - left_family
#         left_desc <- if (left_other > 0) {
#           paste0("(", left_family, "F+", left_other, "O)")
#         } else {
#           paste0("(", left_family, "F)")
#         }
        
#         right_family <- sum(right_labels %in% labels_to_keep)
#         right_other <- length(right_labels) - right_family
#         right_desc <- if (right_other > 0) {
#           paste0("(", right_family, "F+", right_other, "O)")
#         } else {
#           paste0("(", right_family, "F)")
#         }
        
#         flip_desc <- paste(left_desc, "↔", right_desc)
#       }
      
#       # Calculate impact score (how many family labels are affected)
#       family_affected <- sum(c(left_labels, right_labels) %in% labels_to_keep)
      
#       possible_flips[[length(possible_flips) + 1]] <- list(
#         node_id = node_id,
#         description = flip_desc,
#         node = node,
#         family_labels_affected = family_affected,
#         total_labels_affected = length(c(left_labels, right_labels)),
#         left_has_family = left_has_family,
#         right_has_family = right_has_family
#       )
#     }
#   }
  
#   # Sort flips by impact (prioritize flips affecting more family labels)
#   if (length(possible_flips) > 0) {
#     impact_scores <- sapply(possible_flips, function(x) x$family_labels_affected)
#     possible_flips <- possible_flips[order(impact_scores, decreasing = TRUE)]
#   }
  
#   return(possible_flips)
# }

get_potential_flips <- function(dend, labels_to_keep, consider_mixed_flips = TRUE) {
  if (length(labels_to_keep) <= 1) return(NULL)

  idx <- build_node_index_with_labels(dend)
  node_index  <- idx$node_index
  label_cache <- idx$label_cache
  internal_ids <- as.integer(names(node_index))
  possible <- vector("list", length(internal_ids))
  keep_set <- unique(labels_to_keep)

  k <- 0
  for (node_id in internal_ids) {
    node <- node_index[[as.character(node_id)]]
    if (is.null(node) || length(node) != 2) next

    left  <- node[[1]]
    right <- node[[2]]

    left_labels  <- if (is.leaf(left))  attr(left,"label")  else label_cache[[as.character(attr(left,"id"))]]
    right_labels <- if (is.leaf(right)) attr(right,"label") else label_cache[[as.character(attr(right,"id"))]]

    # safety
    if (is.null(left_labels))  left_labels  <- character()
    if (is.null(right_labels)) right_labels <- character()

    left_has  <- any(left_labels  %in% keep_set)
    right_has <- any(right_labels %in% keep_set)
    if (!left_has && !right_has) next

    if (!consider_mixed_flips) {
      left_all  <- length(left_labels)  && all(left_labels  %in% keep_set)
      right_all <- length(right_labels) && all(right_labels %in% keep_set)
      if ((left_all && !right_has) || (right_all && !left_has)) next
    }

    desc <- if (is.leaf(left) && is.leaf(right)) {
      paste(if (attr(left,"label")  %in% keep_set)  paste0("[F]", attr(left,"label"))  else attr(left,"label"),
            "↔",
            if (attr(right,"label") %in% keep_set) paste0("[F]", attr(right,"label")) else attr(right,"label"))
    } else {
      lf <- sum(left_labels  %in% keep_set);  lo <- length(left_labels)  - lf
      rf <- sum(right_labels %in% keep_set); ro <- length(right_labels) - rf
      paste0("(", lf, "F+", lo, "O) ↔ (", rf, "F+", ro, "O)")
    }

    k <- k + 1
    possible[[k]] <- list(
      node_id = node_id,
      description = desc,
      family_labels_affected = sum(c(left_labels, right_labels) %in% keep_set),
      total_labels_affected  = length(c(left_labels, right_labels)),
      left_has_family = left_has,
      right_has_family = right_has
    )
  }
  if (k == 0) return(NULL)
  possible <- possible[seq_len(k)]
  possible[order(vapply(possible, function(x) x$family_labels_affected, 1L), decreasing = TRUE)]
}



#' Alternative version that focuses on a minimal subtree
#' This finds the smallest subtree containing all family labels and works within it
#' @param dend dendrogram
#' @param labels_to_keep labels to keep
#' @param ... additional arguments (not used)
#' @return list of potential flips
get_potential_flips_minimal_subtree <- function(dend, labels_to_keep, ...) {
  if (length(labels_to_keep) <= 1) {
    return(NULL)
  }
  
  # Find the minimal subtree containing all family labels
  minimal_subtree <- find_minimal_subtree(dend, labels_to_keep)
  
  if (is.null(minimal_subtree) || is.leaf(minimal_subtree)) {
    return(NULL)
  }
  
  # Get all internal nodes in the minimal subtree
  node_index <- build_node_index(minimal_subtree)
  internal_nodes <- get_internal_node_ids_fast(minimal_subtree)
  possible_flips <- list()
  
  for (node_id in internal_nodes) {
    node <- get_node_by_id_fast(node_index, node_id)
    if (!is.null(node) && length(node) == 2) {
      left_labels <- labels(node[[1]])
      right_labels <- labels(node[[2]])
      
      # Create detailed flip description
      create_label_summary <- function(labels, target_labels) {
        in_family <- labels[labels %in% target_labels]
        out_family <- labels[!labels %in% target_labels]
        
        if (length(labels) == 1) {
          if (labels %in% target_labels) {
            return(paste0("[F]", labels))
          } else {
            return(labels)
          }
        } else {
          family_count <- length(in_family)
          other_count <- length(out_family)
          
          if (family_count <= 3 && other_count == 0) {
            # Show individual family labels if few
            return(paste0("[F:", paste(in_family, collapse = ","), "]"))
          } else if (other_count > 0) {
            return(paste0("(", family_count, "F+", other_count, "O)"))
          } else {
            return(paste0("(", family_count, "F)"))
          }
        }
      }
      
      left_desc <- create_label_summary(left_labels, labels_to_keep)
      right_desc <- create_label_summary(right_labels, labels_to_keep)
      flip_desc <- paste(left_desc, "↔", right_desc)
      
      # Calculate impact metrics
      family_in_left <- sum(left_labels %in% labels_to_keep)
      family_in_right <- sum(right_labels %in% labels_to_keep)
      
      possible_flips[[length(possible_flips) + 1]] <- list(
        node_id = node_id,
        description = flip_desc,
        node = node,
        family_labels_affected = family_in_left + family_in_right,
        total_labels_affected = length(left_labels) + length(right_labels),
        family_in_left = family_in_left,
        family_in_right = family_in_right,
        other_in_left = length(left_labels) - family_in_left,
        other_in_right = length(right_labels) - family_in_right
      )
    }
  }
  
  return(possible_flips)
}

#' Pre-order dendrogram to group subgroups together
#' @param dend dendrogram
#' @return dendrogram with subgroups clustered together
reorder_by_subgroups <- function(dend) {
  all_labels <- labels(dend)
  subgroup_map <- get_subgroup_membership(all_labels, labels_to_keep)
  
  # Group labels by subgroup
  subgroup_groups <- split(all_labels, subgroup_map)
  
  # Create desired order: subgroups together, maintaining internal structure
  desired_order <- unlist(subgroup_groups)
  
  # Find the order indices
  current_order <- labels(dend)
  order_indices <- match(desired_order, current_order)
  
  # Reorder the dendrogram
  dend_reordered <- reorder(dend, order_indices, agglo.FUN = mean)
  
  return(dend_reordered)
}

#' Get subgroup membership for sequences
#' @param labels sequence labels
#' @return named vector of subgroup assignments
get_subgroup_membership <- function(labels, labels_to_keep) {
  # Since labels_to_keep are already the labels for this specific subgroup/family,
  # we just need to mark which labels belong to this group vs others
  
  subgroups <- character(length(labels))
  names(subgroups) <- labels
  
  # Labels in labels_to_keep belong to the current subgroup
  subgroups[labels_to_keep] <- "CURRENT_SUBGROUP"
  # All other labels belong to different subgroups
  subgroups[!labels %in% labels_to_keep] <- "OTHER_SUBGROUP"
  
  return(subgroups)
}

#' Check if a flip respects subgroup boundaries
#' @param node dendrogram node
#' @param subgroup_map named vector of subgroup assignments
#' @return TRUE if flip respects subgroup boundaries
flip_respects_subgroups <- function(node, subgroup_map) {
  if (is.leaf(node)) return(TRUE)
  
  left_labels <- labels(node[[1]])
  right_labels <- labels(node[[2]])
  
  # Get subgroups for each side
  left_subgroups <- unique(subgroup_map[left_labels])
  right_subgroups <- unique(subgroup_map[right_labels])
  
  # Check if there's overlap in subgroups
  overlap <- intersect(left_subgroups, right_subgroups)
  
  # If there's no overlap, the flip respects boundaries
  return(length(overlap) == 0)
}

#' Modified version of get_potential_flips that respects subgroup boundaries
#' @param dend dendrogram
#' @param labels_to_keep labels to keep
#' @param consider_mixed_flips if TRUE, include flips that mix family and non-family labels
#' @param respect_subgroups if TRUE, only allow flips that don't mix subgroups
#' @return list of potential flips that respect subgroup boundaries
get_potential_flips_subgroup_aware <- function(dend, labels_to_keep, 
                                               consider_mixed_flips = TRUE,
                                               respect_subgroups = TRUE) {
  if (length(labels_to_keep) <= 1) {
    return(NULL)
  }
  
  # Get subgroup membership
  all_labels <- labels(dend)
  subgroup_map <- get_subgroup_membership(all_labels, labels_to_keep)
  
  # Work with the full dendrogram
  node_index <- build_node_index_with_labels(dend)
  internal_nodes <- get_internal_node_ids_fast(dend)
  possible_flips <- list()
  
  for (node_id in internal_nodes) {
    node <- get_node_by_id_fast(node_index, node_id)
    if (!is.null(node) && length(node) == 2) {
      # Check if this flip would affect any family labels
      left_labels <- labels(node[[1]])
      right_labels <- labels(node[[2]])
      
      left_has_family <- any(left_labels %in% labels_to_keep)
      right_has_family <- any(right_labels %in% labels_to_keep)
      
      # Skip if neither side has family labels
      if (!left_has_family && !right_has_family) {
        next
      }
      # Check subgroup boundaries if requested
      if (respect_subgroups && !flip_respects_subgroups(node, subgroup_map)) {
        next  # Skip flips that would mix subgroups
      }
      
      # Option to skip mixed flips (family with non-family)
      if (!consider_mixed_flips) {
        left_all_family <- all(left_labels %in% labels_to_keep)
        right_all_family <- all(right_labels %in% labels_to_keep)
        # Skip if one side is pure family and other is pure non-family
        if ((left_all_family && !right_has_family) || 
            (right_all_family && !left_has_family)) {
          next
        }
      }
      
      # Create flip description (same as before)
      if (is.leaf(node[[1]]) && is.leaf(node[[2]])) {
        left_label <- attr(node[[1]], "label")
        right_label <- attr(node[[2]], "label")
        
        # Mark family labels for clarity
        left_mark <- if (left_label %in% labels_to_keep) paste0("[F]", left_label) else left_label
        right_mark <- if (right_label %in% labels_to_keep) paste0("[F]", right_label) else right_label
        
        flip_desc <- paste(left_mark, "↔", right_mark)
        
      } else if (is.leaf(node[[1]]) && !is.leaf(node[[2]])) {
        left_label <- attr(node[[1]], "label")
        left_mark <- if (left_label %in% labels_to_keep) paste0("[F]", left_label) else left_label
        
        # Summarize right subtree
        right_family <- sum(right_labels %in% labels_to_keep)
        right_other <- length(right_labels) - right_family
        right_desc <- if (right_other > 0) {
          paste0("(", right_family, "F+", right_other, "O)")
        } else {
          paste0("(", right_family, "F)")
        }
        
        flip_desc <- paste(left_mark, "↔", right_desc)
        
      } else if (!is.leaf(node[[1]]) && is.leaf(node[[2]])) {
        right_label <- attr(node[[2]], "label")
        right_mark <- if (right_label %in% labels_to_keep) paste0("[F]", right_label) else right_label
        
        # Summarize left subtree
        left_family <- sum(left_labels %in% labels_to_keep)
        left_other <- length(left_labels) - left_family
        left_desc <- if (left_other > 0) {
          paste0("(", left_family, "F+", left_other, "O)")
        } else {
          paste0("(", left_family, "F)")
        }
        
        flip_desc <- paste(left_desc, "↔", right_mark)
        
      } else {
        # Both are subtrees - summarize each
        left_family <- sum(left_labels %in% labels_to_keep)
        left_other <- length(left_labels) - left_family
        left_desc <- if (left_other > 0) {
          paste0("(", left_family, "F+", left_other, "O)")
        } else {
          paste0("(", left_family, "F)")
        }
        
        right_family <- sum(right_labels %in% labels_to_keep)
        right_other <- length(right_labels) - right_family
        right_desc <- if (right_other > 0) {
          paste0("(", right_family, "F+", right_other, "O)")
        } else {
          paste0("(", right_family, "F)")
        }
        
        flip_desc <- paste(left_desc, "↔", right_desc)
      }
      
      # Calculate impact score (how many family labels are affected)
      family_affected <- sum(c(left_labels, right_labels) %in% labels_to_keep)
      
      possible_flips[[length(possible_flips) + 1]] <- list(
        node_id = node_id,
        description = flip_desc,
        node = node,
        family_labels_affected = family_affected,
        total_labels_affected = length(c(left_labels, right_labels)),
        left_has_family = left_has_family,
        right_has_family = right_has_family
      )
    }
    
    if (exists("left_labels")) rm(left_labels)
    if (exists("right_labels")) rm(right_labels)
  }
  
  rm(internal_nodes)
  # Sort flips by impact (prioritize flips affecting more family labels)
  if (length(possible_flips) > 0) {
    impact_scores <- sapply(possible_flips, function(x) x$family_labels_affected)
    possible_flips <- possible_flips[order(impact_scores, decreasing = TRUE)]
  }
  
  return(possible_flips)
}

#' Generate promising sequences for a dendrogram
#' @param dend dendrogram
#' @param labels_to_keep labels to keep
#' @param n_steps number of steps
#' @param mat_subset matrix subset
#' @param kernel_size kernel size
#' @param cache cache object
#' @param max_sequences maximum number of sequences to generate
#' @param consider_mixed_flips if TRUE, include flips that mix family and non-family labels
#' @param consider_minimal_subtree if TRUE, consider flips that affect a minimal subtree containing all family labels
#' @return list of promising sequences
generate_promising_sequences <- function(dend, labels_to_keep, n_steps, 
                                         mat_subset, kernel_size,
                                         cache = NULL,
                                         max_sequences = 100,
                                         consider_mixed_flips = TRUE,
                                         consider_minimal_subtree = FALSE,
                                         consider_subgroups = TRUE) {
  if (n_steps == 0) return(list())
  
  flip_function <- if (consider_minimal_subtree) get_potential_flips_minimal_subtree else get_potential_flips_subgroup_aware
  
  # Start with single-step sequences
  current_sequences <- list()
  
  # Get initial flips
  possible_flips <- flip_function(dend, labels_to_keep, consider_mixed_flips, consider_subgroups)
  if (length(possible_flips) == 0) return(list())
  
  # Initialize with single flips
  for (flip in possible_flips) {
    current_sequences <- c(current_sequences, list(list(
      sequence = list(flip)
    )))
  }
  rm(possible_flips) # Clean up immediately
  
  # Iteratively extend sequences
  for (step in 2:n_steps) {
    if (length(current_sequences) == 0) break
    
    next_sequences <- list()
    
    # Limit the number of sequences to extend (beam search)
    beam_width <- min(10, length(current_sequences))
    
    # Score current sequences and keep best ones
    if (length(current_sequences) > beam_width) {
      scores <- numeric(length(current_sequences))
      for (i in seq_along(current_sequences)) {
        seq_info <- current_sequences[[i]]
        scored_dend <- apply_flip_sequence(dend, seq_info$sequence)
        new_order <- labels(scored_dend)
        scores[i] <- score_by_order(
          mat_subset = mat_subset,
          row_order = new_order,
          col_order = colnames(mat_subset),
          kernel_size = kernel_size,
          cache = cache
        )
      }
      
      # Keep only the best sequences
      best_indices <- order(scores, decreasing = TRUE)[1:beam_width]
      current_sequences <- current_sequences[best_indices]
      rm(scores, best_indices)
    }
    
    # Extend each current sequence
    for (seq_info in current_sequences) {
      # Get possible flips from current dendrogram
      current_dend <- apply_flip_sequence(dend, seq_info$sequence)
      possible_flips <- flip_function(current_dend, labels_to_keep, 
                                      consider_mixed_flips, consider_subgroups)
      
      if (length(possible_flips) > 0) {
        # Limit extensions per sequence
        max_extensions <- min(3, length(possible_flips))
        
        for (i in 1:max_extensions) {
          flip <- possible_flips[[i]]
          new_sequence <- c(seq_info$sequence, list(flip))
          
          next_sequences <- c(next_sequences, list(list(
            sequence = new_sequence
          )))
          
          # Early exit if we have enough sequences
          if (length(next_sequences) >= max_sequences) {
            rm(possible_flips)
            break
          }
        }
      }
      
      rm(possible_flips)
      
      if (length(next_sequences) >= max_sequences) break
    }
    
    # Clean up current sequences before moving to next
    rm(current_sequences)
    current_sequences <- next_sequences
    
    # Limit total sequences
    if (length(current_sequences) > max_sequences) {
      current_sequences <- current_sequences[1:max_sequences]
    }
    
    # Force garbage collection between steps
    mallinfo::malloc.trim(0L)
  }
  
  # Extract just the sequences (not the dendrograms)
  result <- lapply(current_sequences, function(x) x$sequence)
  rm(current_sequences)
  mallinfo::malloc.trim(0L)
  
  return(result)
}

#' Optimize a dimension of a dendrogram
#' @param dend dendrogram
#' @param labels_to_keep labels to keep
#' @param mat_subset matrix subset
#' @param kernel_size kernel size
#' @param verbose verbose output
#' @param verbose_improvements verbose improvement output
#' @param dimension dimension to optimize
#' @param cache cache object
#' @param consider_mixed_flips if TRUE, include flips that mix family and non-family labels
#' @param consider_minimal_subtree if TRUE, consider flips that affect a minimal subtree containing all family labels
#' @return list with improvement_found, dend, and new_score
try_optimize_dimension_enhanced <- function(dend, labels_to_keep, mat_subset, kernel_size, 
                                            verbose, verbose_improvements, dimension,
                                            subgroup=NULL,
                                            cache = NULL,
                                            consider_mixed_flips = TRUE,
                                            consider_minimal_subtree = FALSE,
                                            consider_subgroups = TRUE,
                                            max_iterations = 50,
                                            max_steps = 4,
                                            max_sequences = 2500,
                                            use_randomization = FALSE) {
  
  calc_score <- function(row_order, col_order) {
    score_by_order(
      mat_subset = mat_subset,
      row_order = row_order,
      col_order = col_order,
      kernel_size = kernel_size,
      cache = cache,
      subgroup = subgroup
    )
  }
  
  # Initialize
  current_dend <- dend
  if (dimension == "rows") {
    current_score <- calc_score(labels(current_dend), colnames(mat_subset))
  } else {
    current_score <- calc_score(rownames(mat_subset), labels(current_dend))
  }
  initial_score <- current_score
  iteration <- 0
  total_improvement_found <- FALSE
  
  if (verbose) cat("    Starting", dimension, "optimization. Initial score:", 
                   sprintf("%.2f", current_score), "\n")
  
  # Keep optimizing until no improvement found
  while (iteration < max_iterations) {
    iteration <- iteration + 1
    iteration_start_score <- current_score
    improvement_found <- FALSE
    best_dend <- current_dend
    best_score <- current_score
    
    if (verbose) cat("      ", dimension, "iteration", iteration, 
                     "- Current score:", sprintf("%.2f", current_score), "\n")
    
    # Get potential flips based on current dendrogram
    if (consider_minimal_subtree) {
      possible_flips <- get_potential_flips_minimal_subtree(current_dend, labels_to_keep)
    } else {
      possible_flips <- get_potential_flips(current_dend, labels_to_keep, consider_mixed_flips)
    }
    
    # Convert to single-step sequences
    possible_one_move_flips <- lapply(possible_flips, function(f) list(f))
    
    # Try all 1-step flips
    if (length(possible_one_move_flips) > 0) {
      
      if (use_randomization) {
        flip_indices <- sample(seq_along(possible_one_move_flips))
      } else {
        flip_indices <- seq_along(possible_one_move_flips)
      }
      
      for (idx in seq_along(flip_indices)) {
        i <- flip_indices[idx]
        flip <- possible_one_move_flips[[i]]
        flipped_dend <- current_dend
        description <- ""
        for (j in seq_along(flip)) {
          flipped_dend <- flip_node_fast(flipped_dend, flip[[j]]$node_id)
          description <- if(j==1) flip[[j]]$description else paste0(description,";",flip[[j]]$description)
        }
        
        new_order <- labels(flipped_dend)
        
        if (dimension == "rows") {
          new_score <- calc_score(new_order, colnames(mat_subset))
        } else {
          new_score <- calc_score(rownames(mat_subset), new_order)
        }
        
        if (verbose && i <= 10) {  # Limit verbose output
          cat("        Flip", i, "(", description, "):", 
              sprintf("%.2f", new_score))
          if (new_score > best_score) cat(" *")
          cat("\n")
        }
        
        if (new_score > best_score) {
          best_score <- new_score
          best_dend <- flipped_dend
          improvement_found <- TRUE
        }
      }
      
      if (improvement_found) {
        # Apply the best flip and continue
        current_dend <- best_dend
        previous_score <- current_score
        current_score <- best_score
        total_improvement_found <- TRUE
        
        improvement <- current_score - previous_score
        improvement_pct <- 100 * improvement / previous_score
        
        if (verbose || verbose_improvements) {
          cat("        Found improvement with 1-step flip:", 
              sprintf("%.2f -> %.2f (+%.2f, %.2f%%)\n", 
                      previous_score, current_score, improvement, improvement_pct))
        }
        
        # Update mat_subset for next iteration based on new ordering
        new_order <- labels(current_dend)
        if (dimension == "rows") {
          mat_subset <- mat_subset[new_order, , drop = FALSE]
        } else {
          mat_subset <- mat_subset[, new_order, drop = FALSE]
        }
        
        rm(possible_one_move_flips, possible_flips)
        mallinfo::malloc.trim(0L)
        # Continue to next iteration
        next
      }
    }
    
    # If no 1-step improvement, try multi-step flips
    if (!improvement_found) {
      if (verbose) cat("        No 1-step improvement, trying multi-step flips...\n")
      
      for (n_steps in 2:max_steps) {
        if (verbose || verbose_improvements) cat("        Trying", n_steps, "-step flips...\n")
        
        all_nstep_sequences <- generate_promising_sequences(
          current_dend, labels_to_keep, n_steps, mat_subset, kernel_size, cache, max_sequences,
          consider_mixed_flips, consider_minimal_subtree, consider_subgroups
        )        
        
        if (length(all_nstep_sequences) > 0) {
          best_seq_score <- best_score
          best_seq_dend <- best_dend
          
          if (use_randomization) {
            sequence_indices <- sample(seq_along(all_nstep_sequences))
          } else {
            sequence_indices <- seq_along(all_nstep_sequences)
          }
          
          for (seq_idx_pos in seq_along(sequence_indices)) {
            seq_idx <- sequence_indices[seq_idx_pos]
            sequence <- all_nstep_sequences[[seq_idx]]
            
            test_dend <- current_dend
            description <- ""
            for (step_idx in seq_along(sequence)) {
              flip <- sequence[[step_idx]]
              test_dend <- flip_node_fast(test_dend, flip$node_id)
              description <- paste0(description, " -> ", flip$node_id)
            }
            
            new_order <- labels(test_dend)
            
            if (dimension == "rows") {
              new_score <- calc_score(new_order, colnames(mat_subset))
            } else {
              new_score <- calc_score(rownames(mat_subset), new_order)
            }
            
            if (verbose) {
              cat("          Sequence", seq_idx, "score:", sprintf("%.2f", new_score))
              cat(" ", description)
              if (new_score > best_seq_score) cat(" *")
              cat("\n")
            }
            
            if (new_score > best_seq_score) {
              best_seq_score <- new_score
              best_seq_dend <- test_dend
              improvement_found <- TRUE
            }
          }
          
          
          if (improvement_found) {
            # Apply the best multi-step sequence
            current_dend <- best_seq_dend
            previous_score <- current_score
            current_score <- best_seq_score
            total_improvement_found <- TRUE
            
            improvement <- current_score - previous_score
            improvement_pct <- 100 * improvement / previous_score
            
            if (verbose || verbose_improvements) {
              cat("        Found improvement with", n_steps, "-step sequence:", 
                  sprintf("%.2f -> %.2f (+%.2f, %.2f%%)\n", 
                          previous_score, current_score, improvement, improvement_pct))
            }
            
            # Update mat_subset for next iteration
            new_order <- labels(current_dend)
            if (dimension == "rows") {
              mat_subset <- mat_subset[new_order, , drop = FALSE]
            } else {
              mat_subset <- mat_subset[, new_order, drop = FALSE]
            }
            
            rm(all_nstep_sequences, test_dend)
            mallinfo::malloc.trim(0L)
            
            break  # Exit n_steps loop since we found improvement
          }
        } else {
          if (verbose) cat("        No valid", n_steps, "-step sequences found\n")
          break
        }
      }
    }
    
    # If still no improvement found after trying multi-step, we're done
    if (!improvement_found) {
      if (verbose) cat("      No further", dimension, "improvements found. Stopping.\n")
      break
    }
  }
  
  # Final summary for this dimension
  if (total_improvement_found) {
    total_improvement <- current_score - initial_score
    total_improvement_pct <- 100 * total_improvement / initial_score
    
    if (verbose || verbose_improvements) {
      cat("    ", dimension, "optimization complete after", iteration, "iterations.\n")
      cat("    Total", dimension, "improvement:", 
          sprintf("%.2f -> %.2f (+%.2f, %.2f%%)\n", 
                  initial_score, current_score, total_improvement, total_improvement_pct))
    }
  } else {
    if (verbose) cat("    No", dimension, "improvements found.\n")
  }
  
  rm(possible_flips)
  if (exists("mat_subset")) rm(mat_subset) 
  mallinfo::malloc.trim(0L)
  
  if (exists("all_nstep_sequences")) {
    rm(all_nstep_sequences)
  }
  if (exists("test_dend")) {
    rm(test_dend)  
  }
  
  return(list(
    improvement_found = total_improvement_found, 
    dend = current_dend, 
    new_score = current_score,
    iterations = iteration
  ))
}

# Add this function before optimize_order_enhanced
#' Final cleanup to separate subgroups as much as possible
#' @param dend dendrogram to clean up
#' @param family_groups list of family groups with their labels
#' @param dimension "rows" or "columns" 
#' @param verbose verbose output
#' @return cleaned dendrogram
final_subgroup_separation <- function(dend, family_groups, dimension = "rows", verbose = FALSE) {
  if (verbose) cat("  Final", dimension, "subgroup separation phase...\n")
  
  current_dend <- dend
  total_separations <- 0
  max_iterations <- 100
  iteration <- 0
  
  # Global optimization loop
  repeat {
    iteration <- iteration + 1
    if (iteration > max_iterations) {
      if (verbose) cat("    Max iterations reached\n")
      break
    }
    
    # Calculate current fragmentation score
    current_order <- labels(current_dend)
    current_fragmentation <- calculate_fragmentation_score(current_order, family_groups)
    
    if (verbose && iteration == 1) {
      cat("    Initial fragmentation score:", current_fragmentation, "\n")
    }
    
    # Find the best flip
    node_index <- build_node_index(current_dend)
    internal_nodes <- get_internal_node_ids_fast(current_dend)
    best_flip <- NULL
    best_improvement <- 0
    best_new_fragmentation <- current_fragmentation
    
    for (node_id in internal_nodes) {
      node <- get_node_by_id_fast(node_index, node_id)
      if (!is.null(node) && length(node) == 2) {
        
        # Test the flip
        test_dend <- flip_node_fast(current_dend, node_id)
        test_order <- labels(test_dend)
        test_fragmentation <- calculate_fragmentation_score(test_order, family_groups)
        
        # Lower fragmentation = better
        improvement <- current_fragmentation - test_fragmentation
        
        if (improvement > best_improvement) {
          best_improvement <- improvement
          best_flip <- node_id
          best_new_fragmentation <- test_fragmentation
        }
      }
    }
    
    # Apply best flip if found
    if (!is.null(best_flip) && best_improvement > 0) {
      current_dend <- flip_node_fast(current_dend, best_flip)
      total_separations <- total_separations + 1
      
      if (verbose) {
        cat("    Iteration", iteration, "- Applied flip\n")
        cat("      Fragmentation:", current_fragmentation, "->", best_new_fragmentation,
            "(improvement:", sprintf("%.1f", best_improvement), ")\n")
      }
    } else {
      # No improvement found
      if (verbose) {
        cat("    No more improvements found after", total_separations, "flips\n")
        cat("    Final fragmentation score:", current_fragmentation, "\n")
      }
      break
    }
  }
  
  # Show final family statistics
  if (verbose) {
    cat("    Final family contiguity:\n")
    final_order <- labels(current_dend)
    for (fam_name in names(family_groups)) {
      fam_labels <- family_groups[[fam_name]]
      longest_block <- longest_contiguous_block(final_order, fam_labels)
      total_members <- sum(fam_labels %in% final_order)
      contiguity_pct <- 100 * longest_block / total_members
      cat("      ", fam_name, ": longest block =", longest_block, "/", total_members,
          sprintf("(%.1f%%)\n", contiguity_pct))
    }
  }
  
  return(current_dend)
}

#' Count the number of "interruptions" in family sequence
#' @param order current order of labels
#' @param family_labels labels belonging to this family
#' @return interruption score
calculate_family_interruptions <- function(order, family_labels) {
  positions <- match(family_labels, order)
  positions <- positions[!is.na(positions)]
  
  if (length(positions) <= 1) return(0)
  
  positions_sorted <- sort(positions)
  
  # Count how many non-consecutive jumps there are
  interruptions <- sum(diff(positions_sorted) > 1)
  
  # Also penalize the size of gaps
  gap_sizes <- diff(positions_sorted) - 1
  gap_penalty <- sum(gap_sizes^2)  # Square penalty for large gaps
  
  return(interruptions + 0.1 * gap_penalty)
}

#' Find the longest contiguous block within a family
#' @param order current order of labels  
#' @param family_labels labels belonging to ONE family
#' @return length of longest contiguous block
longest_contiguous_block <- function(order, family_labels) {
  positions <- match(family_labels, order)
  present_positions <- sort(positions[!is.na(positions)])
  
  if (length(present_positions) <= 1) return(length(present_positions))
  
  # Find longest sequence of consecutive positions
  max_block <- 1
  current_block <- 1
  
  for (i in 2:length(present_positions)) {
    if (present_positions[i] == present_positions[i-1] + 1) {
      current_block <- current_block + 1
    } else {
      max_block <- max(max_block, current_block)
      current_block <- 1
    }
  }
  max_block <- max(max_block, current_block)
  
  return(max_block)
}

#' Simple fragmentation score: prefer longer contiguous blocks
#' @param order current order of labels
#' @param family_groups list of family groups with their labels
#' @return fragmentation score
calculate_fragmentation_score <- function(order, family_groups) {
  total_score <- 0
  
  for (fam_name in names(family_groups)) {
    fam_labels <- family_groups[[fam_name]]
    longest_block <- longest_contiguous_block(order, fam_labels)
    total_members <- sum(fam_labels %in% order)
    
    fragmentation_penalty <- total_members - longest_block
    total_score <- total_score + fragmentation_penalty
  }
  
  return(total_score) 
}



#' Optimize the order of a matrix and dendrograms
#' @param mat matrix to optimize
#' @param row_dend row dendrogram
#' @param col_dend column dendrogram
#' @param kernel_size kernel size
#' @param max_iter_family maximum number of iterations per family
#' @param consider_mixed_flips if TRUE, include flips that mix family and non-family labels
#' @param consider_minimal_subtree if TRUE, consider flips that affect a minimal subtree containing all family labels
#' @param consider_subgroups if TRUE, consider flips that affect subgroups
#' @param use_randomization if TRUE, use randomization to select flips
optimize_order_enhanced <- function(mat, row_dend, col_dend,
                                    kernel_size = 5,
                                    max_iter_family = 50,
                                    max_iter_global = 200,
                                    convergence_threshold = 0.001,
                                    max_steps = 4,
                                    verbose = FALSE,
                                    verbose_improvements = TRUE,
                                    use_cache = TRUE,
                                    consider_mixed_flips = TRUE,
                                    consider_minimal_subtree = FALSE,
                                    consider_subgroups = TRUE,
                                    use_randomization = FALSE,
                                    optimize_subgroups = FALSE) {
  
  initial_score <- calc_blockiness_optimized(mat[
    labels(row_dend),
    labels(col_dend)
  ], kernel_size)
  
  if (verbose) {
    cat("=== Enhanced Two-Stage Blockiness Maximization ===\n")
    cat("Initial blockiness:", initial_score, "\n")
    cat("Settings: cache =", use_cache, ", randomization =", use_randomization, "\n\n")
  }
  
  # Initialize cache if requested
  blockiness_cache <- if (use_cache) create_blockiness_cache() else NULL
  
  # Working copies of dendrograms with IDs
  working_row_dend <- add_node_ids(row_dend)
  working_col_dend <- add_node_ids(col_dend)
  
  # Get family groups
  col_labels <- labels(working_col_dend)
  family_groups_cols <- split(col_labels, get_family(col_labels))
  
  # Identify rows for each family
  family_groups_rows <- lapply(names(family_groups_cols), function(fam) {
    cols <- family_groups_cols[[fam]]
    if (length(cols) > 1) {
      rows <- which(rowSums(mat[, cols, drop = FALSE]) > 0)
    } else {
      rows <- which(mat[, cols] > 0)
    }
    return(rownames(mat)[rows])
  })
  names(family_groups_rows) <- names(family_groups_cols)
  
  # Track global progress
  global_improvements <- numeric()
  
  # Optimize each family
  for (fam in names(family_groups_cols)) {
    if (verbose_improvements) cat("Optimizing family", fam, "\n")
    
    fam_cols <- family_groups_cols[[fam]]
    fam_rows <- family_groups_rows[[fam]]
    
    # Skip if family is too small
    if (length(fam_cols) <= 1 || length(fam_rows) <= 1) {
      if (verbose) cat("  Skipping family - too small\n")
      next
    }
    
    mat_subset <- mat[labels(working_row_dend), labels(working_col_dend)]
    mat_subset[!labels(working_row_dend) %in% fam_rows, !labels(working_col_dend) %in% fam_cols] <- 0
    
    family_initial_score <- calc_blockiness_optimized(mat_subset, kernel_size)
    
    if (verbose) cat("  Initial family score:", family_initial_score, "\n")
    
    # Optimization loop for this family
    iteration <- 0
    current_score <- family_initial_score
    improvement_found <- TRUE
    consecutive_small_improvements <- 0
    score_history <- numeric()
    
    while (improvement_found && iteration < max_iter_family) {
      improvement_found <- FALSE
      iteration <- iteration + 1
      previous_score <- current_score
      
      if (verbose) cat("  Iteration", iteration, "- Current score:", current_score, "\n")
      
      # Try to improve rows
      if (verbose) cat("    Trying row optimization...\n")
      
      # Get current order and subset
      current_row_order <- labels(working_row_dend)
      current_col_order <- labels(working_col_dend)
      
      # instead of subsetting we will change the rest to zero
      mat_subset_current <- mat[current_row_order, current_col_order]
      mat_subset_current[!current_row_order %in% fam_rows, !current_col_order %in% fam_cols] <- 0
      
      row_improvement <- try_optimize_dimension_enhanced(
        dend = working_row_dend,
        labels_to_keep = fam_rows,
        mat_subset = mat_subset_current,
        kernel_size = kernel_size,
        verbose = verbose,
        verbose_improvements = verbose_improvements,
        dimension = "rows",
        cache = blockiness_cache,
        consider_mixed_flips = consider_mixed_flips,
        consider_minimal_subtree = consider_minimal_subtree,
        max_iterations = max_iter_family,
        max_steps = max_steps,
        max_sequences = if (nrow(mat_subset_current) > 150L || ncol(mat_subset_current) > 80L) 300L else 1200L,
        consider_subgroups = consider_subgroups,
        use_randomization = use_randomization,
        subgroup = fam
      )
      
      if (row_improvement$improvement_found) {
        working_row_dend <- row_improvement$dend
        current_score <- row_improvement$new_score
        improvement_found <- TRUE
        if (verbose_improvements) {
          cat("    Row improvement found! New score:", current_score, "\n")
        }
      }
      
      # Try to improve columns
      if (verbose) cat("    Trying column optimization...\n")
      
      # Update matrix subset with new row order
      current_row_order <- labels(working_row_dend)
      current_col_order <- labels(working_col_dend)
      mat_subset_updated <- mat[current_row_order, current_col_order]
      mat_subset_updated[!current_row_order %in% fam_rows, !current_col_order %in% fam_cols] <- 0
      
      col_improvement <- try_optimize_dimension_enhanced(
        dend = working_col_dend,
        labels_to_keep = fam_cols,
        mat_subset = mat_subset_updated,
        kernel_size = kernel_size,
        verbose = verbose,
        verbose_improvements = verbose_improvements,
        dimension = "columns",
        cache = blockiness_cache,
        consider_mixed_flips = consider_mixed_flips,
        consider_minimal_subtree = consider_minimal_subtree,
        consider_subgroups = consider_subgroups,
        max_sequences = if (nrow(mat_subset_updated) > 150L || ncol(mat_subset_updated) > 80L) 300L else 1200L,
        use_randomization = use_randomization,
        subgroup = fam
      )
      
      if (col_improvement$improvement_found) {
        working_col_dend <- col_improvement$dend
        current_score <- col_improvement$new_score
        improvement_found <- TRUE
        if (verbose_improvements) {
          cat("    Column improvement found! New score:", current_score, "\n")
        }
      }
      
      # Check for convergence
      if (improvement_found) {
        improvement_ratio <- abs(current_score - previous_score) / max(previous_score, 1)
        score_history <- c(score_history, current_score)
        
        if (improvement_ratio < convergence_threshold) {
          consecutive_small_improvements <- consecutive_small_improvements + 1
          if (consecutive_small_improvements >= 3) {
            if (verbose) {
              cat("  Converged after", iteration, "iterations",
                  "(improvement ratio:", round(improvement_ratio, 6), ")\n")
            }
            break
          }
        } else {
          consecutive_small_improvements <- 0
        }
      } else {
        if (verbose) cat("    No improvement found in either dimension\n")
      }
      
      # Clear cache periodically to manage memory
      if (use_cache && !is.null(blockiness_cache) && blockiness_cache$size() > 5000) {
        if (verbose) cat("  Clearing cache (size:", blockiness_cache$size(), ")\n")
        blockiness_cache$clear()
        mallinfo::malloc.trim(0L)
      }
    }
    
    # Report family results
    improvement <- current_score - family_initial_score
    improvement_pct <- 100 * improvement / family_initial_score
    
    cat(sprintf("  Final family (%s) score: %.1f (+%.1f, %.1f%% improvement) after %d iterations\n", 
                fam, current_score, improvement, improvement_pct, iteration))
    
    global_improvements <- c(global_improvements, improvement_pct)
  }
  
  # Final reporting
  if (verbose_improvements) {
    final_score <- calc_blockiness_optimized(mat[
      labels(working_row_dend),
      labels(working_col_dend)
    ], kernel_size)
    
    cat("\n=== Optimization Complete ===\n")
    cat("Final blockiness:", final_score, "\n")
    cat("Overall improvement:", 
        sprintf("%.1f%%", 100 * (final_score - initial_score) / initial_score), "\n")
    cat("Average family improvement:", 
        sprintf("%.1f%%", mean(global_improvements)), "\n")
    
    if (use_cache && !is.null(blockiness_cache)) {
      cat("Cache statistics: final size =", blockiness_cache$size(), "\n")
    }
  }
  
  if(optimize_subgroups){
    # === FINAL SUBGROUP SEPARATION PHASE ===
    if (verbose_improvements) cat("\n=== Final Subgroup Separation Phase ===\n")
    
    # Apply final separation for rows
    if (verbose_improvements) cat("Cleaning up row subgroup mixing...\n")
    working_row_dend <- final_subgroup_separation(
      dend = working_row_dend,
      family_groups = family_groups_rows,
      dimension = "rows",
      verbose = verbose_improvements
    )
    
    # Apply final separation for columns  
    if (verbose_improvements) cat("Cleaning up column subgroup mixing...\n")
    working_col_dend <- final_subgroup_separation(
      dend = working_col_dend, 
      family_groups = family_groups_cols,
      dimension = "columns", 
      verbose = verbose_improvements
    )
    
    # Calculate final score after cleanup
    final_cleanup_score <- calc_blockiness_optimized(mat[
      labels(working_row_dend),
      labels(working_col_dend)
    ], kernel_size)
    
    if (verbose_improvements || verbose) {
      cat("\n=== Optimization Complete ===\n")
      cat("Final blockiness:", final_cleanup_score, "\n")
      cat("Score after subgroup separation:", sprintf("%.2f", final_cleanup_score), "\n")
    }
  }
  
  # Clear cache before returning
  if (use_cache && !is.null(blockiness_cache)) {
    blockiness_cache$clear()
    blockiness_cache$finalize()
    blockiness_cache <- NULL
    rm(blockiness_cache)
    mallinfo::malloc.trim(0L)
  }
  
  rm(mat_subset, mat_subset_current, mat_subset_updated)
  # Return the optimized dendrograms
  mallinfo::malloc.trim(0L)
  return(list(
    row_dend = working_row_dend, 
    col_dend = working_col_dend,
    improvements = global_improvements
  ))
}


######################################

collapse_clade <- function(tree_, node) {
  ape::drop.tip(tree_, tip = node[-1])
}

# Union–find over iuis_group edges
uf <- function(edges) {
  nodes <- unique(c(edges$iuis_group.x, edges$iuis_group.y))
  parent <- setNames(as.list(nodes), nodes)
  find <- function(x) { while(parent[[x]] != x) x <- parent[[x]]; x }
  union <- function(x, y) { rx <- find(x); ry <- find(y); if (rx != ry) parent[[ry]] <<- rx }
  for (i in seq_len(nrow(edges))) union(edges$iuis_group.x[i], edges$iuis_group.y[i])
  comp <- sapply(nodes, find)
  data.table(iuis_group = names(comp), cluster = comp)
}

# Build allele→merged-group map for a given locus (e.g., "IGHV" or "IGHD")
build_allele_group_map <- function(husa_path, locus, keep_labels, allow.cartesian=FALSE) {
  hs <- fread(husa_path)
  if (!"iuis_allele" %in% names(hs)) {
    if ("husa" %in% names(hs)) {
      hs[, iuis_allele := husa]
    } else {
      hs[, iuis_allele := NA_character_]
    }
  }
  hs <- hs[gene_type == locus, .(gene_type, allele, iuis_allele)]
  hs <- hs[!is.na(iuis_allele) & iuis_allele != ""]
  hs <- hs[!duplicated(paste(gene_type, allele, iuis_allele))]
  hs <- hs[, .(iuis_allele_dup = unlist(strsplit(iuis_allele, ","))),
           by = .(gene_type, allele, iuis_allele)]
  hs <- hs[!is.na(iuis_allele_dup) & iuis_allele_dup != ""]
  hs[, iuis_group := alakazam::getGene(iuis_allele_dup, strip_d = FALSE, omit_nl = FALSE)]
  hs <- hs[iuis_group %in% paste0(substr(locus, 1, 3), keep_labels)]
  if (nrow(hs) == 0L) {
    return(setNames(character(), character()))
  }
  
  links <- merge(hs, hs, by = "iuis_allele", allow.cartesian=allow.cartesian)[iuis_group.x != iuis_group.y,
                                             .(iuis_group.x, iuis_group.y)]
  if (nrow(links) == 0L) {
    # No overlaps, trivial mapping
    out <- hs[, .(allele, iuis_group)]
    setnames(out, "iuis_group", "all_iuis_groups")
    return(setNames(out$all_iuis_groups, out$allele))
  }
  
  clusters <- uf(links)
  clusters[, all_iuis_groups := paste(sort(unique(iuis_group)), collapse = ","), by = cluster]
  result <- merge(hs, clusters[, .(iuis_group, all_iuis_groups)], by = "iuis_group")
  result <- result[!duplicated(paste(allele, all_iuis_groups))]
  
  alleles_non_dup <- hs[!allele %in% unique(result$allele)]
  alleles_non_dup <- alleles_non_dup[!duplicated(allele, iuis_group)]
  
  setNames(
    c(result$all_iuis_groups, alleles_non_dup$iuis_group),
    c(result$allele,         alleles_non_dup$allele)
  )
}

# Build collapsed usage objects for a locus from repertoire table and allele-group map
build_usage_objects <- function(repertoire_csv, locus, call_col, group_map, keep_labels) {
  ggrep <- fread(repertoire_csv)
  subj_col <- "vdjbase_subject"
  
  # Map calls to merged iuis groups
  ggrep[, merged := group_map[get(call_col)]]
  
  usage <- ggrep[!is.na(merged),
                 .(count = .N),
                 by = .(subject = get(subj_col), merged)]
  usage[, total := sum(count), by = subject]
  usage[, rel_usage := count / total]
  
  # Long explode: iuis_group from merged
  usage_dup <- usage[, .(
    iuis_group = gsub("IG[KLH]", "", unlist(strsplit(merged, ",")))
  ), by = .(subject, merged, count, total, rel_usage)]
  usage_dup <- usage_dup[iuis_group %in% keep_labels]
  usage_dup[, iuis_group := factor(iuis_group, levels = keep_labels)]
  setorder(usage_dup, subject, iuis_group)
  collapsed <- unique(usage_dup[, .(subject, merged, rel_usage)])
  collapsed[,merged:=gsub("IG[KLH]","",merged)]
  # Matrix for collapsed boxplots: one column per merged v_iuis_group
  m_box <- dcast.data.table(unique(collapsed[, .(subject, merged, rel_usage)]),
                            subject ~ merged, value.var = "rel_usage", fill = 0)
  m_box[, subject := NULL]
  m_box <- as.matrix(m_box)
  
  # Map iuis_group → merged, aligned to heatmap column order later
  map <- unique(usage_dup[, .(iuis_group, merged)])
  list(usage_dup = usage_dup, m_box = m_box, map = map)
}

# Build an anno_link that draws exactly one boxplot per merged group, centered over its columns
# align_to_vec must be in the same order as the heatmap columns
build_collapsed_box_anno <- function(m_box, align_to_vec, ylab = NULL, show_axis_at = 1L,
                                     axis_fontsize = 12, ylab_fontsize = 12, ylab_x = -2.2,
                                     ylab_x_npc = TRUE, ylab_y = 0.5) {
  rg <- range(m_box)
  HeatmapAnnotation(
    usage = anno_link(
      align_to = align_to_vec,
      which    = "column",
      height   = unit(4, "cm"),
      link_height = unit(0.1, "cm"),
      link_gp  = gpar(col = "white", lwd = 0),
      panel_fun = local({
        atv <- align_to_vec
        function(index, nm) {
          # Determine the merged-group name for this block
          grp <- unique(atv[index])
          vals <- as.numeric(m_box[, grp, drop = TRUE])
          pushViewport(viewport(xscale = c(0, 1), yscale = rg))
          if (!is.null(ylab) && is.numeric(show_axis_at) && length(index) > 0) {
            if (min(index) == show_axis_at) {
              grid.yaxis(gp = gpar(fontsize = axis_fontsize))
              grid.text(ylab, x = if (ylab_x_npc) unit(ylab_x, "npc") else unit(ylab_x, "mm"),
                        y = unit(ylab_y, "npc"), rot = 90, just = "centre",
                        gp = gpar(fontsize = ylab_fontsize))
            }
          }
          grid.boxplot(vals, pos = 0.5, direction = "vertical")
          popViewport()
        }
      })
    )
  )
}

create_divergent_colors <- function(breaks, scheme = "blue_red") {
  schemes <- list(
    # Current scheme (orange to red)
    "current" = c("white", "#F39300", "#7D0025", "#330000"),
    
    # Blue to Red (classic divergent)
    "blue_red" = c("#2166AC", "#92C5DE", "#F4A582", "#D6604D"),
    
    # Purple to Green (colorblind friendly)
    "purple_green" = c("#762A83", "#C2A5CF", "#A6DBA0", "#008837"),
    
    # Blue to Yellow to Red (more gradual)
    "blue_yellow_red" = c("#4575B4", "#74ADD1", "#FEE090", "#F46D43"),
    
    # Green to Yellow to Red (traffic light)
    "traffic_light" = c("#1A9850", "#91CF60", "#FEE08B", "#FC8D59"),
    
    # Purple to Orange (high contrast)
    "purple_orange" = c("#5E4FA2", "#9E0142", "#F46D43", "#FDAE61"),
    
    # Blue to White to Red (centered divergent)
    "blue_white_red" = c("#053061", "#4393C3", "#FDDBC7", "#67001F"),
    
    # Viridis-like (modern)
    "viridis_like" = c("#440154", "#31688E", "#35B779", "#FDE725"),
    
    # Plasma-like (warm)
    "plasma_like" = c("#0D0887", "#7E03A8", "#CC4778", "#F89441"),
    
    # Cool to Warm (scientific)
    "cool_warm" = c("#3B4CC0", "#6B93D3", "#F4A582", "#D6604D")
  )
  
  colors <- schemes[[scheme]]
  
  # If breaks start with zero, use white as the first color
  if (breaks[1] == 0) {
    if (length(breaks) == 3) {
      colors <- c("white", colors[c(3, 4)]) # white + last 2 colors
    } else {
      colors <- c("white", colors[1:(length(breaks) - 1)]) # white + all other colors
    }
  } else if (length(breaks) == 3) {
    colors <- colors[c(1, 3, 4)] # Use 3 colors for 3 breaks (original logic)
  }
  
  return(circlize::colorRamp2(breaks, colors))
}

generate_consensus <- function(rss_data_segment, vec_onehot) {
  consensus <- consensusString(DNAStringSet(rss_data_segment), ambiguityMap="?")
  consensus_onehot <- sapply(strsplit(consensus, "")[[1]], function(x) {
    vec_onehot[[x]]
  })
  list(consensus = consensus, consensus_onehot = consensus_onehot)
}

create_row_labels <- function(row_order, consensus, nucleotide_colors, nonamer_start, rss=T, fontsize_px = 26) {
  consensus_split <- strsplit(consensus, "")[[1]]
  sapply(
    row_order, function(s) {
      s <- strsplit(s, "")[[1]]
      paste0(sapply(seq_along(s), function(i) {
        nucleotide <- s[i]
        if (nucleotide == consensus_split[i]) {
          nuc_cons <- '.'
          text <- "<span style='color: @{nucleotide_colors[nucleotide]}; font-family: mono; font-size:@{fontsize_px}px;'>@{nuc_cons}</span>"
        } else {
          text <- "<span style='color: @{nucleotide_colors[nucleotide]}; font-family: mono; font-size:@{fontsize_px}px;'>@{nucleotide}</span>"
        }
        if (i == 8 & rss==TRUE) {
          text <- paste0("<span style='color: @{nucleotide_colors[nucleotide]}; font-family: mono; font-size:@{fontsize_px}px;'> </span>", text)
        }
        if (i == nonamer_start) {
          text <- paste0(text, "<span style='color: @{nucleotide_colors[nucleotide]}; font-family: mono; font-size:@{fontsize_px}px;'> </span>")
        }
        GetoptLong::qq(text)
      }), collapse = "")
    }
  )
}

create_consensus_label <- function(consensus, nucleotide_colors, nonamer_start, rss=T, fontsize_px = 26, box_positions = integer(0)) {
  consensus_split <- strsplit(consensus, "")[[1]]
  paste0(sapply(seq_along(consensus_split), function(i) {
    x <- consensus_split[i]
    # gridtext spans do not support background/border, so changed positions are marked
    # with bold red glyphs (the robust equivalent of "boxing" the changed positions).
    box_css <- if (i %in% box_positions) " color:#d40000; font-weight:bold;" else ""
    text <- "<span style='color: @{nucleotide_colors[x]}; font-family: mono; font-size:@{fontsize_px}px;@{box_css}'>@{x}</span>"
    if (i == 8 && rss==TRUE) {
      text <- paste0("<span style='color: @{nucleotide_colors[x]}; font-family: mono; font-size:@{fontsize_px}px;'> </span>", text)
    }
    if (i == nonamer_start) {
      text <- paste0(text, "<span style='color: @{nucleotide_colors[x]}; font-family: mono; font-size:@{fontsize_px}px;'> </span>")
    }
    GetoptLong::qq(text)
  }), collapse = "")
}


subgroups_annotation_custom <- function(val, colors_rss_subgroup){
  families <- strsplit(val, ",")[[1]]
  families <- trimws(families)
  
  if (length(families) == 1) {
    # Single family: draw square
    function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = colors_rss_subgroup[families], col = NA))
    }
  } else if (length(families) == 2) {
    # Two families: draw two triangles
    function(x, y, w, h) {
      grid.polygon(
        x = unit.c(x - 0.5*w, x + 0.5*w, x + 0.5*w),
        y = unit.c(y + 0.5*h, y + 0.5*h, y - 0.5*h),
        gp = gpar(fill = colors_rss_subgroup[families[1]], col = NA)
      )
      grid.polygon(
        x = unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
        y = unit.c(y + 0.5*h, y - 0.5*h, y - 0.5*h),
        gp = gpar(fill = colors_rss_subgroup[families[2]], col = NA)
      )
    }
  } else {
    # fallback if >2 (optional)
    function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "black", col = NA))
    }
  }
}


colors15 <- setNames(
  c(
  "#E41A1C", # red
  "#377EB8", # blue
  "#4DAF4A", # green
  "#984EA3", # purple
  "#FF7F00", # orange
  "#FFFF33", # yellow
  "#A65628", # brown
  "#F781BF", # pink
  "#999999", # gray
  "#66C2A5", # teal
  "#FC8D62", # salmon
  "#8DA0CB", # light blue
  "#E78AC3", # light pink
  "#A6D854", # lime green
  "#FFD92F"  # gold
), paste0("V", 1:15))

### create the rss heatmap ordered by the ch gene order
vec_onehot <- list(
  "A" = c(1,0,0,0,0),
  "T" = c(0,1,0,0,0),
  "C" = c(0,0,1,0,0),
  "G" = c(0,0,0,1,0),
  "-" = c(0,0,0,0,1),
  "?" = c(0,0,0,0,0)
)

# Color scheme and nucleotide colors
COLOR_SCHEME <- "blue_red"
col_fun <- create_divergent_colors(c(0, 3, 5, 10, 15), COLOR_SCHEME)
cl <- c("seagreen", "darkblue", "darkorange3", "firebrick4", "firebrick4","black","black","darkmagenta")
nucleotide_colors <- setNames(cl, c("A", "C", "G", "T", "U",".","?","-"))

########################################
