################################################################################
# Concentration parameter has little bearing
# 3_1 still inpreted as 2_1
# ==============================================================================
rm(list = ls())
par_default <- par(no.readonly = TRUE)
load(sprintf("../data/Stranger.rda"))
attached <- search() # Check no Half_siblings already attached
if(any(grepl("siblings", attached))) stop('detach("Full_siblings")')

attach(Stranger)
n_repeats <- length(ys_store[[1]])
c_params <- names(ys_store)
n_markers <- as.numeric(names(ps_store[[1]][[1]]))
cols <- RColorBrewer::brewer.pal(n = n_repeats, "Paired") # Colours for repeats
cols_light <- sapply(cols, adjustcolor, alpha = 0.25)

# ==============================================================================
# Plot results generated given equifrequent alleles across all marker counts
# ==============================================================================
# Plot the posterior relapse probability trajectories
plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,1), bty = "n", las = 1, xaxt = "n",
     xlab = "Marker count", ylab = "Posterior reinfection probability")
axis(side = 1, at = seq(0, max(n_markers), 10))
for(i in 1:n_repeats){
  lines(x = 1:max(n_markers), y = ps_store_all_ms[[as.character(i)]], col = cols[i], lwd = 2)
}

# ==============================================================================
# Process results results for different allele frequency types and for different
# marker counts
# ==============================================================================

# Extract probability of relapse
post_I <- sapply(ps_store, function(X) {
  sapply(X, function(XX) sapply(XX, function(post) post$marg[,"I"]))
}, simplify = F)

# Plots posterior reinfection probabilities
par(mfrow = c(3,1))
for(c in c_params){
  plot(NULL, xlim = range(n_markers)+c(-10,10), ylim = c(0,1),
       xaxt = "n", bty = "n", panel.first = grid(nx = NA, ny = NULL),
       ylab = "Posterior reinfection probability",
       xlab = "Number of markers",
       main = sprintf("Concentration parameter: %s", c))
  axis(side = 1, at = n_markers)
  for(i in 1:n_repeats) {
    lines(y = post_I[[as.character(c)]][,i], x = n_markers,
          lwd = 1.5, col = cols_light[i], pch = 20)
    points(y = post_I[[as.character(c)]][,i], x = n_markers,
           lwd = 1.5, col = cols[i], pch = 20)
  }
}

# Aside: check graphs are all ordered the same [make into a unit test?]
# Extract graph summary (i.e., discard logp)
justRGs <- sapply(ps_store, function(XX) {
    sapply(XX, function(XXX) {
      sapply(XXX, function(post) {
        sapply(post$RGs, function(RG) c(RG$clone, RG$sib))
      }, simplify = F)
  }, simplify = F)
}, simplify = F)


# Check all the graphs are returned in the same order
justRG <- justRGs[[1]][[1]][[1]]
RGcheck <- sapply(c_params, function(c) {
    sapply(1:n_repeats, function(i) {
      sapply(n_markers, function(m) {
        identical(justRG, justRGs[[as.character(c)]][[i]][[as.character(m)]])
      })
  })
})

if (!all(RGcheck)) stop("graphs not returned in the same order")
example_y <- ys_store[[1]][[1]]
example_MOIs <- determine_MOIs(example_y)
ts_per_gs <- rep(1:length(example_y), example_MOIs)
gs <- paste0("g", 1:sum(example_MOIs))

# Extract probability of the data given the relationship graph
llikeRGs <- sapply(ps_store, function(XX) {
    sapply(XX, function(XXX) {
      sapply(XXX, function(post) sapply(post$RGs, function(RG) RG$logp))
    }, simplify = F)
}, simplify = F)

# Re-order colours to graph type
graph_cols <- RColorBrewer::brewer.pal(n = 9, "Paired")
graph_plot_order <- c(4, 2, 6, 8, 5, 9, 1, 3, 7)
names(graph_cols) <- graph_plot_order

# Plot graphs
par(mfrow = c(3,3))
for(g in graph_plot_order) {
  ps <- ps_store[[1]][[1]][[1]]
  RG <- ps$RGs[[g]]
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  igraphRG <- RG_to_igraph(RG, gs, ts_per_gs) # Convert to igraph object
  igraphRG <- igraph::set_vertex_attr(igraphRG, "name", value = NA) # Remove genotype names
  plot_RG(RG =  igraphRG, vertex_palette = "Greys", labels = NA)
  box(col = graph_cols[as.character(g)], lwd = 3)
}

# For each m, c combination, plot the graph likelihood and data
for(c in 100){ # Just focus on one since concentration parameter has little bearing
    par(mfcol = c(n_repeats,length(c_params)), mar = c(0,0,0,0))
    for(m in n_markers){
      for(i in 1:n_repeats) {
        x <- llikeRGs[[as.character(c)]][[as.character(i)]][,as.character(m)][graph_plot_order]
        x[x == -Inf] <- NA # Mask -Inf
        x <- x - min(x, na.rm = T) # Re-scale before exponentiating (otherwise all 0)
        x <- exp(x)/sum(exp(x), na.rm = T) # Exponentiate and normalise
        pie(x[!is.na(x)], col = graph_cols[!is.na(x)], labels = NA, border = NA)
        symbols(x = 0, y = 0, circles = c(0.25), inches = FALSE,
                bg = "white", fg = cols[i], add = TRUE, lwd = 2)
        text(label = round(post_I[[as.character(c)]][as.character(m),i], 2), x = 0, y = 0, cex = 0.75)
      }
      mtext(text = sprintf("concentration %s, marker count %s", c, m), side = 1, cex = 0.5, line = -1)
  }
}

detach("Stranger")






