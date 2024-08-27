################################################################################
# Plot results for an initial infection of two meiotic siblings (and three in a
# sibling case) and a MOI=1 case recurrence which is either a stranger, clone,
# regular sibling (i.e., a full sibling drawn independently), or meiotic sibling
# (i.e. a full sibling drawn dependently).
#
# For each case, plot results when alleles are equifrequent and data are
# available on all marker counts, and when alleles are not equifrequent and data
# are on a subset of marker counts.
#
# For equifrequent alleles across all marker counts, plot the trajectory of the
# posterior probability of the expected state (line plot), and of all states
# (simplex plot), the latter showing where the posterior probability goes when
# it evades the expected state.
#
# For different allele frequencies (generated using different Dirichlet
# concentration parameters) and a subset of marker counts, plot the
# trajectory of the posterior probability of the expected state (line plot), and
# the relative likelihood of the relationship graphs.
#
# Bulk prevalence data from three meiotic siblings are identical to bulk data
# from parents (furthermore, bulk relative abundance data from four meiotic
# siblings are identical to bulk relative abundance data from parents). As such,
# infections with three or more meiotic siblings are liable to be classified as
# infections of strangers with a clonal edge to a sibling relapse, whether or
# not the relapsing sibling is a full sibling or a meiotic sibling. Aside:
# single cell data are needed for superinfection vs co-transmission.
################################################################################

# Set working directory to source file location
rm(list = ls())
par_default <- par(no.readonly = TRUE)
cases <- c("Stranger", "Clone", "Regular_sibling", "Meiotic_sibling")

for(case in cases[4]){

  par(par_default)
  attached <- search()
  if(exists("ys_store")) detach("output")
  load(sprintf("./%s.rda", case))
  attach(output)
  cols <- RColorBrewer::brewer.pal(n = n_repeats, "Paired") # Colour repeats
  cols_light <- sapply(cols, adjustcolor, alpha = 0.25)

  # ==============================================================================
  # Plots given equifrequent alleles across all marker counts
  # ==============================================================================
  par(mfrow = c(2,1))

  # PLot trajaectories
  plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,1), bty = "n", las = 1,
       xaxt = "n", xlab = "Marker count", ylab = "Posterior expected state probability")
  axis(side = 1, at = seq(0, max(n_markers), 10))
  title(main = gsub("_", " ", case))
    for(MOIs in MOIs_per_infection) {
    LTY = ifelse(MOIs == "2_1", 1, 2)
    for(i in 1:n_repeats){
      lines(x = 1:max(n_markers),
            y = sapply(ps_store_all_ms[[MOIs]][[as.character(i)]], function(x) x[,exp_state]),
            col = cols[i], lwd = 2, lty = LTY)
    }
  }
  legend("topright", col = cols, lwd = 3, title = "Repeats", box.col = "white",
         legend = 1:n_repeats)
  legend("bottomright", lty = c(1,2), lwd = 2, title = "MOIs", box.col = "white",
         legend = gsub("_", " and ", MOIs_per_infection))

  # Plot simplex (helpful when the posterior evades the expected state)
  par(mar = c(0,0,0,0))
  V_labels <- c("Recrudescence", "Relapse", "Reinfection")
  plot_simplex(v_labels =  V_labels, classifcation_threshold = 0.5)
  for(MOIs in MOIs_per_infection) {
    LTY = ifelse(MOIs == "2_1", 1, 2)
    for(i in n_repeats){
      xy_post <- cbind(c(0,0), apply(do.call(rbind, ps_store_all_ms[[MOIs]][[as.character(i)]]), 1, project2D))
      lines(x = xy_post["x",], y = xy_post["y",], lty = LTY)
      points(x = xy_post["x",], y = xy_post["y",], pch = "-")}
    }

  # ==============================================================================
  # Plots for different allele frequencies and for marker counts in n_markers
  # ==============================================================================
  par(par_default)
  par(mfrow = c(3,1))
  for(c in c_params){
    plot(NULL, xlim = range(n_markers)+c(-10,10), ylim = c(0,1),
         xaxt = "n", bty = "n", panel.first = grid(nx = NA, ny = NULL),
         ylab = "Posterior expected state probability",
         xlab = "Number of markers",
         main = sprintf("Concentration parameter: %s", c))
    axis(side = 1, at = n_markers)
    for(MOIs in MOIs_per_infection) {
      for(i in 1:n_repeats) {
        lines(y = post_S[[as.character(c)]][[MOIs]][,i], x = n_markers,
              lwd = 1.5, col = cols_light[i], pch = 20,
              lty = ifelse(MOIs == "2_1", 1, 2))
        points(y = post_S[[as.character(c)]][[MOIs]][,i], x = n_markers,
               lwd = 1.5, col = cols[i], pch = 20)
      }
    }
  }

  example_y <- ys_store[[1]][[1]][[1]]
  example_MOIs <- determine_MOIs(example_y)
  ts_per_gs <- rep(1:length(example_y), example_MOIs)
  gs <- paste0("g", 1:sum(example_MOIs))

  graph_cols <- RColorBrewer::brewer.pal(n = 9, "Paired")
  graph_plot_order <- c(4, 2, 6, 8, 5, 9, 1, 3, 7)
  names(graph_cols) <- graph_plot_order

  # Plot graphs
  par(mfrow = c(3,3))
  for(g in graph_plot_order) {
    ps <- ps_store[[1]][[1]][[1]][[1]]
    RG <- ps$RGs[[g]]
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    igraphRG <- RG_to_igraph(RG, gs, ts_per_gs) # Convert to igraph object
    igraphRG <- igraph::set_vertex_attr(igraphRG, "name", value = NA) # Remove names
    plot_RG(RG =  igraphRG, vertex_palette = "Greys", labels = NA)
    box(col = graph_cols[as.character(g)], lwd = 3)
  }

  # For each m, c combination, plot the graph likelihood
  for(c in c_params){
    for(MOIs in MOIs_per_infection){
      par(mfcol = c(n_repeats,length(c_params)), mar = c(0,0,0,0))
      for(m in n_markers){
        for(i in 1:n_repeats) {
          x <- llikeRGs[[as.character(c)]][[MOIs]][[as.character(i)]][,as.character(m)][graph_plot_order]
          x[x == -Inf] <- NA # Mask -Inf
          x <- x - min(x, na.rm = T) # Re-scale before exponentiating (otherwise all 0)
          x <- exp(x)/sum(exp(x), na.rm = T) # Exponentiate and normalise
          pie(x[!is.na(x)], col = graph_cols[!is.na(x)], labels = NA, border = NA)
          symbols(x = 0, y = 0, circles = c(0.25), inches = FALSE,
                  bg = "white", fg = cols[i], add = TRUE, lwd = 2, lty = ifelse(MOIs == "2_1", 1, 2))
          text(label = round(post_S[[as.character(c)]][[MOIs]][as.character(m),i], 2), x = 0, y = 0, cex = 0.75)
        }
        mtext(text = sprintf("concentration %s, marker count %s", c, m), side = 1, cex = 0.5, line = -1)
      }
    }
  }
  detach("output")
}





