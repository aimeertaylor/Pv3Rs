################################################################################
#
# ==============================================================================
rm(list = ls())
par_default <- par(no.readonly = TRUE)
relapsing_parasites <- c("Stranger", "Clone", "Full_sibling", "Meiotic_sibling")

for(relapsing_parasite in relapsing_parasites){

  par(par_default)
  attached <- search()
  if(exists("ys_store")) detach("output")
  load(sprintf("./%s.rda", relapsing_parasite))
  attach(output)
  cols <- RColorBrewer::brewer.pal(n = n_repeats, "Paired") # Colours for repeats
  cols_light <- sapply(cols, adjustcolor, alpha = 0.25)

  # ==============================================================================
  # Plots given equifrequent alleles across all marker counts
  # ==============================================================================
  plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,1), bty = "n", las = 1,
       xaxt = "n", xlab = "Marker count", ylab = "Posterior expected state probability")
  axis(side = 1, at = seq(0, max(n_markers), 10))
  title(main = relapsing_parasite)

  legend("bottom", col = cols, lwd = 3, inset = 0, legend = 1:n_repeats, horiz = TRUE, bty = "n")
  for(MOIs in MOIs_per_infection) {
    LTY = ifelse(MOIs == "2_1", 1, 2)
    for(i in 1:n_repeats){
      lines(x = 1:max(n_markers),
            y = ps_store_all_ms[[MOIs]][[as.character(i)]], col = cols[i], lwd = 2, lty = LTY)
    }
  }

  # ==============================================================================
  # Plots for different allele frequencies and for marker counts in n_markers
  # ==============================================================================
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





