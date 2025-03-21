################################################################################
# Plot results for an initial infection of two or three meiotic siblings and a
# recurrent stranger, clone, regular sibling, or meiotic sibling.
#
# For each case, plot results when alleles are equifrequent and data are
# available on all marker counts (trajectories and simplex plots), and when
# alleles are not equifrequent and data are on a subset of marker counts (line
# plots and pie charts).
#
# For equifrequent alleles across all marker counts, plot the trajectory of the
# posterior probability of the expected state (line plot), and of all states
# (simplex plot), the latter showing where the posterior probability goes when
# it evades the expected state.
#
# Bulk prevalence data from three or four meiotic siblings are identical to bulk
# data from parents. As such, infections with three meiotic siblings are liable
# to be classified as infections of strangers with a clonal edge to a sibling
# relapse (regular or meiotic). And also show what happens in the two selfed
# followed by one recombinant parent-child-like sibling case.
################################################################################

# Set working directory to source file location
rm(list = ls())
par_default <- par(no.readonly = TRUE)
cases <- c("Stranger", "Clone", "Regular_sibling", "Meiotic_sibling")
load("../../data/output_Cln.Str.Sib.rda")

for(case in cases){

  attached <- search()
  if("output" %in% attached) detach("output")
  if(exists("output")) rm("output")
  output <- output_Cln.Str.Sib[[case]]
  attach(output)

  par(par_default)
  cols <- RColorBrewer::brewer.pal(n = n_repeats, "Paired") # Colour repeats

  # Compute effective cardinality cumulatively in m_rorder
  cum_card_eff <- sapply(fs_store, function(fs) {
    cumsum(sapply(fs[m_rorder], function(x) 1/sum(x^2)))})

  # Determine which fs in cum_card_eff are equifrequent
  equifs <- which(as.numeric(colnames(cum_card_eff)) > c_cutoff)

  # ==============================================================================
  # Plots given equifrequent alleles across all marker counts
  # ==============================================================================
  par(mfrow = c(2,1))

  # Plot trajectories
  plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,1), bty = "n", las = 1,
       xaxt = "n", xlab = "", ylab = "Posterior expected state probability")


  axis_at <- c(1, seq(0, max(n_markers), 10)[-1])
  axis(side = 1, at = axis_at, cex.axis = 0.7) # Marker count
  axis(side = 1, line = 1, at = axis_at, cex.axis = 0.7, # Effective cardinality
       tick = F, labels = sprintf("(%s)", cum_card_eff[,equifs][axis_at]))
  title(main = gsub("_", " ", case))
  title(xlab = "Marker count (cumulative effective cardinality)", line = 3.5)

  for(MOIs in MOIs_per_infection) {
    LTY = ifelse(MOIs == "2_1", 1, 2)
    for(i in 1:n_repeats){
      lines(x = 1:max(n_markers),
            y = sapply(ps_store_all_ms_uniform[[MOIs]][[as.character(i)]], function(x) x[,exp_state]),
            col = cols[i], lwd = 2, lty = LTY)
    }
  }
  legend("topright", col = cols, lwd = 3, title = "Repeats", box.col = "white",
         legend = 1:n_repeats)
  if(length(MOIs_per_infection) > 2) {
    legend("bottomright", lty = c(1,2), lwd = 2, title = "MOIs", box.col = "white",
           legend = gsub("_", " and ", MOIs_per_infection))
  }


  # Plot simplex (helpful when the posterior evades the expected state)
  par(mar = c(0,0,0,0))
  V_labels <- c("Recrudescence", "Relapse", "Reinfection")
  plot_simplex(v_labels = V_labels, v_cutoff = 0.5)
  for(MOIs in MOIs_per_infection) {
    LTY = ifelse(MOIs == "2_1", 1, 2)
    for(i in 1:n_repeats){
      xy_post <- cbind(c(0,0), apply(do.call(rbind, ps_store_all_ms_uniform[[MOIs]][[as.character(i)]]), 1, project2D))
      lines(x = xy_post["x",], y = xy_post["y",], lty = LTY, col = cols[i])
      points(x = xy_post["x",], y = xy_post["y",], pch = "-")
    }
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
         xlab = "Marker count (effective cardinality)")
    title(main = sprintf("Concentration parameter: %s", c), line = 0)
    axis(side = 1, at = n_markers, cex.axis = 0.7)
    axis(side = 1, line = 1, at = n_markers, cex.axis = 0.7, tick = F,
         labels = sprintf("(%s)", round(cum_card_eff[n_markers,as.character(c)])))
    for(MOIs in MOIs_per_infection) {
      for(i in 1:n_repeats) {
        lines(y = post_S[[as.character(c)]][[MOIs]][,i], x = n_markers,
              lwd = 1.5, col = cols[i], pch = 20,
              lty = ifelse(MOIs == "2_1", 1, 2))
        points(y = post_S[[as.character(c)]][[MOIs]][,i], x = n_markers,
               lwd = 1.5, col = cols[i], pch = 20)
      }
    }
  }

  for(MOIs in MOIs_per_infection) {

    if (provide_correct_MOIs) {
      MOIs_num <- as.numeric(strsplit(MOIs, "_")[[1]])
    } else {
      MOIs_num <- determine_MOIs(ys_store[[as.character(c)]][[MOIs]][[as.character(i)]])
    }

    ts_per_gs <- rep(1:length(MOIs_num), MOIs_num)
    gs <- paste0("g", 1:sum(MOIs_num))
    RGs <- enumerate_RGs(MOIs = MOIs_num, igraph = T)
    graph_cols <- colorRamp(RColorBrewer::brewer.pal(n = 9, "Paired"))(seq(0, 1, length.out = length(RGs)))
    graph_cols <- apply(round(graph_cols), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    graph_plot_order <- 1:length(RGs) #c(4, 2, 6, 8, 5, 9, 1, 3, 7)
    names(graph_cols) <- graph_plot_order

    # Plot graphs
    par(mfrow = c(ceiling(sqrt(length(RGs))),ceiling(sqrt(length(RGs)))))

    for(g in graph_plot_order) {
      RG <- RGs[[g]]
      par(mar = c(0.5, 0.5, 0.5, 0.5))
      igraphRG <- RG_to_igraph(RG, gs, ts_per_gs) # Convert to igraph object
      plot_RG(RG =  igraphRG, vertex_palette = "Greys", vertex.label = NA)
      box(col = graph_cols[as.character(g)], lwd = 3)
    }

    # For each m, c combination, plot the graph likelihood
    for(c in c_params){
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
}





