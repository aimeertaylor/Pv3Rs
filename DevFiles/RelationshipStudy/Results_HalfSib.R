################################################################################
# As the number of markers increases, the frequency of relapse classification
# increases, unless intra-episode matches are biased towards rare alleles.
# ==============================================================================
rm(list = ls())

attached <- search()
if("output" %in% attached) detach(output)
if(exists("output")) rm("output")
load("../../data/output_HalfSib.PCLikeSib.rda")
output <- output_HalfSib.PCLikeSib[["Half"]]
attach(output)

par_default <- par(no.readonly = TRUE)
if (n_repeats != 10) stop("Plots assume 10 repeats")
cols <- RColorBrewer::brewer.pal(n = n_repeats, "Paired") # Colours for repeats

# Compute effective cardinality cumulatively
cum_card_eff <- sapply(fs_store, function(fs) {
  cumsum(sapply(fs[marker_rorder], function(x) 1/sum(x^2)))})

# Determine which fs in cum_card_eff are equifrequent
equifs <- which(as.numeric(colnames(cum_card_eff)) > c_cutoff)

# Function to compute the odds of relapse vs reinfection (eq. 4 of halfsibs.tex)
odds_LI = function(m, m2, m3, m4) {
  (2/9) * (1 + 2^(log2(5/2)*m4 + 1)/(2^(m - 2*m2) + 2^(2*m3)))
}

# ==============================================================================
# Process results generated given equifrequent alleles, parents from a single
# population, and all marker counts from one onwards.
# ==============================================================================
# Compute summary statistics of the data: locus types
locus_types <- sapply(1:n_repeats, function(i) {
  y <- ys_store[[as.character(tail(c_params, 1))]][["admixture_FALSE"]][[i]]
  y_summary <- locus_type_summary(y)
})
rownames(locus_types) <- paste0("m", 1:max(n_markers))


# Compute summary statistics of the data: locus type proportions
locus_types_props <- sapply(1:n_repeats, function(i) {
  types <- c("All diff.", "All match", "Intra-match", "Inter-match")
  exp_locus_type_props <- setNames(rep(0,4), types)
  sapply(1:max(n_markers), function(m) {
    x <- table(locus_types[marker_rorder[1:m], i])/m
    exp_locus_type_props[names(x)] <- x
    return(exp_locus_type_props)
  })
}, simplify = F)

# Generate all possible half sib alleles
halfsib_alleles <- enumerate_halfsib_alleles(n_alleles)

# Format into a list to pass to locus_type_summary
halfsib_y <- list(init = apply(halfsib_alleles[,1:2], 1, unique),
                  recur = as.list(halfsib_alleles[,3]))

# Generate locus types
halfsib_locus_types <- locus_type_summary(y = halfsib_y)

# Compute expected locus type proportions
exp_locus_type_props <- table(halfsib_locus_types)/nrow(halfsib_alleles)

# Find first "All diff.' from which to compute odds: should be one for all
m_min_diff <- sapply(1:n_repeats, function(i) {
  min(which(locus_types[marker_rorder,i] == "All diff."))})

# Are m-2m2, 2m3, m4 close to (4/9)*m for large m? Yes
exp_4on9m <- sapply(1:n_repeats, function(i) {
  x <- locus_types[marker_rorder,i]
  m2 <- sum(x == "All match")
  m3 <- sum(x == "Intra-match")
  m4 <- sum(x == "Inter-match")
  m <- length(x)
  c((m - 2*m2), 2*m3, m4)
})
rownames(exp_4on9m) <- c("m-2m2", "2m3", "m4")
dotchart(t(exp_4on9m), color = cols, pch = 20, pt.cex = 1.5, cex = 0.75)
abline(v = (4/9)*max(n_markers))

# Compute odds of relapse:reinfection based on locus types
odds_equ4 <- sapply(1:n_repeats, function(i) {
  sapply(1:max(n_markers), function(m) {
    x <- locus_types[marker_rorder[1:m],i]
    m2 <- sum(x == "All match")
    m3 <- sum(x == "Intra-match")
    m4 <- sum(x == "Inter-match")
    odds_LI(m, m2, m3, m4)
  })
})

# Posterior probabilities based on locus types
post_equ4 <- odds_equ4 / (1 + odds_equ4)

# Compute odds of relapse:reinfection based on posterior probabilities
odds_pv3r <- sapply(1:n_repeats, function(i) {
  odds_model <- sapply(ps_store_all_ms[[as.character(i)]], function(x) x[,"L"]) /
    sapply(ps_store_all_ms[[as.character(i)]], function(x) x[,"I"])
  odds_model[m_min_diff[i]:max(n_markers)]
})

# Denominator condition that 2^(2*m3) >> 2^(m - 2*m3)
denm_equ4 <- sapply(1:n_repeats, function(i) {
  sapply(1:max(n_markers), function(m) {
    x <- locus_types[marker_rorder[1:m],i]
    m2 <- sum(x == "All match")
    m3 <- sum(x == "Intra-match")
    x <- c(2^(2*m3), 2^(m - 2*m3))
    (x/sum(x))[1]
  })
})

# Check no positive recrudescence probability
any(sapply(1:n_repeats, function(i){
  sapply(ps_store_all_ms[[as.character(i)]], function(x) x[,"C"])
}) > 0)


# ==============================================================================
# Plot results generated given equifrequent alleles, parents from a single
# population, and all marker counts from one onwards.
# ==============================================================================
layout(mat = matrix(c(1,2,3,4,4,4), ncol = 2)) # Layout for plots
axis_at <- c(1, seq(0, max(n_markers), 50)[-1])

# Plot the posterior relapse probability trajectories
plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,1), bty = "n", las = 1,
     xaxt = "n", xlab = "Marker count (effective cardinality)",
     ylab = "Posterior relapse probability: Pv3Rs")
abline(h = 2/11, lty = "dashed")
text(x = max(n_markers), y = 2/11, labels = "2/11", pos = 3) # based on 2/9 odds
legend("bottom", col = cols, lwd = 3, inset = 0, legend = 1:n_repeats, bty = "n", cex = 0.5)
axis(side = 1, at = axis_at, cex.axis = 0.7) # Marker count
axis(side = 1, line = 1, at = axis_at, cex.axis = 0.6, # Effective cardinality
     tick = F, labels = sprintf("(%s)", cum_card_eff[,equifs][axis_at]))
for(i in n_repeats:1){
  lines(x = 1:max(n_markers),
        y = sapply(ps_store_all_ms[[as.character(i)]], function(x) x[,"L"]),
        col = cols[i], lwd = 2)
}

# Plot the posterior relapse computed from the odds
plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,1), bty = "n", las = 1,
     xaxt = "n", xlab = "Marker count (effective cardinality)",
     ylab = "Posterior relapse probability: eq. 4")
abline(h = 2/9, lty = "dashed")
text(x = max(n_markers), y = 2/11, labels = "2/11", pos = 3) # based on 2/9 odds
legend("bottom", col = cols, lwd = 3, inset = 0, legend = 1:n_repeats, bty = "n", cex = 0.5)
axis(side = 1, at = axis_at, cex.axis = 0.7) # Marker count
axis(side = 1, line = 1, at = axis_at, cex.axis = 0.6, # Effective cardinality
     tick = F, labels = sprintf("(%s)", cum_card_eff[,equifs][axis_at]))
for(i in n_repeats:1){
  lines(x = 1:max(n_markers),
        y = post_equ4[,i],
        col = cols[i], lwd = 2)
}

# Plot inter-to-intra-match ratio when the first element of
# (2^2m3, 2^(m - 2m2)) normalised exceeds 0.75
plot(NULL, xlim = c(1,max(n_markers)), ylim = c(0,2), bty = "n", las = 1,
     xaxt = "n", xlab = "Marker count",
     ylab = "Intra-to-inter match ratio where 2^2m3 > 2^(m - 2m2)")
abline(h = 0.5*log2(5/2), lty = "dashed") # See latex file
text(x = max(n_markers)-10, y = 0.5*log2(5/2), labels = "0.5*log2(5/2)", pos = 1) # See latex file
axis(side = 1, at = axis_at, cex.axis = 0.7) # Marker count
axis(side = 1, line = 1, at = axis_at, cex.axis = 0.6, # Effective cardinality
     tick = F, labels = sprintf("(%s)", cum_card_eff[,equifs][axis_at]))
for(i in n_repeats:1){
  ratio <- locus_types_props[[i]]["Intra-match", ]/locus_types_props[[i]]["Inter-match", ]
  lines(x = (1:max(n_markers))[denm_equ4[,i] > 0.75],
        y = ratio[denm_equ4[,i] > 0.75], col = cols[i], lwd = 2)
}


# Plot locus type proportion given all markers
locus_types_all <- apply(locus_types, 2, function(x) table(x)/length(x))
dotchart(t(locus_types_all), color = cols, pch = 20, pt.cex = 1.5, cex = 0.75)
title(main = sprintf("Type proportion after %s markers", max(n_markers)))
segments(x0 = exp_locus_type_props["All diff."],
         x1 = exp_locus_type_props["All diff."],
         y0 = 37, y1 = 46, lty = "dashed")
segments(x0 = exp_locus_type_props["All match"],
         x1 = exp_locus_type_props["All match"],
         y0 = 25, y1 = 34, lty = "dashed")
segments(x0 = exp_locus_type_props["Inter-match"],
         x1 = exp_locus_type_props["Inter-match"],
         y0 = 13, y1 = 22, lty = "dashed")
segments(x0 = exp_locus_type_props["Intra-match"],
         x1 = exp_locus_type_props["Intra-match"],
         y0 = 1, y1 = 10, lty = "dashed")

# The proportion of intra match proportions with the number of markers
par(mfrow = c(1,1))
intra_match_props <- sapply(1:max(n_markers), function(m) {
  apply(locus_types, 2, function(x) {
    z <- table(x[1:m])/length(x[1:m])
    ifelse(is.na(z["Intra-match"]), 0, z["Intra-match"])
  })
})
colnames(intra_match_props) <- sprintf("%s markers", 1:max(n_markers))
dotchart(intra_match_props[,seq(10,150,20)], color = cols, pch = 20,
         pt.cex = 1.5, cex = 0.75, labels = rep("", n_repeats),
         main = "Intra-match proportion for different marker counts")
abline(v = exp_locus_type_props["Intra-match"])

# Plot simplex
par(par_default)
par(mar = c(0,0,0,0))
V_labels <- c("Recrudescence", "Relapse", "Reinfection")
plot_simplex(v_labels = V_labels, v_cutoff = 0.5)
for(i in 1:n_repeats){
  xy_post <- cbind(c(0,0), apply(do.call(rbind, ps_store_all_ms[[as.character(i)]]), 1, project2D))
  lines(x = xy_post["x",], y = xy_post["y",], col = cols[i])
  points(x = xy_post["x",], y = xy_post["y",], pch = 20)}

# ==============================================================================
# Process results results for different allele frequency types, for a migrant
# parent (admixture TRUE) as well as parents from a single population, and for
# different marker counts
# ==============================================================================
# Aside: check graphs are all ordered the same [make into a unit test?]
# Extract graph summary (i.e., discard logp)
justRGs <- sapply(ps_store, function(X) {
  sapply(X, function(XX) {
    sapply(XX, function(XXX) {
      sapply(XXX, function(post) {
        sapply(post$RGs, function(RG) c(RG$clone, RG$sib))
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)


# Check all the graphs are returned in the same order
justRG <- justRGs[[1]][[1]][[1]][[1]]
RGcheck <- sapply(c_params, function(c) {
  sapply(c(TRUE, FALSE), function(admixture) {
    sapply(1:n_repeats, function(i) {
      sapply(n_markers, function(m) {
        identical(justRG, justRGs[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]][[as.character(m)]])
      })
    })
  })
})

if (!all(RGcheck)) stop("graphs not returned in the same order")
example_y <- ys_store[[1]][[1]][[1]]
example_MOIs <- determine_MOIs(example_y)
ts_per_gs <- rep(1:length(example_y), example_MOIs)
gs <- paste0("g", 1:sum(example_MOIs))

#-------------------------------------------------------------------------------
# Loop over migrant and non-migrant scenarios
#-------------------------------------------------------------------------------
for(admixture in c("admixture_FALSE", "admixture_TRUE")) {

  par(par_default)

  # Extract probability of relapse
  post_L <- sapply(ps_store, function(X) {
    sapply(X[[admixture]], function(XX) sapply(XX, function(post) post$marg[,"L"]))
  }, simplify = F)

  # Extract exceedance proportion
  prop_L <- sapply(post_L, function(post_L_c) {
    apply(post_L_c, 1, function(post_L_c_m) mean(post_L_c_m > 0.5))
  })

  # Extract probability of the data given the relationship graph
  llikeRGs <- sapply(ps_store, function(X) {
    sapply(X[[admixture]], function(XX) {
      sapply(XX, function(post) sapply(post$RGs, function(RG) RG$logp))
    }, simplify = F)
  }, simplify = F)

  # Plots posterior relapse probabilities
  par(mfrow = c(4,1))
  for(c in c_params){
    plot(NULL, xlim = range(n_markers)+c(-10,10), ylim = c(0,1),
         xaxt = "n", bty = "n", panel.first = grid(nx = NA, ny = NULL),
         ylab = "Posterior relapse probability",
         xlab = "Marker count (effective cardinality)",
         main = sprintf("Concentration parameter: %s", c))
    axis(side = 1, at = n_markers, cex.axis = 0.7)
    axis(side = 1, line = 1, at = n_markers, cex.axis = 0.7, tick = F,
         labels = sprintf("(%s)", round(cum_card_eff[,as.character(c)][n_markers])))
    for(i in 1:n_repeats) {
      lines(y = post_L[[as.character(c)]][,i], x = n_markers,
            lwd = 1.5, col = cols[i], pch = 20)
      points(y = post_L[[as.character(c)]][,i], x = n_markers,
             lwd = 1.5, col = cols[i], pch = 20)
    }
  }

  # Plot exceedance proportion
  par(par_default)
  fields::image.plot(prop_L, axes = FALSE,
                     ylab = "Dirchlet concentration parameter",
                     xlab = "Number of markers",
                     main = "Posterior relapse probability 0.5 exceedence proportion",
                     col = RColorBrewer::brewer.pal(n = 11, "RdBu"),
                     breaks = seq(0,1,length.out = 12))
  axis(at = seq(0,1,length.out = length(n_markers)),
       side = 1, line = -0.5,
       labels = n_markers,
       cex.axis = 1, tick = F)
  axis(at = seq(0,1,length.out = length(c_params)),
       side = 2, line = -0.5, las = 1,
       labels = c_params,
       cex.axis = 1, tick = F)

  # Re-order colours to graph type
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
    igraphRG <- igraph::set_vertex_attr(igraphRG, "name", value = NA) # Remove genotype names
    plot_RG(RG =  igraphRG, vertex_palette = "Greys", labels = NA)
    box(col = graph_cols[as.character(g)], lwd = 3)
  }

  # For each m, c combination, plot the graph likelihood and data
  for(c in c_params){
    par(mfcol = c(n_repeats,length(n_markers)), mar = c(0,0,0,0))
    for(m in n_markers){
      for(i in 1:n_repeats) {
        x <- llikeRGs[[as.character(c)]][[i]][,as.character(m)][graph_plot_order]
        x[x == -Inf] <- NA # Mask -Inf
        x <- x - min(x, na.rm = T) # Re-scale before exponentiating (otherwise all 0)
        x <- exp(x)/sum(exp(x), na.rm = T) # Exponentiate and normalise
        pie(x[!is.na(x)], col = graph_cols[!is.na(x)], labels = NA, border = NA)
        symbols(x = 0, y = 0, circles = c(0.25), inches = FALSE, bg = "white", fg = "white", add = TRUE)
        text(label = round(post_L[[as.character(c)]][as.character(m),i], 2), x = 0, y = 0, cex = 0.5)
      }
      mtext(text = sprintf("concentration %s, marker count %s", c, m), side = 1, cex = 0.5, line = -1)
    }
  }
}

detach("output")


