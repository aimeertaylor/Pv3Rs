#' Plots the data
#'
#' Plots the alleles (colours), which are observed in different episodes (rows),
#' on different markers (columns), where episodes are grouped by patient. The
#' patients and per-patient episodes are plotted from bottom to top. If more
#' than one allele is detected per episode per marker, the corresponding
#' row-column entry is subdivided into different colours. The legend depicts the
#' alleles of the markers as the markers appear from left to right in the main
#' plot. Otherwise stated, the legend is ordered by the order of markers stated
#' on on the horizontal axis of the main plot. The colour scheme is adaptive. It
#' is designed to visually differentiate the alleles as much as possible: the
#' maximum range of qualitative scheme, with contrast of hue between adjacent
#' colours, is always used; the adjacent colours are interpolated only if a
#' given marker has more than 12 alleles. The names of the alleles are printed
#' on top of their colours if marker_annotate.
#'
#' @param ys A nested list of per-patient, per-episode, per-marker allelic data.
#'   Specifically, a per-patient list of a per-episode list of a per-marker list
#'   of character vectors of observed alleles.
#' @param fs A per-marker list of numeric vectors of allele frequencies. If NULL
#'   (default), for a given marker, each allele is represented equally in the
#'   marker legend. If specified, legend areas are proportional to allele
#'   frequencies; i.e., common alleles have relatively large legend areas, and
#'   rare alleles have relatively small legend areas.
#' @param marker_alleles A per-marker list of character vectors of all possible
#'   alleles. If NULL (default), only the alleles present in the data are
#'   represented in the legend. Because the colour scheme is adaptive (see
#'   introduction), the same allele will have a different colour in a plot of an
#'   alternative data list if more or fewer alleles are observed at the given
#'   marker across the alternative data list. Specify marker_alleles to fix the
#'   colour of a given allele across plots of different data list, thereby
#'   facilitating cross-comparison.
#'@param marker_annotate Logical. If true (default), the names of the alleles
#'  are printed on top of their colours in the legend.
#'
#' @examples
#' # Examples
#'
#' # Plot example Plasmodium vivax data set
#' plot_data(ys = ys_VHX_BPD)
#' plot_data(ys = ys_VHX_BPD, fs = fs_VHX_BPD)
#' plot_data(ys = ys_VHX_BPD, fs = fs_VHX_BPD, marker_annotate = F)
#'
#' # Demonstrating the adaptive nature of the colour scheme:
#' ys <- ys_VHX_BPD["VHX_52"] # A single patient with
#' plot_data(ys) # Legend adapts to alleles detected in VHX_52 only
#' plot_data(ys, fs = fs_VHX_BPD, marker_alleles = lapply(fs_VHX_BPD, names))
#'
#'
#' @export
plot_data = function(ys,
                     fs = NULL,
                     marker_alleles = NULL,
                     marker_annotate = TRUE){

  # Function to create ramped colours based on "Paired" brewer.pal
  cols <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))

  # Extract patients and episode counts and names
  n_patients <- length(ys)
  n_episodes <- sum(unlist(lapply(ys, length)))
  patient_names <- names(ys)
  episode_names <- unlist(lapply(ys, names), use.names = FALSE)

  # Name patients if unnamed or not uniquely named
  if (is.null(patient_names) | !identical(unique(patient_names), patient_names)) {
    names(ys) <- paste0("p", 1:n_patients)
  }

  # Name episodes if unnamed or not uniquely named
  if (is.null(episode_names) | !identical(unique(episode_names), episode_names)) {
    lapply(patient_names, function(patient_name) {
      y <- ys[[patient_name]]
      names(y) <- paste0(patient_name, 1:length(y))
    })
  }

  # Extract patient and episode unique identifiers
  pids <- names(ys)
  eids <- unlist(lapply(ys, names), use.names = FALSE)

  # Remove patient-level list - redundant since episodes now have unique identifiers
  ys_by_epi <- unlist(ys, recursive = F, use.names = FALSE)
  names(ys_by_epi) <- eids

  # Extract marker info.
  markers_in_data <- unique(unlist(lapply(ys_by_epi, function(epi) names(epi)), use.names = F))
  marker_alleles_in_data <- lapply(markers_in_data, function(marker) {
    sort(
      unique(
        unlist(
          lapply(ys_by_epi, function(epi) epi[[marker]]),
          use.names = FALSE)))
  })
  names(marker_alleles_in_data) <- markers_in_data

  if (is.null(marker_alleles)) {
    markers <- markers_in_data
    marker_alleles <- marker_alleles_in_data
  } else {
    markers <- names(marker_alleles)
    markers_agree <- all(markers_in_data %in% markers)
    alleles_agree <- all(sapply(markers, function(marker) all(marker_alleles_in_data[[marker]] %in% marker_alleles[[marker]])))
    if (!markers_agree | !alleles_agree) stop("Marker information disagreement")
  }

  marker_cardinalities <- sapply(marker_alleles, length)
  n_markers = length(markers)

  # Extract MOIs
  mois <- unlist(sapply(ys_by_epi, function(epi) {
    moi <- max(sapply(epi, function(marker) {
      length(unique(marker[!is.na(marker)]))
    }))}), use.names = FALSE)
  names(mois) <- eids
  maxMOI =  max(mois) # Maximum MOI

  # Reformat data s.t. columns house alternative alleles in polyclonal samples
  factorial_maxMOI = factorial(maxMOI) # For expanding to wide format
  markers_wide = paste(rep(markers, each = factorial_maxMOI), rep(1:factorial_maxMOI, n_markers), sep = '_')
  n_markers_wide = length(markers_wide) # equal to n_markers * factorial_maxMOI
  marker_data_wide_format = array(dim = c(n_episodes, n_markers_wide), dimnames = list(eids, markers_wide))

  for(epi in eids){
    for(marker in markers){
      epi_marker_alleles = sort(unique(ys_by_epi[[epi]][[marker]]))
      if (all(is.na(epi_marker_alleles))) next # TRUE if typed but all NA or untyped
      moi_marker = length(epi_marker_alleles[!is.na(epi_marker_alleles)])
      to_fill_list = lapply(seq(factorial_maxMOI, 1, -factorial_maxMOI/moi_marker), function(x){1:x})
      for(k in 1:moi_marker){
        to_fill_list_k <- to_fill_list[[k]]
        to_fill = paste(rep(marker, each = max(to_fill_list_k)), to_fill_list_k, sep = '_')
        marker_data_wide_format[epi, to_fill] = rep(epi_marker_alleles[k], each = max(to_fill_list_k))
      }
    }
  }

  # Re-colour to maximise contrast per marker
  marker_data_wide_factor = array(dim = dim(marker_data_wide_format), dimnames = dimnames(marker_data_wide_format))
  for(marker in markers){
    ind = grepl(paste(marker, '_', sep = ''), colnames(marker_data_wide_format))
    marker_data_wide_factor[,ind] = apply(marker_data_wide_format[, ind, drop = FALSE], 2, function(x){
      y = as.numeric(factor(x, ordered = TRUE, levels = marker_alleles[[marker]]))
      return(y) # The level rather than the allele
    })
  }

  # Create an index to visually group individuals
  IDs_factor <- as.factor(rep(pids, sapply(ys, length)))
  IDs = as.numeric(IDs_factor)
  ID_01 = c(0, rep(NA, n_episodes-1))
  if (n_episodes > 1) {
    for(i in 2:n_episodes){
      if(IDs[i] == IDs[i-1]){
        ID_01[i] = ID_01[i-1]
      } else { # switch
        ID_01[i] = abs(ID_01[i-1] - 1)
      }
    }
  }

  # Y axis
  ID_midpoints = rep(NA, length(unique(IDs)))
  names(ID_midpoints) = unique(IDs)
  for(ID in unique(IDs)){
    inds = which(IDs == ID)
    inds_mapped_01 = (inds - 0)/(n_episodes + 1 - 0)
    ID_midpoints[as.character(ID)] = median(inds_mapped_01)
  }

  # Plot
  par(mfrow = c(1,1),
      fig = c(0,0.75,0.03,0.97), # fig = c(x1, x2, y1, y2) (to leave room for the legend)
      mar = c(4,4,0,0.5)) # mar = c(bottom, left, top, right)

  # Add IDs to each vertical edge
  Empty_array = array(dim = dim(marker_data_wide_factor))
  matrix_to_plot = t(rbind(-1, cbind(ID_01, Empty_array, ID_01), -1))
  image(matrix_to_plot, ylim = c(0,1), col = grey(c(0,0.35,0.75)), axes = FALSE)

  # Add marker data
  COLs = seq(1, n_markers_wide, factorial_maxMOI)
  for(i in 1:length(COLs)){
    Y = marker_data_wide_factor
    Y[,-(COLs[i]:(COLs[i]+(factorial_maxMOI)-1))] = NA
    if (all(is.na(Y))) next # If no data skip to next maker
    image(t(rbind(-1, cbind(NA, Y, NA), -1)), col = c(grey(0), cols(marker_cardinalities[i])),
          axes = FALSE, add = TRUE, zlim = c(-1, marker_cardinalities[i]))
  }

  # Add axes
  xaxis = seq(0, 1, length.out = n_markers_wide+2)
  xaxis_diff = (xaxis[2] - xaxis[1])
  if(maxMOI == 1){
    xaxis_at = xaxis[2:(n_markers+1)]
  } else {
    xaxis_at = xaxis[seq((1 + factorial_maxMOI/2), n_markers_wide, factorial_maxMOI)] + ifelse((factorial_maxMOI %% 2) == 0, xaxis_diff/2, 1)
  }
  axis(at = xaxis_at, side = 1, line = -0.5, labels = markers, cex.axis = 0.5, tick = F)
  axis(at = xaxis_at, side = 1, line = 0.5, labels = paste('(', marker_cardinalities, ')', sep = ''), cex.axis = 0.5, tick = F)
  axis(side = 2, at = ID_midpoints, labels = levels(IDs_factor)[unique(IDs)], cex.axis = 0.5, las = 2, lwd.ticks = 0.25, lwd = 0)
  axis(side = 1, at = c(0, tail(xaxis, 1)), line = -0.5, labels = rep('Grouping', 2), las = 2, cex.axis = 0.5, tick = F)
  title(ylab = bquote(.(n_episodes) ~ 'episodes (one row per episode, grouped by patient ID)'),
        line = 3.3, cex.lab = 0.5)
  title(xlab = 'Marker (number of distinct alleles)', line = 3, cex.lab = 0.5)

  # Add legend
  par(fig = c(0,1,0,1)) # Critical to avoid odd placement of first legend column
  legendwidth = (1-0.76)/n_markers
  for(marker in markers){ # For each marker in turn
    i = which(marker == markers)
    if(i == 1){
      SMALLplot = c(0.75, 0.75 + legendwidth, 0.03, 0.97)
    } else {
      SMALLplot = c(0.75 + legendwidth * (i-1), 0.75 + legendwidth * i, 0.03, 0.97)
    }
    if (marker_cardinalities[marker] == 1) {
      Breaks <- c(0,1)
      At <- 0.5
    } else {
      if(is.null(fs)) {
        Breaks <- NULL
        At <- seq(0, 1, length.out = marker_cardinalities[marker])
      } else {
        # Only use the frequencies whose alleles are in marker_alleles
        Breaks <- c(0, cumsum(fs[[marker]][as.character(marker_alleles[[marker]])]))
        At <- Breaks[-length(Breaks)] + diff(Breaks)/2
      }
    }
    matrix_to_plot[is.na(matrix_to_plot)] <- 0
    if (marker_annotate) {
      Labels <- marker_alleles[[marker]]
    } else {
      Labels <- rep(NA, marker_cardinalities[marker])
    }
    fields::image.plot(abs(matrix_to_plot), # Two z values: 0 and 1
                       col = cols(marker_cardinalities[marker]),
                       breaks = Breaks,
                       legend.only = TRUE,
                       add = TRUE,
                       axis.args = list(at = At,
                                        labels = Labels,
                                        cex.axis = 0.5,
                                        tick = FALSE,
                                        line = -1.5,
                                        hadj = 0.5),
                       smallplot = SMALLplot, legend.mar = 0)
  }
}
