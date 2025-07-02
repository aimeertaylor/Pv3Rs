#' Plots the data
#'
#' Plots the alleles (colours), which are observed in different episodes
#' (columns), on different markers (rows), where episodes are grouped by
#' patient. The per-patient episodes are plotted from left to right in
#' chronological order. If more than one allele is detected per episode per
#' marker, the corresponding row is subdivided into rows of different colours.
#' The order of markers follow the allele frequencies, if provided. Otherwise,
#' markers are ordered lexicographically.
#'
#' The legend depicts the alleles of the markers in the same vertical order as
#' the main plot. The default colour scheme is adaptive. It is designed to
#' visually differentiate the alleles as much as possible: the maximum range of
#' qualitative scheme, with contrast of hue between adjacent colours.
#' Interpolation is used to make different colour palettes for markers with
#' different numbers of possible alleles. The names of the alleles are printed
#' on top of their colours if \code{marker.annotate} is set to \code{TRUE}.
#'
#'@param ys A nested list of per-patient, per-episode, per-marker allelic data.
#'  Specifically, a per-patient list of a per-episode list of a per-marker list
#'  of character vectors of observed alleles.
#'@param fs A per-marker list of numeric vectors of allele frequencies. If NULL
#'  (default), for a given marker, only the alleles present in the data are
#'  represented in the legend, and each allele is represented equally. Because
#'  the colour scheme is adaptive (see introduction), the same allele will have
#'  a different colour in a plot of an alternative data list if more or fewer
#'  alleles are observed at the given marker across the alternative data list.
#'  If fs is specified, all possible alleles are represented and legend areas
#'  are proportional to allele frequencies; i.e., common alleles have relatively
#'  large legend areas, and rare alleles have relatively small legend areas.
#'  Specify fs to fix the colour of a given allele across plots of different
#'  data lists, thereby facilitating cross-comparison.
#'@param patient.vert Logical. If true, patient IDs are printed vertically;
#'  otherwise they are printed horizontally (default).
#'@param mar A vector which gives the number of lines of margin for the main
#'  plot; see the entry for `mar` in \code{\link[graphics]{par}}.
#'@param gridlines Logical. If true (default), white gridlines separating
#'  patients and markers will be drawn.
#'@param palette Colour palette for alleles, see return value of
#'  \code{\link[RColorBrewer]{brewer.pal}}. Generally, the colours are
#'  interpolated. If a marker has \code{d} possible alleles, then the colours
#'  used are the \code{1/(d+1), ..., d/(d+1)} quantiles of the colours in
#'  \code{palette}. This is done such that markers with different number of
#'  possible alleles would generally use different colours.
#'@param marker.annotate Logical. If true (default), the names of the alleles
#'  are printed on top of their colours in the legend.
#'@param legend.lab Label for the axis of the color legend, defaults to "Allele frequencies".
#'  Set to \code{NA} if not desired. If this is the case, consider adjusting
#'  \code{legend.ylim} to use more plotting space.
#'@param legend.line Distance in units of character height of the legend label
#'  from the colour bar, defaults to 1.5.
#'@param legend.ylim Vector of two floats representing the y-coordinate limits
#'  of the legend in device coordinates (between 0 and 1). Defaults to c(0.05, 0.2).
#'@param cex.maj Float representing the font scaling of major axis labels.
#'@param cex.min Float representing the font scaling of minor axis labels.
#'@param cex.text Float representing the font scaling for the allele labels.
#'@param x.line Distance between top x-axis and x-axis label, defaults to 0.2.
#'@param y.line Distance between left y-axis and y-axis label, defaults to 2.5.
#'
#' @examples
#'
#' # Plot example Plasmodium vivax data set
#' mar <- c(2, 3.5, 1.5, 1) # extra vertical margin for vertical patient labels
#' plot_data(ys = ys_VHX_BPD, patient.vert = TRUE, mar = mar, legend.lab = NA)
#' plot_data(ys = ys_VHX_BPD, fs = fs_VHX_BPD, patient.vert = TRUE, mar = mar,
#'           legend.lab = NA)
#' plot_data(ys = ys_VHX_BPD, fs = fs_VHX_BPD, marker.annotate = FALSE,
#'           patient.vert = TRUE, mar = mar, legend.lab = NA)
#'
#' # Demonstrating the adaptive nature of the colour scheme:
#' ys <- ys_VHX_BPD["VHX_52"] # A single patient
#' plot_data(ys, fs = fs_VHX_BPD) # Colours and the legend match plots above
#' plot_data(ys) # Colours and the legend adapt to alleles detected in VHX_52
#'
#'
#' @export
plot_data = function(ys,
                     fs = NULL,
                     patient.vert = FALSE,
                     mar = c(1.5, 3.5, 1.5, 1),
                     gridlines = TRUE,
                     palette = RColorBrewer::brewer.pal(12, "Paired"),
                     marker.annotate = TRUE,
                     legend.lab = "Allele frequencies",
                     legend.line = 0.2,
                     legend.ylim = c(0.05, 0.2),
                     cex.maj = 0.7, cex.min = 0.5, cex.text = 0.5,
                     x.line = 0.2, y.line = 2.5
){
  # Function to create ramped colours based on "Paired" brewer.pal
  cols <- grDevices::colorRampPalette(palette)

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
    ys <- sapply(patient_names, function(patient_name) {
      y <- ys[[patient_name]]
      names(y) <- paste(patient_name, 1:length(y), sep = "_")
      return(y)
    }, USE.NAMES = T, simplify = F)
  }

  # Extract patient and episode unique identifiers
  pids <- names(ys)
  eids <- unlist(lapply(ys, names), use.names = FALSE)

  # Remove patient-level list - redundant since episodes now have unique identifiers
  ys_by_epi <- unlist(ys, recursive = F, use.names = FALSE)
  names(ys_by_epi) <- eids

  # Extract marker info.
  markers_in_data <- sort(unique(unlist(lapply(ys_by_epi, function(epi) names(epi)), use.names = F)))
  marker_alleles_in_data <- lapply(markers_in_data, function(marker) {
    sort(
      unique(
        unlist(
          lapply(ys_by_epi, function(epi) epi[[marker]]),
          use.names = FALSE)))
  })
  names(marker_alleles_in_data) <- markers_in_data

  if (is.null(fs)) {
    markers <- markers_in_data
    marker_alleles <- marker_alleles_in_data
  } else {
    marker_alleles = lapply(fs, names)
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
      seglen <- factorial_maxMOI/moi_marker
      for(k in 1:moi_marker){
        seg <- (k*seglen-seglen+1):(k*seglen)
        to_fill = paste(rep(marker, seglen), seg, sep = '_')
        marker_data_wide_format[epi, to_fill] = rep(epi_marker_alleles[k], seglen)
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
  IDs <- as.vector(IDs_factor)
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

  # episode axis, note that patient IDs are unique
  ID_midpoints <- rep(NA, length(pids))
  names(ID_midpoints) <- pids
  for(ID in pids){
    inds = which(IDs == ID)
    ID_midpoints[as.character(ID)] <- (stats::median(inds)-0.5)/n_episodes
  }

  # Plot
  graphics::par(mfrow = c(1,1),
                fig = c(0,1,legend.ylim[2]+0.01,1), # fig = c(x1, x2, y1, y2) (to leave room for the legend)
                mar = mar) # mar = c(bottom, left, top, right)
  figsize <- graphics::par("fin")
  n_x <- nrow(marker_data_wide_factor)
  n_y <- ncol(marker_data_wide_factor)
  x_coords <- seq(1, 2*n_x-1, 2)/(2*n_x)
  y_coords <- seq(1, 2*n_y-1, 2)/(2*n_y)

  # Add IDs to each vertical edge
  rim_ratio <- 40
  rim_coords <- seq(-1, 2*rim_ratio+1, 2)/(2*rim_ratio)
  Empty_array <- array(dim = c(rim_ratio, n_x))
  matrix_to_plot <- t(rbind(ID_01, Empty_array, ID_01))
  graphics::image(
    x_coords, rim_coords, matrix_to_plot, xlab="", ylab="",
    col = grDevices::grey(c(0.4,0.7)), axes = FALSE,
    ylim=c(1+1/rim_ratio,-1/rim_ratio))

  # Add marker data
  COLs = seq(1, n_markers_wide, factorial_maxMOI)
  for(i in 1:length(COLs)){
    Y = marker_data_wide_factor
    Y[,-(COLs[i]:(COLs[i]+(factorial_maxMOI)-1))] = NA
    if (all(is.na(Y))) next # If no data skip to next maker
    n_alleles <- marker_cardinalities[i]
    #graphics::image(t(rbind(-1, cbind(NA, Y, NA), -1)),
    graphics::image(x_coords, y_coords, Y,
                    col = c(grDevices::grey(0), cols(n_alleles+2)[(1:n_alleles)+1]),
                    axes = FALSE, add = TRUE, zlim = c(-1, n_alleles))
  }

  # Add gridlines
  if(gridlines) {
    for (i in 1:n_episodes) {
      x <- i/n_episodes
      graphics::segments(x, -1/rim_ratio, x, 0, lwd=0.5, col="white")
      graphics::segments(x, 1, x, 1+1/rim_ratio, lwd=0.5, col="white")
    }
    graphics::abline(v=c(0, cumsum(unlist(lapply(ys, length))))/n_episodes, lwd=0.5, col="white")
    graphics::abline(h=(0:n_markers)/n_markers, lwd=0.5, col="white")
    # graphics::box(lwd=0.5)
  }

  # Add axes
  graphics::mtext(text = paste0(markers,'\n','(',marker_cardinalities,')'), side = 2, at = ((1:n_markers)-0.5)/(n_markers),
                  line = 0.2, cex = cex.min, las = 1)
  graphics::mtext(text = rep('Grouping', 2), side = 2, at = c(-0.5,rim_ratio+0.5)/rim_ratio,
                  line = 0.2, cex = cex.min, las = 2)
  graphics::mtext(text = pids, side = 1, at = ID_midpoints, cex = cex.min,
                  las = ifelse(patient.vert, 3, 1), line = ifelse(patient.vert, 0.2, 0))
  xlabel <- ifelse(n_patients > 1,
                   sprintf("%s episodes in %s people (one column per episode, grouped by person; for each person, episodes are ordered chronologically from left to right)", n_episodes, n_patients),
                   sprintf("%s episodes in 1 person (one column per episode, ordered chronologically from left to right)", n_episodes))
  graphics::mtext(text = xlabel, line = x.line, cex = cex.maj)
  graphics::title(ylab = 'Marker (number of distinct alleles)',
                  line = y.line, cex.lab = cex.maj)

  # Add legend
  graphics::par(fig = c(0,1,0,1)) # Critical to avoid odd placement of first legend column
  legendsize <- legend.ylim[2] - legend.ylim[1]
  rowsize = legendsize/n_markers
  for(i in 1:n_markers){ # For each marker in turn
    marker <- markers[i]
    SMALLplot = c(0.05, 0.95,
                  legend.ylim[1] + rowsize*(n_markers-i), legend.ylim[1] + rowsize*(n_markers-i+1))
    if (marker_cardinalities[marker] == 1) {
      Breaks <- c(0,1)
      At <- 0.5
    } else {
      if(is.null(fs)) {
        Breaks <- NULL
        At <- ((1:marker_cardinalities[marker]) - 0.5) / marker_cardinalities[marker]
      } else {
        # Only use the frequencies whose alleles are in marker_alleles
        Breaks <- c(0, cumsum(fs[[marker]][as.character(marker_alleles[[marker]])]))
        At <- Breaks[-length(Breaks)] + diff(Breaks)/2
      }
    }
    if (any(diff(Breaks) < .Machine$double.eps)) {
      stop("Some near-zero frequencies: set fs to NULL to enable plotting")
    }
    matrix_to_plot[is.na(matrix_to_plot)] <- 0
    if (marker.annotate) {
      Labels <- marker_alleles[[marker]]
      # let n be the number of characters in the allele label
      min_spaces <- 0.01*sapply(Labels, nchar)
      # labels must be vertically apart by at least x% of the legend height
      pos_vec <- c(-min_spaces[1]/2, 1 + utils::tail(min_spaces, 1)/2) # and x/2% away from the top and bottom borders
      marker_order <- ifelse(is.null(fs), 1:length(Labels), order(fs[[marker]], decreasing=TRUE))
      for(j in marker_order) {
        if (min(abs(pos_vec - At[j])) < min_spaces[j]) { # label would be too close to an existing label
          Labels[j] <- NA
        } else {
          pos_vec <- c(pos_vec, At[j])
        }
      }
    } else {
      Labels <- rep(NA, marker_cardinalities[marker])
    }
    n_alleles <- marker_cardinalities[marker]
    fields::image.plot(abs(matrix_to_plot), # Two z values: 0 and 1
                       col = cols(n_alleles+2)[(1:n_alleles)+1],
                       breaks = Breaks,
                       legend.only = TRUE,
                       legend.lab = ifelse(i == n_markers, legend.lab, NA),
                       legend.line = legend.line,
                       legend.cex = cex.maj,
                       horizontal = TRUE,
                       add = TRUE,
                       axis.args = list(cex.axis = 0.5, tick = FALSE, labels = FALSE),
                       smallplot = SMALLplot, legend.mar = 0)
    y_ndc <- mean(SMALLplot[3:4])
    x_ndc <- SMALLplot[1] + (SMALLplot[2]-SMALLplot[1])*At
    text_ndc(x_ndc, y_ndc, Labels, adj=c(0.5,0.5), cex=cex.text)
  }
}

# Plot text in normalised device coordinates (NDC)
text_ndc <- function(x_ndc, y_ndc, labels, ...) {
  x_user <- graphics::grconvertX(x_ndc, from = "ndc", to = "user")
  y_user <- graphics::grconvertY(y_ndc, from = "ndc", to = "user")
  graphics::text(x = x_user, y = y_user, labels = labels, xpd = NA, ...)
}
