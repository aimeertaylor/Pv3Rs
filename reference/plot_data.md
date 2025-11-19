# Plots the data

Plots allelic data as a grid of coloured rectangles.

## Usage

``` r
plot_data(
  ys,
  fs = NULL,
  person.vert = FALSE,
  mar = c(1.5, 3.5, 1.5, 1),
  gridlines = TRUE,
  palette = RColorBrewer::brewer.pal(12, "Paired"),
  marker.annotate = TRUE,
  legend.lab = "Allele frequencies",
  legend.line = 0.2,
  legend.ylim = c(0.05, 0.2),
  cex.maj = 0.7,
  cex.min = 0.5,
  cex.text = 0.5,
  x.line = 0.2,
  y.line = 2.5
)
```

## Arguments

- ys:

  Nested list of per-person, per-episode, per-marker allelic data; see
  **Examples** and
  [`compute_posterior()`](https://aimeertaylor.github.io/Pv3Rs/reference/compute_posterior.md)
  for the expected per-person structure.

- fs:

  A per-marker list of numeric vectors of allele frequencies. If `NULL`
  (default), only the alleles present in `ys` are shown in the the
  legend, with all per-marker alleles represented equally. Because the
  colour scheme is adaptive, the same allele may have different colours
  given different `ys`. If `fs` is specified, all alleles in `fs`
  feature in the legend with areas proportional to allele frequencies,
  so that common alleles occupy larger areas than rarer alleles. Specify
  `fs` to fix the allele colour scheme across plots of different `ys`.

- person.vert:

  Logical. If `TRUE` (default), person IDs are printed vertically;
  otherwise, they are printed horizontally.

- mar:

  Vector of numbers of lines of margin for the main plot; see `mar`
  entry of [`par`](https://rdrr.io/r/graphics/par.html).

- gridlines:

  Logical. If true (default), white grid lines separating people and
  markers are drawn.

- palette:

  Colour palette for alleles, see the **Value** section of
  [`brewer.pal`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).
  Generally, colours are interpolated: if a marker has `d` possible
  alleles, then the colours used are the `1/(d+1), ..., d/(d+1)`
  quantiles of the palette to ensure that markers with different allele
  counts use different colours.

- marker.annotate:

  Logical. If true (default), the names of the alleles are printed on
  top of their colours in the legend.

- legend.lab:

  Label for the axis of the legend. Defaults to "Allele frequencies".
  Set to `NA` to omit the label; if so, consider adjusting `legend.ylim`
  to use more plotting space.

- legend.line:

  Distance (in character heights) from the colour bar to the legend
  label (defaults to `1.5`).

- legend.ylim:

  Vector specifying the y-coordinate limits of the legend in device
  coordinates (between 0 and 1). Defaults to c(0.05, 0.2).

- cex.maj:

  Numeric; font scaling of major axis labels.

- cex.min:

  Numeric; font scaling of minor axis labels.

- cex.text:

  Numeric; font scaling of the allele labels.

- x.line:

  Distance between top x-axis and x-axis label, defaults to 0.2.

- y.line:

  Distance between left y-axis and y-axis label, defaults to 2.5.

## Value

None

## Details

This function plots alleles (colours), which are observed in different
episodes (columns), on different markers (rows), with episodes grouped
by person. Per-person episodes are plotted from left to right in
chronological order. If multiple alleles are detected for a marker
within an episode, the corresponding grid element is subdivided
vertically into different colours.

By default, markers are ordered lexicographically. If `fs` is provided,
markers are ordered to match the order within `fs`.

The legend depicts the alleles for each marker in the same vertical
order as the main plot. The default colour scheme is adaptive, designed
to visually differentiate the alleles as clearly as possible by
maximizing hue contrast within a qualitative palette. Interpolation is
used to make different colour palettes for markers with different
numbers of possible alleles. The names of the alleles are printed on top
of their colours if `marker.annotate` is set to `TRUE`.

## Examples

``` r
oldpar <- par(no.readonly = TRUE) # Store user's options before plotting

# Running example (runs across compute_posterior, plot_data and plot_simplex)
# based on real data from chloroquine-treated participant 52 of the Vivax
# History Trial (Chu et al. 2018a, https://doi.org/10.1093/cid/ciy319)
ys <- ys_VHX_BPD["VHX_52"] # ys is a list of length one (one participant)
plot_data(ys, fs = fs_VHX_BPD, marker.annotate = FALSE)


# Full data set:
mar <- c(2, 3.5, 1.5, 1) # extra vertical margin for vertical person labels
plot_data(ys = ys_VHX_BPD, person.vert = TRUE, mar = mar, legend.lab = NA)

plot_data(ys = ys_VHX_BPD, person.vert = TRUE, mar = mar, legend.lab = NA,
          fs = fs_VHX_BPD)

plot_data(ys = ys_VHX_BPD, person.vert = TRUE, mar = mar, legend.lab = NA,
          fs = fs_VHX_BPD, marker.annotate = FALSE)


# Demonstrating the adaptive nature of the colour scheme:
y <- ys_VHX_BPD["VHX_52"] # A single person
# Compared to first example, colours now involve only the alleles detected in VHX_52
plot_data(y)


par(oldpar) # Restore user's options
```
