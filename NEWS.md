# Pv3Rs (development version)

# Pv3Rs 1.0.0

* CRAN submission post Bioinformatics peer review. 

* README now refers to Bioinformatics article "Pv3Rs: Plasmodium vivax relapse, 
recrudescence, and reinfection statistical genetic inference".

* License is now GPL (>=3) for compatibility with bundled code macros.Rd. This
macro was written by Duncan Murdoch and copied from the rgl package. The rgl
package has a GPL-2 | GPL-3 [expanded from: GPL] license.

* Documentation of `plot_data()`, `compute_posterior()` and `plot_simplex()` now 
features an illustrative example based on real data.

* `plot_simplex()` now features `lim.par` for more control of the white space 
around the simplex.

* In `plot_simplex()`, `p_coords` is now `p.coords`. 

# Pv3Rs 0.0.2

* Users' par options are now restored in functions, vignettes and examples

* CRAN re-submission. First attempt failed automatic checks: par was not 
preceded by graphics:: Second attempt successful.

# Pv3Rs 0.0.1

* Initial CRAN submission. Failed manual checks: intra-function use of par.

