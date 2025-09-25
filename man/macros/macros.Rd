% Copied from macros.Rd of the rgl package, which has a GPL-2 | GPL-3 [expanded from: GPL] license.
% Author of original code: Duncan Murdoch
% Link to original code:
% https://github.com/dmurdoch/rgl/blob/d81b3686e00f2571f5a0e44200cc1ea90e723387/man/macros/macros.Rd#L10-L17

% Macro to create link to an HTML vignette (not an article) from a help page.
% In this version the vignette must be in the same package as the help page
% arg 1 is the name of the vignette (without .html)
% arg 2 is an optional anchor within the vignette (without the leading #)
% arg 3 is the text to show on the link (without quotes)
% Use the RdVignettePath environment variable if vignettes aren't stored
% in ../doc/, e.g. pkgdown puts them in ../articles/.
\newcommand{\HTMLVignette}{\Sexpr[stage=build,results=rd]{
   local({anchor <- "#2"
          if (nchar(anchor))
            anchor <- paste0("#", anchor)
          vigpath <- Sys.getenv("RdVignettePath", unset = "../doc/")
          paste0("\\\\\\\\href{", vigpath, "#1.html",
                 anchor, "}{#3}")
          })}}

