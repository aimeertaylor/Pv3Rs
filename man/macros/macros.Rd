% Macro to create link to an HTML vignette from a help page.
% arg 1: name of vignette (without .html)
% arg 2: (optional) anchor (without leading #)
% arg 3: text to display on the link
% arg 4: (optional) override path to vignette (e.g., "https://yourpkg.github.io/yourpkg/articles/")
% If arg 4 is empty, defaults to Sys.getenv("RdVignettePath", unset = "../doc/")

\newcommand{\HTMLVignette}[4][]{\Sexpr[stage=build,results=rd]{
  local({
    anchor <- "#2"
    if (nchar(anchor))
      anchor <- paste0("#", anchor)
    vigpath <- if (nchar("#1")) {
      "#1"
    } else {
      Sys.getenv("RdVignettePath", unset = "../doc/")
    }
    paste0("\\\\href{", vigpath, "#2.html",
           anchor, "}{#3}")
  })
}}


% Macro to link PDF vignettes, the format is as above
\newcommand{\PDFVignette}{\Sexpr[stage=build,results=rd]{
   local({anchor <- "#2"
          if (nchar(anchor))
            anchor <- paste0("#", anchor)
          vigpath <- Sys.getenv("RdVignettePath", unset = "../doc/")
          paste0("\\\\\\\\href{", vigpath, "#1.pdf",
                 anchor, "}{#3}")
          })}}

% Macro to create link to an HTML vignette in another package on CRAN
% from a help page.
% arg 1 is the package name
% arg 2 is the name of the vignette (without .html)
% arg 3 is an optional anchor within the vignette (without the leading #)
% arg 4 is the text to show on the link.

\newcommand{\HTMLVignetteCRAN}{\Sexpr[stage=build,results=rd]{
  local({anchor <- "#3"
         if (nchar(anchor))
           anchor <- paste0("#", anchor)
         paste0("\\\\\\\\href{https://cran.r-project.org/package=#1/vignettes/#2.html",
                anchor, "}{#4}")
         })}}
