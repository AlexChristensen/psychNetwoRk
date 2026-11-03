.onload <- function(libname, pkgname)
{
    library.dynam("psychNetwoRk",package=pkgname,lib.loc=libname)
}

.onAttach <- function(libname, pkgname)
{
    msg <- styletext(styletext(paste("\npsychNetwoRk (version ", utils::packageVersion("psychNetwoRk"), ")", sep = ""), defaults = "underline"), defaults = "bold")
    msg <- paste(msg,"\n\nFor bugs and errors, submit an issue to <https://github.com/AlexChristensen/psychNetwoRk/issues>")

    packageStartupMessage(msg)
}
