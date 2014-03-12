manuscript.code <- function() {
    fig.code <- file.path(system.file(package="LakeMetabolizer"), "manuscript.script.R")
    file.show(fig.code)
}
