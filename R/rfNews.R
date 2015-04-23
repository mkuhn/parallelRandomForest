rfNews <- function() {
    newsfile <- file.path(system.file(package="parallelRandomForest"), "NEWS")
    file.show(newsfile)
}
