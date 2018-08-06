#' Hello Country

#' @param countries string vertices countaining countries.
#' @examples
#' hello_country(countries)
#' @export
hello_country <- function(countries) {
    hellos = paste('Hello',countries)
    for (i in hellos ){
        print(i)
    }
}

