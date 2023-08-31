#' Wrapper for the 'laplace.test' function from the 'lawstat' package
#'
#' This function is a wrapper for the 'laplace.test' function from the 'lawstat' package.
#' The advantage of using this function instead of the 'laplace.test' function of the
#' 'lawstat' package is that this function shows the tables with the critical values
#' of the statistics provided by Puig and Stephens (2000). This makes interpretation of the results easier.
#'
#' @param y A numeric vector containing the sample data.
#' @param print_tables A boolean variable that indicates whether tables with critical values will be printed or not.
#' @returns The result of the laplace.test function.
#' @importFrom lawstat laplace.test
#' @import greekLetters
#' @export
#'
#' @seealso
#' \code{\link[lawstat]{laplace.test}} Goodness-of-fit Test Statistics for the Laplace Distribution from package \code{lawstat}.
#'
#' @references Puig, P. and Stephens, M. A. (2000). Tests of fit for the Laplace distribution, with applications.
#' \emph{Technometrics}, \strong{42}(4), 417-424. \doi{10.2307/1270952}.
#'
#' @examples
#' \donttest{
#' normal_sample = rnorm(100, 0, 10)
#' laplace_sample = rlaplace(100, 0, 10)
#' laplace.dist.test(normal_sample)
#' laplace.dist.test(laplace_sample)
#' }
#'
laplace.dist.test = function(y, print_tables = TRUE) {

  #if(is.null(print_tables)){
  #  print_tables = TRUE
  #}
  if(print_tables){

    # Criar um novo ambiente local
    local_env <- new.env()

    data_file <- system.file("extdata", file = "tabelas_critical_points.rda", package = "diagL1")
    # Carregar o dataframe no ambiente local
    load(data_file, envir = local_env)

    tabela_A2 = local_env$tabela_A2
    tabela_W2 = local_env$tabela_W2
    tabela_U2 = local_env$tabela_U2
    tabela_sqrtn_D = local_env$tabela_sqrtn_D
    tabela_sqrtn_V = local_env$tabela_sqrtn_V

    # Print tables with critical values
    names(tabela_A2) = c( "n", paste(greekLetters::greeks("alpha")," = 0.5", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.25", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.1", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.05", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.025", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.01", sep =""))
    message("Critical values of statistic A2, provided by Puig and Stephens (2000).", "\n")
    print(tabela_A2)
    cat("\n")

    names(tabela_W2) = c( "n", paste(greekLetters::greeks("alpha")," = 0.5", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.25", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.1", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.05", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.025", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.01", sep =""))
    message("Critical values of statistic W2, provided by Puig and Stephens (2000).", "\n")
    print(tabela_W2)
    cat("\n")

    names(tabela_U2) = c( "n", paste(greekLetters::greeks("alpha")," = 0.5", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.25", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.1", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.05", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.025", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.01", sep =""))
    message("Critical values of statistic U2, provided by Puig and Stephens (2000).", "\n")
    print(tabela_U2)
    cat("\n")

    names(tabela_sqrtn_D) = c( "n", paste(greekLetters::greeks("alpha")," = 0.5", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.25", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.1", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.05", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.025", sep =""),
                          paste(greekLetters::greeks("alpha")," = 0.01", sep =""))
    message("Critical values of statistic D, provided by Puig and Stephens (2000).", "\n")
    print(tabela_sqrtn_D)
    cat("\n")

    names(tabela_sqrtn_V) = c( "n", paste(greekLetters::greeks("alpha")," = 0.5", sep =""),
                               paste(greekLetters::greeks("alpha")," = 0.25", sep =""),
                               paste(greekLetters::greeks("alpha")," = 0.1", sep =""),
                               paste(greekLetters::greeks("alpha")," = 0.05", sep =""),
                               paste(greekLetters::greeks("alpha")," = 0.025", sep =""),
                               paste(greekLetters::greeks("alpha")," = 0.01", sep =""))
    message("Critical values of statistic V, provided by Puig and Stephens (2000).", "\n")
    print(tabela_sqrtn_V)
    cat("\n\n")
  }

  result = lawstat::laplace.test(y)
  return(result)
}
