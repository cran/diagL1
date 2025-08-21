#main author: Kévin Allan Sales Rodrigues

#' Predictive Influence Function
#'
#' @param y A vector with response variables.
#' @param x A matrix with a single explanatory variable.
#' @param w Explanatory variables vector for z, usually explanatory variables averages.
#' @param num_cores Number of cores you want to use for parallel processing (default = 2).
#'
#' @returns
#' A vector with predictive influence function value for each observation.
#'
#' @importFrom utils capture.output
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom cubature hcubature
#' @importFrom foreach foreach %dopar%
#' @export
#' @details Please, install and load the "foreach" package to use this function. For more details see Rodrigues (2024).
#'
#' @references Rodrigues, K. A. S. (2024). \strong{Analysis of the adjustment of the L1 regression model}.
#' Phd dissertation, University of São Paulo, BR.
#'
#' Rodrigues, K. A. S., Elian, S. N., & Pereira, G. H. A. (2025). Homoscedasticity tests for L1 regression and their performance evaluation through simulations. \emph{Statistics}, Advance online publication. https://doi.org/10.1080/02331888.2025.2536097
#'
#' Rodrigues, K. A. S., & Elian, S. N. (2025). Influence measures for L1 regression: an analysis with the R package diagL1. \emph{Journal of Applied Statistics}, Advance online publication. https://doi.org/10.1080/02664763.2025.2510691
#'
#'
# #' @examples

### Aplicações da tese

# dados sobre Bile
#yy = bile[,2]
#xx = bile[,1]
#PIF(yy, xx, w = c(15.24828), num_cores = 4)

# dados sobre Incêndios
#yyy = Fire[,3]
#xxx = Fire[,5:6]
#PIF(yyy, xxx, w = c(32.3617, 10.69583), num_cores = 4)

# Observações
# os valores de W foram obtidos por meio das médias das variáveis explicativas.
# parallel::detectCores() verifica o número de núcleos do processador do computador.

# Próximos passos, realizar o "grosso" da computação em uma linguagem compilada, como C++.
# Realizar a computação paralela em C++.

PIF = function(y, x, w, num_cores = 2){

  ##### constant, FVP, FVP_j [functions]

  ### function of constants with respect to z
  constant = function(n){2^(-n-1)*(1/(n+1))^(-n-1)*exp(-n-1)}

  ### numeric FVPP (all observations)
  FVP = function(y, z_values, explanatory_variables, W){
    N = length(y)

    FVPP = cbind(z_values,0) #initializing matrix with FVPPs
    i=1 # "z_values" vector index counter

    for(z in z_values){

      capture.output(                  # to avoid polluting the screen
        if (!is.null(dim(explanatory_variables)) && dim(explanatory_variables)[2] > 1){ # multivariate case
          FVPP[i,2] <- constant(N)*(
            sum(abs(
              regL1(c(y, z)       #fitting model with z observation
                    ~ rbind(explanatory_variables, W) )$residuals
            ))
          )^(-N-1)

        }else{ # univariate case
          FVPP[i,2] <- constant(N)*(
            sum(abs(
              regL1(c(y, z)       #fitting model with z observation
                    ~ c(explanatory_variables, W) )$residuals
            ))
          )^(-N-1)
        })

      i = i+1 # increase counter
    }
    return(FVPP)
  }

  ### numeric FVPP_j (without j-th observation)
  FVP_j = function(y, z_values, explanatory_variables, W){
    N = length(y)
    tam_z = length(z_values) # number of values for z
    FVP_j_final = NULL  # final matrix
    excluded_index = NULL # vector of excluded indices

    for(j in 1:N){

      excluded_index = c(excluded_index, rep(j, tam_z))

      FVPP_ = cbind(z_values,0) #initializing matrix with FVPPs
      i=1 # index counter

      for(z in z_values){

        capture.output(
          if (!is.null(dim(explanatory_variables)) && dim(explanatory_variables)[2] > 1){ # multivariate case
            FVPP_[i,2] <- constant(N-1)*(
              sum(abs(
                regL1(c(y[-j], z)          #fitting model with observation z and without j
                      ~rbind(explanatory_variables[-j,], W) )$residuals
              ))
            )^(-N-1+1)

          }else{ # univariate case
            FVPP_[i,2] <- constant(N-1)*(
              sum(abs(
                regL1(c(y[-j], z)          #fitting model with observation z and without j
                      ~c(explanatory_variables[-j], W) )$residuals
              ))
            )^(-N-1+1)
          }
        )

        i = i+1 # increase counter
      }
      FVP_j_final = rbind(FVP_j_final, FVPP_)
    }

    FVP_j_final = cbind(FVP_j_final, excluded_index) # adding index excluded
    colnames(FVP_j_final) = c("z_values", "FVPPj", "j")
    return(FVP_j_final)
  }

  #####BEGIN
  options(warn = -1)

  # detectando números de núcleos
  numCores = num_cores #parallel::detectCores() #verifica o número de núcleos do processador
  # Criar um cluster com o número especificado de núcleos
  cluster = parallel::makeCluster(numCores)
  # Registrar o cluster para uso com doParallel
  doParallel::registerDoParallel(cluster)

  #####################Calcular FVPPjs -> trocar dados

  #valores dos hiperparâmetros (número de observações)
  number_of_observations = 1:length(y)
  k = number_of_observations

  const_int <- foreach::foreach(k=number_of_observations,
                                #.export = c("FVP_j", "y", "x"),
                                .combine='rbind', .packages = c("diagL1", "cubature")) %dopar% {

         ##uso as médias das variáveis explicativas como vetor W.
         #x_col_mean_old = as.numeric(colMeans(cbind(1,x))[-1])
         #x_col_mean = as.numeric(sprintf("%.5f",x_col_mean_old)) # 5 casas decimais
         #result = tryCatch({

         integrandoF_ = function(z){
           FVP_j(y = as.matrix(y),
                 z_values = z, explanatory_variables = as.matrix(x),
                 W = w )[k,2]}


          # integrando vetorizado
          integrandoF4_ = Vectorize(integrandoF_, SIMPLIFY = TRUE, USE.NAMES = TRUE);

          # constant de integração para W fixado
          constant_integracao_aux = cubature::hcubature(integrandoF4_, -Inf, Inf, tol= 1e-10, fDim=1, maxEval=10000)$integral;

          #res <- constant_integracao_aux
          #}, error = function(e) {
          #res <- 0
          #})
        }

  #const_int

  # vetor com constants de integração
  constant_integracao_ = as.numeric(const_int)
  num_loop = 1 #ordem do loop

  #constants estão aqui
  #constant_integracao_ # teste via print

  ####################################### Cálculo das integrais DKL
  # lembrar que é um integral dupla

  #####################Calcular FVPP com todos os dados
  #Com  pacote cubature e adaptIntegrate(...)

  ##### Obtendo constant de integração
  # integrando para Y (respostas), X (variáveis explicativas) e W fixados - z livre
  integrandoF = function(z){
    FVP(y = as.matrix(y), z_values = z,
        explanatory_variables = as.matrix(x), W = w )[1,2]};

  # integrando vetorizado
  integrandoF4 = Vectorize(integrandoF, SIMPLIFY = TRUE, USE.NAMES = TRUE)

  # constant de integração para W fixado
  constant_integracao = cubature::hcubature(integrandoF4, -Inf, Inf, tol= 1e-10, fDim=1, maxEval=10000)$integral

  #constant_integracao # teste via print

  ############################################ fim

  # carregar constants de integração de FVPP e FVPPjs

  ###### Partindo para o cálculo do DKL por observação

  #valores dos hiperparâmetros (número de observações)
  number_of_observations = 1:length(y)


  DKL_num <- foreach::foreach(k=number_of_observations,
                              #.export = c("integrandoDKL", "FVP_j", "y", "x"),
                              .combine='rbind', .packages = c("diagL1", "cubature")) %dopar% {

              #uso as médias das variáveis explicativas como vetor W.
              #result = tryCatch({

              # Cálculo da DKL
              integrandoDKL = function(z){
                FVP(y = as.matrix(y), z_values = z,
                    explanatory_variables = as.matrix(x),
                    W = w )[1,2]*log(FVP(y = as.matrix(y),
                    z_values = z, explanatory_variables = as.matrix(x),
                    W = w )[1,2]/FVP_j(y = as.matrix(y),
                    z_values = z, explanatory_variables = as.matrix(x),
                                                          W = w )[k,2])
                                }

                  # integrando vetorizado
                  integrandoDKL_V = Vectorize(integrandoDKL, SIMPLIFY = TRUE, USE.NAMES = TRUE);

                  # constant de integração para W fixado

                  #-1000, 1000
                  resultado_integracao_DKL = cubature::hcubature(integrandoDKL_V, -Inf, Inf, tol= 1e-5, fDim=1, maxEval=10000)$integral;

                  DKL_final = log(constant_integracao_[k]/constant_integracao) +(constant_integracao^-1)*resultado_integracao_DKL

                  #res <- DKL_final
                  #}, error = function(e) {
                  #res <- 0
                  #})
              }

  #as.numeric(DKL_num)
  #as.numeric(DKL_num)[as.numeric(DKL_num)<0]
  #as.numeric(const_int)

  result = as.numeric(DKL_num)

  # Stop cluster
  parallel::stopCluster(cluster)
  options(warn = 0)
  #####END

  return(result)
}






