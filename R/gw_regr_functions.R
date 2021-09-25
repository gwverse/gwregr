# from ff_play_v3_gwCore_functions.R
# gw_do_local_lm - undertakes the local weighted regression in an `apply` function. 
# gw_get_lm_eval - returns an evaluation function to evaluate the GWR results 
# from gw_verse_v2.RMD
# gw_single_bw_gwr - creates a function to evaluate a single GWR bandwidth
# gw_regr - returns a gwr function once thebandwidth has been specified

# see https://ourcodingclub.github.io/tutorials/writing-r-package/
# library(devtools)
# setwd("/Users/geoaco/Desktop/my_docs_mac/leeds_work/research/gw_verse/gw_verse_R/gw_/R)
# load_all(".")
# setwd("/Users/geoaco/Desktop/my_docs_mac/leeds_work/research/gw_verse/gw_verse_R/gw_")
# library(roxygen2); # Read in the roxygen2 R package
# roxygenise();


#' Defines an evaluation function to be used in determining the GWR bandwidth
#' 
#' This returns an evaluation function for use in specifying GWR bandwidths for a given regression model formula. The options are "AIC" for a corrected AIC evaluation, or "CV" for a leave one out cross validation. It uses a matrix approach to assess the regression specified in the formula.   
#' @param eval Either "AIC" or "CV". 
#' @return The evaluation function.
#' @examples
#' library(sf)
#' data(georgia)
#' # define a distance matrix, a location and an adaptive bandwidth
#' dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
#' # create the nearby function - see the help for `gw_get_nearby` in the gw_package
#' nearby_func = gw_get_nearby(adaptive = TRUE)
#' # create the weighting function - see the help for `gw_get_weight` in the gw_ package
#' weight_func = gw_get_weight(kernel = "bisquare", adaptive = TRUE)
#' # adaptive r fixed bandwidth
#' adaptive = TRUE
#' # specify a bandwidth
#' bw = 50
#' # create an index of observations
#' indexMat  = matrix(1:nrow(georgia), ncol = 1) 	
#' # determine nearby locations for all observation (list for adaptive, matrix for fixed bandwidth) 
#' nearbyMat = apply(indexMat, 1, function(x) 
#'   nearby_func(x, dist_mat, bw))
#' # apply the weight function to the nearby locations 
#' weightMat = apply(indexMat, 1, function(x) 
#'   gw_do_weight(x, bw, nearbyMat, dist_mat, weight_func))
#' # define and run the evaluation function
#' eval_func = gw_get_lm_eval(eval = "AIC")
#' as.vector(eval_func(st_drop_geometry(georgia) , formula, indexMat, weightMat))
#' @export
gw_get_lm_eval = function(eval){
  if (eval == "AIC") {
    return(function(data, formula, indexMat, weightMat){
      Y = as.matrix(data[,all.vars(formula)[1]])
      X = as.matrix(data[,all.vars(formula)[-1]])
      X = cbind(1,X)
      nobs = nrow(data)
      # need to remember to transpose the apply result!
      S = t(apply(indexMat, 1, function(x) 	{
        X[x,] %*% solve(t(X*weightMat[,x])%*%X)%*%{t(X*weightMat[,x])} 
      } ) )
      tr.S = sum(diag(S))
      RSS = t(Y)%*%t(diag(nobs)-S)%*%(diag(nobs)-S)%*%Y
      sigma.hat = RSS/nobs
      AICc = nobs*log(sigma.hat) + nobs*log(2*pi) + nobs*((nobs + tr.S) / (nobs - 2 - tr.S))
      AICc
    }		
    )
  }
  if (eval == "CV") {
    return(function(data, formula, indexMat, weightMat){
      Y = as.matrix(data[,all.vars(formula)[1]])
      X = as.matrix(data[,all.vars(formula)[-1]])
      X = cbind(1,X)
      nobs = nrow(data)
      CV = apply(indexMat, 1, function(x)  {
        weightMat[,x][x] = 0
        Y[x] - (X[x,] %*% solve(t(X*weightMat[,x])%*%X)%*%{t(X*weightMat[,x])%*%Y})
      } )
      CV = t(CV) %*% CV
      CV	
    }
    )
  }
}


#' Evaluate a single GWR bandwidth
#' 
#' This returns a function that evaluates a geographically weighted regression under a given bandwidth. It uses functions from the core gw_verse package (gw_) to specify local observation and selection and weighting functions (get_neraby_func, gw_get_weight) and requires an evaluation function to be specified (see get_lm_eval_fun in this package). The returned function can be used evaluate different bandwidths, for a given regression model formula.   
#' @param sf_data A point or polygon spatial dataset in sf format, containing the attributes to modelled. 
#' @param adative A logical value TRUE or FALSE to indicate whether an adaptive or fixed bandwidth distance is being used.
#' @param kernel The type of distance weighting to be used: one of "bisquare", "gaussian", "exponential", "tricube" or "boxcar".
#' @param eval An evaluation method for the local model, either "AIC" for a corrected AIC evaluation, or "CV" for a leave one out corss validation.  
#' @return The evaluation measure, to be minimised in bandwidth evaluation.
#' @examples
#' library(sf)
#' # load data and define a model formula
#' data(georgia)
#' formula = as.formula(MedInc ~ PctBach + PctEld)
#' ## 1. A single ADAPTIVE bandwidth
#' # create the function
#' gwr_bw_func = gw_single_bw_gwr(georgia,adaptive=TRUE,kernel="bisquare",eval="AIC")
#' # apply it
#' gwr_bw_func(bw = 100, formula)	
#' ## 2. Determining an optimal bandwidth  with `optimse` 
#' optimise(gwr_bw_func,c(10,nrow(georgia)),formula=formula,maximum=FALSE)
#' ## 3. Determining an adaptive bandwidth with a linear search
#' # slower but guarantees local minima are not returned
#' bwsa = 10:nrow(georgia)
#' # apply the function to all adaptive bandwidths
#' res = sapply(bwsa, function(x) gwr_bw_func(x, formula))
#' # find the minimum
#' bwsa[which.min(res)]
#' ## 4. Determining an optimal FIXED bandwidth 
#' gwr_bw_func = gw_single_bw_gwr(georgia,adaptive=FALSE,kernel="tricube",eval="CV")
#' gwr_bw_func(bw =130000, formula)
#' dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
#' bwsf = seq(50000, (ceiling(max(dist_mat)/5000) * 5000), 5000)
#' # 4a) by optimising
#' optimise(gwr_bw_func,bwsf,formula=formula,maximum=FALSE)
#' # 4b) with a linear search
#' res = sapply(bwsf, function(x) gwr_bw_func(x, formula))
#' bwsf[which.min(res)]
#' @export
gw_single_bw_gwr= function(sf_data, adaptive, kernel, eval) {
  ## 1. Input data 
  dist_mat = as.matrix(dist(st_coordinates(st_centroid(sf_data)), upper = T, diag = T))
  df = st_drop_geometry(sf_data)
  ## 2. Core - from gw_ 
  nearby_func = gw_get_nearby(adaptive)
  weight_func = gw_get_weight(kernel, adaptive)
  ## 3. Application - gw_regr 
  eval_func = gw_get_lm_eval(eval)
  # Tidy un-needed objects
  rm(sf_data)
  # define returned function
  function(bw, formula) {
    indexMat  = matrix(1:nrow(df), ncol = 1) 	
    nearbyMat = apply(indexMat, 1, function(x) 
      nearby_func(x, dist_mat, bw))
    weightMat = apply(indexMat, 1, function(x) 
      gw_do_weight(x, bw, nearbyMat, dist_mat, weight_func))
    res = as.vector(eval_func(df, formula, indexMat, weightMat))
    return(res)
  }
}

#' Undertakes a single weighted regression 
#'
#' Returns the local GWR coefficients for a given bandwidth. This is usually run after the optimal bandwidth has been determined, within the gw_regr function and not as a standalone function. 
#' @param w A vector of weights with length equal to the number of observations in data. 
#' @param formula A formula to be used in the regression, with terms contained in the variables names of the input data.
#' @param data A flat data table in data.frame, matrix or tibble format
#' @return A vector of coefficient estimates in the order Intercept, Variable 1, Variable 2, etc
#' @examples
#' library(sf)
#' # load data and define a model formula
#' data(georgia)
#' formula = as.formula(MedInc ~ PctBach + PctEld)
#' ##  1. how the function works
#' # define a distance matrix, a location and an adaptive bandwidth
#' dist_mat = as.matrix(dist(st_coordinates(st_centroid(georgia)), upper = T, diag = T))
#' obs_index = 70
#' bw = 30
#' # create the nearby function - see the help for `gw_get_nearby` in the gw_ package
#' nearby_func = gw_get_nearby(adaptive = TRUE)
#' # apply to get an index of locations and get a vector of distances
#' index = nearby_func(obs_index, dist_mat, bw)
#' dists = dist_mat[obs_index,index]
#' # create the weighting function and weight the nearby locations 
#' #  - see the help for `gw_get_weight` in the gw_ package
#' weight_func = gw_get_weight(kernel = "bisquare", adaptive = TRUE)
#' # apply the weight function to the distances and extend to a vector of length observations
#' w = weight_func(bw, dists)
#' w_vec = rep(0, nrow(georgia))
#' w_vec[index] = w
#' # extract the coefficients for this location
#' gw_do_local_lm(w_vec, formula, georgia)
#' ## 2. Using the function operationally
#' 
#' @export
gw_do_local_lm = function(w, formula, data){
  data$w = w
  res = lm(formula = formula, data = data, weights = w)
  return(as.vector(coef(res)))
}

#' Undertakes a GWR  
#'
#' Returns a function to undertake GWR once the optimal bandwidth has been defined. 
#' @param formula A formula to be used in the regression, with terms contained in the variables names of the input spatial data.
#' @param sf_data A point or polygon spatial dataset in sf format, containing the attributes to modelled. 
#' @param adative A logical value TRUE or FALSE to indicate whether an adaptive or fixed bandwidth distance is being used.
#' @param kernel The type of distance weighting to be used: one of "bisquare", "gaussian", "exponential", "tricube" or "boxcar".
#' @param bw The bandwidth (fixed or adaptive), either manually specified or preferably as determined through an evaluation of bandwidths using gw_single_bw_gwr.
#' @return A matrix of coefficient estimates, with each row representing the corresponding location in the input data, and each column the coefficient estimates in the order Intercept, Variable 1, Variable 2, etc as specified in the formula. 
#' @examples
#' library(sf)
#' library(tmap)
#' #load data and define a model formula
#' data(georgia)
#' formula = as.formula(MedInc ~ PctBach + PctEld)
#' # define the bandwidth function
#' gwr_bw_func = gw_single_bw_gwr(georgia,adaptive=TRUE,kernel="bisquare",eval="AIC")
#' # define a vector of adaptive bandwidths
#' bwsa = 10:nrow(georgia)
#' # apply the bandwidth function to all bandwidths
#' res = sapply(bwsa, function(x) gwr_bw_func(x, formula))
#' # return the best result
#' bw = bwsa[which.min(res)]
#' # define the final GWR model and run it
#' gw_regr_final = gw_regr (formula, georgia, adaptive = TRUE, kernel = "bisquare", bw = bw)
#' coef_mat = gw_regr_final(formula)
#' # examine the result
#' head(coef_mat)
#' summary(coef_mat)
#' # rename and then map the coefficient estimates
#' colnames(coef_mat) = c("Intercept", paste0(all.vars(formula)[-1], "CE"))
#' georgia = cbind(georgia, coef_mat)
#' tm_shape(georgia)+tm_polygons(c("Intercept","PctBachCE","PctEldCE")) +
#' tm_layout(legend.position = c("right","top"), frame = F)
#' @export
gw_regr  = function(formula, sf_data, adaptive, kernel, bw) {
  ## 1. Input data related
  dist_mat = as.matrix(dist(st_coordinates(st_centroid(sf_data)), upper = T, diag = T))
  df = st_drop_geometry(sf_data)
  ## 2 core - ie gw_ related
  nearby_func = gw_get_nearby(adaptive)
  weight_func = gw_get_weight(kernel, adaptive)
  ## 3. application - gw_regr related
  # gw_do_local_lm in the returned function below
  # Tidy un-needed objects
  rm(sf_data)
  function(formula) {
    indexMat  = matrix(1:nrow(df), ncol = 1) 	
    nearbyMat = apply(indexMat, 1, function(x) 
      nearby_func(x, dist_mat, bw))
    weightMat = apply(indexMat, 1, function(x) 
      gw_do_weight(x, bw, nearbyMat, dist_mat, weight_func))
    coefMat = t(apply(weightMat, 2, function(x) 
      gw_do_local_lm(x, formula, df)))
    return(coefMat)
  }
}
