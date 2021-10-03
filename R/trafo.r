
#' Log-transform child ages
#'
#' @param age Chronological age in years.
#' @param adult.age (Default: 20).
#' @return Ages transformed so that
#' children are below zero and log-transformed
#' and adults are above zero.
#' 
#' @export 
trafo <- function(age,adult.age=20) {
    ifelse(age <= adult.age,
           log(age+1)-log(adult.age+1),
           (age-adult.age)/(adult.age+1))
}

#' Inverse log-transform of child ages
#' 
#' @param x Age estimate where child ages have been log-transformed.
#' @param adult.age (Default: 20)
#' @return Transformed age estimates transformed back to years.
#' 
#' @export 
anti.trafo <- function(x,adult.age=20) {
    ifelse(x<0,
           (1+adult.age)*exp(x)-1,
           (1+adult.age)*x+adult.age)
}
