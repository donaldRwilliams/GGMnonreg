#' @importFrom stats na.omit quantile cov2cor pnorm qnorm sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
symm_mat <- function (x) {
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  x
}

z2r <- function (z) {
  (exp(2 * z) - 1)/(1 + exp(2 * z))
}

fisher_z <- function(rho){
  .5 * log(( 1 + rho )/ ( 1 - rho ))
}

power_z <- function(r, n, c = 0,
                    type = "pearson",
                    compare = FALSE,
                    alpha = 0.05){

  # abs r
  r <- abs(r)

  # fisher z
  z <- fisher_r_to_z(r)

  if(type == "pearson"){
    # variance of z
    var_z <- 1 / (n - c - 3)

  } else if(type == "spearman"){

    var_z <- (1 + r^2/2) / (n - c - 3)

  } else{

    stop("invalid type (must be pearson or spearman)")
  }


  # differnece ?
  if(compare == TRUE){

    var_z <- var_z * 2

  }

  # z score
  z_score <- z/ sqrt(var_z)

  # quantile
  q <- stats::qnorm(1 - alpha / 2)

  # power
  1 - stats::pnorm(q - z_score)

}


csws_labels <- ifelse(1:35 %in% c(7,10,16,24,29),
                      "Family Support",
                      ifelse(1:35 %in% c(3,12,20,25,35),
                             "Competition",
                             ifelse(1:35 %in% c(1,4,17,21,30),
                                    "Appearence",
                                    ifelse(1:35%in%c(2,8,18,26,31),
                                           "God's Love",
                                           ifelse(1:35 %in% c(13, 19, 22, 27,  33),
                                                  "Academic Competence",
                                                  ifelse(1:35 %in% c(5, 11, 14, 28, 34),
                                                         "Virtue", "Approval From Others"))))))

tas_labels <- ifelse(1:20 %in% c(1,3,6,7,9,13,14),
                     "Difficulty\nIdentifying Feelings",
                     ifelse(1:20 %in% c(2,4,11,12,17),
                            "Difficulty\nDescribing Feelings",
                            "Externally\nOriented Feelings"))

iri_labels <- ifelse(1:28 %in% c(3, 8, 11, 15, 21, 25, 28),
                     "Perspective Taking",
                     ifelse(1:28 %in% c(2, 4, 9, 14, 18, 20, 22),
                            "Empathic Concern",
                            ifelse(1:28 %in% c(1, 5, 7, 12, 16, 23, 26), "Fantasy",
                                   "Personal Distress")))

rsa_labels <- ifelse(1:33 %in% c(1, 4, 5, 32),
                     "Planned Future",
                     ifelse(1:33 %in% c(2, 11, 17, 25, 31, 33),
                            "Perception of Self",
                            ifelse(1:33 %in% c(3, 7, 13, 16, 24, 29),
                                   "Family Cohesion",
                                   ifelse(1:33  %in% c(6, 9, 10, 12, 15, 19, 27),
                                          "Social Resources",
                                          ifelse(1:33 %in% c(8, 14, 18, 21, 22, 26),
                                                 "Social Competence", "Structured Style")))))


globalVariables(c("p", "n"))

