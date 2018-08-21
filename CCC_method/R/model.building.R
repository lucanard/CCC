#' @title Build your own model for the CCC_method
#' @description the function builds the statistical models to predict the Y dependent variables of theCCC_method  from a chemical standards dataset.   
#' @param X1Y the dataset built by the dataset building function
#' @param ny the Y dependent variable you want to build the model on (select a value between 1 and 9). Default option includes all the Y variables used in the CCC_method
#' @param models the statistical model you want to use to build the model. The options inclued logistic ("logistic"), lasso ("lasso"), ridge ("ridge"), pls ("pls"), and penalized-pls ("b_ppls"). The default option builds the models according to the CCC_method.  
#' @usage model.building(X1Y, ny, models)
#' @return a list of statistical models to perform the CCC approach 
#' @import glmnet
#' @import ppls
#' @export "model.building"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' STD_RP <- read.csv(system.file("extdata", "STD_RP.csv", package = "CCC"), row.names = 1, stringsAsFactors = FALSE)
#' X1Y <- dataset.building(STD_RP)
#' models <- model.building(X1Y)
model.building <- function(X1Y, ny = c(1:9), models = NULL) {
  options(warn=-1)
modelli <- list(length=length(ny))
if (is.null(models)) {
for (i in ny) {
  if (i == 1) {model <- "logistic"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 2) {model <- "pls"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 6) {model <- "pls"
  modelli [[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 4) {model <- "logistic"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 5) {model <- "logistic"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 8) {model <- "ridge"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 9) {model <- "logistic"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 3) {model <- "pls"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  if (i == 7) {model <- "lasso"
  modelli[[i]] <- building(X1Y, ny = i, model = model)
  }
  }
} else {model <- models
for (i in ny) {
modelli[[i]] <- building(X1Y, ny = i, model = model[i])
}
}
return(modelli)
}