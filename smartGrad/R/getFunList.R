#' getFunList Function
#'
#' This function allows you to
#' @param #love Do you love cats? Defaults to TRUE.
#' @keywords Model
#' @export
#' @examples
#' getFunList()

#--------------------------------------------------------------------->


getFunList<- function(fn = "")
{
  func.list <- c("Rosenbrock_Banana",
                 "Extended_Trigonometric_dim5",
                 "Extended_Freudenstein_Roth",
                 "Perturbed_Quadratic",
                 "Raydan_1",
                 "Raydan_2",
                 "FLETCHCR",
                 "COSINE",
                 "Generalized_Quartic",
                 "Schumer_Steiglitz",
                 "ackley2",
                 "ackley3")

  print(func.list)

}
