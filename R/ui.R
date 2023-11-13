#' Create the dashboard UI
#' 
#' @return A shiny dashboard UI
#' @author Jared Andrews
#' 
#' @importFrom shiny fluidRow column numericInput
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody box
#' @importFrom igvShiny igvShinyOutput
#' @importFrom DT DTOutput
#' 
#' @rdname INTERNAL_create_ui
.create_ui <- function() {
    dashboardPage(
        dashboardHeader(title = "ASOspotter"),
        dashboardSidebar(
            numericInput("wing_size", "Wing size", value = 250, min = 1, max = 10000, step = 1)
        ),
        dashboardBody(
            fluidRow(
                column(width = 12,
                    box(
                        title = "Viewer", width = 12, color = "olive",
                        igvShinyOutput("igv")
                    ),
                    box(
                        title = "Variants", width = 12, solidHeader = TRUE, color = "black",
                        DTOutput("variants")
                    )
                )
            )
        )
    )
}