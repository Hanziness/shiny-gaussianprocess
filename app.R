# (c) Imre Gera, 2019
# www.stud.u-szeged.hu/Gera.Imre.1

library(shiny)
library(ggplot2)
library(data.table)
library(dplyr)
library(e1071)

# Kernel functions

kernel_randomPlane <- function(x, y, ...) {
    return(x %*% t(y))
}

kernel_brown <- function(x, y, ...) {
    len = length(x)
    return(pmin(x, y))
}

kernel_brownBridge <- function(x, y, ...) {
    # x[, 1] is the original x vector & y is just x transposed (so the original vector is also y[1, ])
    return(pmin(x, y) - (x[, 1] %*% t(y[1, ])))
}

kernel_squareExp1D <- function(x, y, alpha = 1, ...) {
    # this actually uses the 2-norm, but we can avoid it since it has no effect
    # on 1-dimensional vectors
    return(exp(-alpha * (x - y)^2))
}

kernel_ornstein_uhlenbeak <- function(x, y, alpha = 1, ...) {
    return(exp(-alpha * abs(x - y)))
}

kernel_periodic <- function(x, y, alpha = 1, beta = 1) {
    return(exp( -alpha * sin(beta * pi * (x - y))^2 ))
}

# Creates trajectory data
createTrajectoryData <- function (input, seed = NULL) {
    # set seed if it was given
    if (is.null(seed)) {
        # if parameter is NULL, use "rndseed" from the input object
        set.seed(input$rndseed)
    } else {
        set.seed(seed)
    }

    # sampling points of the process
    x = seq(from = input$xrange[1], to = input$xrange[2], length.out = input$samplingFreq)
    n = input$samplingFreq
    C = matrix(0, nrow = n, ncol = n)
    covF = get(input$kernelF) # covF will become a function this way

    # Create the covariance matrix (using vectorized inputs)
    if (input$kernelF %in% c("kernel_randomPlane")) {
        # the random plane kernel achieves the effect by multiplying two vectors
        # so we don't convert x into a matrix here!
        C = covF(x, x)
    } else {
        # the rest can deal with two matrices -> repeat x n times (becomes a matrix)
        # -> the functions will create the covariance matrix from mX and mX'
        mX = matrix(rep(x, n), nrow = n)
        C = covF(mX, t(mX), alpha = input$alpha, beta = input$beta)
    }

    # sample from the Gaussian process
    u = matrix(rnorm(n, mean = 0, sd = 1), ncol = 1)      # u ~ N(0, 1)
    dRes = svd(C)                                         # C = USV
    y = dRes$u %*% diag(sqrt(dRes$d)) %*% u               # U S^0.5 u

    # Prepare data for plotting
    plotData = data.frame(x = x, y = y)
    return(plotData)
}

# Function to simulate a Gauissan process' time to reach a certain state
simulate_ReachTime <- function (input, numRuns = 200, absCeilFloor = 1) {
    results = data.table(matrix(numeric(), nrow = numRuns * 2, ncol = 3, dimnames = list(NULL, c("point", "reachTime", "runId"))))
    for (i in 1:numRuns) {
        thisRun = createTrajectoryData(input, seed = input$rndseed + (i - 1))
        firstNegative = thisRun[head(which(thisRun$y <= (-absCeilFloor)), n = 1), ]$x
        firstPositive = thisRun[head(which(thisRun$y >= (absCeilFloor)), n = 1), ]$x

        if (length(firstNegative) == 0) {
            firstNegative = NA
        }

        if (length(firstPositive) == 0) {
            firstPositive = NA
        }

        results[i * 2 - 1] <- list(-absCeilFloor, firstNegative, i)
        results[i * 2] <- list(absCeilFloor, firstPositive, i)
    }

    return(results)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    navbarPage("Gaussian Process Simulation",
        tabPanel("Trajectory",
        # Sidebar with a slider input for number of bins
        sidebarLayout(
            sidebarPanel(
                selectInput(
                    "kernelF",
                    "Kernel function to use",
                    choices = c(
                        "Random plane" = "kernel_randomPlane",
                        "Brownian motion" = "kernel_brown",
                        "Brownian bridge" = "kernel_brownBridge",
                        "Squared exponential" = "kernel_squareExp1D",
                        "Ornstein-Uhlenbeck" = "kernel_ornstein_uhlenbeak",
                        "Periodic" = "kernel_periodic"
                    )
                ),
                conditionalPanel(
                    condition = "['kernel_brown', 'kernel_brownBridge', 'kernel_ornstein_uhlenbeak'].indexOf(input.kernelF) >= 0",
                    helpText(
                        "The Brownian motion and Ornstein-Uhlenbeck kernels should",
                        "be used with the x axis minimum set to 0. The Brownian bridge",
                        "kernel needs the interval to be set to [0, 1]."
                    )
                ),
                sliderInput(
                    "xrange",
                    "Range on x axis",
                    min = -2,
                    max = 30,
                    step = 0.1,
                    value = c(0, 1)
                ),
                numericInput(
                    "samplingFreq",
                    "Number of samples to take",
                    min = 2,
                    max = 500,
                    step = 1,
                    value = 30
                ),
                numericInput(
                    "rndseed",
                    "Random seed",
                    min = 0,
                    step = 1,
                    value = 1997
                ),
                checkboxInput(
                    "plot_addLines",
                    "Connect points with a line?",
                    value = FALSE
                ),
                conditionalPanel(
                    condition = "['kernel_squareExp1D', 'kernel_ornstein_uhlenbeak', 'kernel_periodic'].indexOf(input.kernelF) >= 0",
                    numericInput(
                        "alpha",
                        "Alpha",
                        min = 0.001,
                        step = 0.1,
                        value = 1
                    )
                ),
                conditionalPanel(
                    condition = "input.kernelF === 'kernel_periodic'",
                    numericInput(
                        "beta",
                        "Beta",
                        min = 0.001,
                        step = 0.1,
                        value = 1
                    )
                )
            ),


            # Show a plot of the generated trajectory
            mainPanel(
               plotOutput("gaussianTrajectory")
            )
        )),
        tabPanel("Time to reach state",
                 sidebarLayout(
                     sidebarPanel = sidebarPanel(
                         sliderInput(
                             "numSimulations",
                             "Number of simulations",
                             min = 1,
                             max = 1500,
                             step = 1,
                             value = 300
                         ),
                         sliderInput(
                             "targetStateAbs",
                             "Target state (absolute)",
                             min = 0,
                             max = 10,
                             step = 0.1,
                             value = 1
                         )
                     ),
                    mainPanel = mainPanel(
                        tabsetPanel(type = "tabs",
                            tabPanel("Mean time", plotOutput("gaussianReachTimePlot")),
                            tabPanel("Histograms",
                                     plotOutput("gaussianReachPositiveHist"),
                                     plotOutput("gaussianReachNegativeHist"))
                        )
                    )
                 )),
        tabPanel("Brownian bridge",
                 sidebarLayout(
                     sidebarPanel = sidebarPanel(
                         numericInput(
                             "bridgeseed",
                             "Seed for the Brownian bridge",
                             min = 0,
                             step = 1,
                             value = 1997
                         ),
                         sliderInput(
                             "bridgefreq",
                             "Frequency",
                             min = 2,
                             max = 2000,
                             value = 1000
                         ),
                         sliderInput(
                             "bridgeSimulations",
                             "Number of simulations",
                             min = 1,
                             max = 2000,
                             value = 500
                         )
                     ),
                     mainPanel = mainPanel(
                        tabsetPanel(type = "tabs",
                            tabPanel("Trajectory", plotOutput("brownianBridgeTrajectory")),
                            tabPanel("Reach", plotOutput("brownianBridgeReach"))
                        )
                     )
                 ))
        )
)



# Define server logic
server <- function(input, output) {

    plotRes <- 120

    # Only runs when the reach time-related input change
    simData_reach = reactive({
        simulate_ReachTime(input = input, numRuns = input$numSimulations, absCeilFloor = input$targetStateAbs)
    })

    # Simulate a Gaussian process' trajectory.
    output$gaussianTrajectory <- renderPlot({
        plotData = createTrajectoryData(input)

        plt = ggplot(plotData, aes(x = x, y = y)) +
            theme_minimal() +
            geom_point(size = 1.2) +
            labs(x = "t") +
            coord_fixed(xlim = c(input$xrange[1], input$xrange[2]))


        if (input$plot_addLines) {
            plt = plt + geom_line()
        }

        print(plt)

    }, res = plotRes)

    # Creates a bar plot showing on average when the process reached the
    # positive and negative boundaries.
    output$gaussianReachTimePlot <- renderPlot({
        finalData = simData_reach() %>%
            mutate(gr = cut_interval(point, n = 2)) %>%
            group_by(gr) %>%
            summarize(m = mean(reachTime, na.rm = TRUE))

        plt = ggplot(finalData, aes(x = gr, y = m, group = gr, fill = gr)) +
            theme_minimal() +
            geom_bar(stat = "identity") +
            geom_text(aes(label=m), vjust=-0.3, size=3.5) +
            scale_x_discrete(labels = c(as.character(-input$targetStateAbs), as.character(input$targetStateAbs))) +
            guides(fill = FALSE) +
            labs(x = "Boundary", y = "Average time to reach boundary")

        print(plt)
    }, res = plotRes)

    # Creates a histogram and overlays a density curve on when the process
    # reached the positive boundary.
    output$gaussianReachPositiveHist <- renderPlot({
        finalData = simData_reach() %>%
            na.omit() %>%
            filter(point >= 0)

        ggplot(finalData, aes(x = reachTime, y = ..density..)) +
            geom_histogram(binwidth = 0.1, fill = "#00796B", alpha = 0.75) +
            geom_density(aes(y=..density..), colour = "#E91E63", fill = "#E91E63", size = 1.0, alpha = 0.2) +
            labs(x = "Time", y = "Density", title = "Time to reach positive boundary", subtitle = paste0("Data points: ", nrow(finalData))) +
            theme_minimal()
    }, res = plotRes)

    # Same as the previous one, but with the negative boundary.
    output$gaussianReachNegativeHist <- renderPlot({
        finalData = simData_reach() %>%
            na.omit() %>%
            filter(point <= 0)

        ggplot(finalData, aes(x = reachTime, y = ..density..)) +
            geom_histogram(binwidth = 0.1, fill = "#00796B", alpha = 0.75) +
            geom_density(aes(y=..density..), colour = "#1976D2", fill = "#1976D2", size = 1.0, alpha = 0.2) +
            labs(x = "Time", y = "Density", title = "Time to reach negative boundary", subtitle = paste0("Data points: ", nrow(finalData))) +
            theme_minimal()
    }, res = plotRes)

    # Displays a trajectory of a Brownian bridge using rbridge from package e1071
    output$brownianBridgeTrajectory <- renderPlot({
        set.seed(input$bridgeseed)
        plotFrame = data.table(x = seq(from = 0, to = 1, length.out = input$bridgefreq), y = as.numeric(rbridge(end = 1, frequency = input$bridgefreq)))

        plt <- ggplot(plotFrame, aes(x = x, y = y)) +
            geom_hline(yintercept = 0.0, size = 0.7, colour = "#37474F", alpha = 0.8) +
            geom_line(size = 0.8) +
            geom_hline(yintercept = min(plotFrame$y), colour = "#0D47A1", size = 0.4, alpha = 0.6) +
            geom_hline(yintercept = max(plotFrame$y), colour = "#E65100", size = 0.4, alpha = 0.6) +
            scale_x_continuous(limits = c(0, 1)) +
            labs(x = "t", title = "Trajectory of a Brownian bridge") +
            theme_minimal()

        print(plt)
    }, res = plotRes)

    # histogram of the maximum reaches for the brownian bridge
    output$brownianBridgeReach <- renderPlot({
        runs = seq(from = 1, to = input$bridgeSimulations, by = 1)

        y = sapply(runs, function (x) {
            set.seed(input$bridgeseed + x)
            max(abs(rbridge(frequency = input$bridgefreq)))
        })

        plotFrame <- data.table(x = runs, y = y)

        plt <- ggplot(plotFrame, aes(x = y, y = ..density..)) +
            geom_histogram(fill = "#00796B", alpha = 0.75, bins = 60) +
            geom_density(fill = "#EF6C00", colour = "#EF6C00", alpha = 0.2, size = 1.0) +
            labs(x = "Maximum absolute reach", y = "Density", title = "Maximum reaches of the Brownian bridge") +
            theme_minimal()

        print(plt)

    }, res = plotRes)
}

# Run the application
shinyApp(ui = ui, server = server)
