library(shiny)
library(ggplot2)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "input_tabs",
        tabPanel("Current Tephrostratigraphy",
                 numericInput("total_sites", "Total Sites", value = 30),
                 numericInput("total_ash_layers", "Total Ash Layers", value = 20)
        ),
        tabPanel("Proportions of Current Layers",
                 numericInput("cat_a_proportion", "Cat A Proportion", value = 0.25),
                 numericInput("cat_b_proportion", "Cat B Proportion", value = 0.65),
                 numericInput("cat_c_proportion", "Cat C Proportion", value = 0.1)
        ),
        tabPanel("Probability Options",
                 numericInput("cat_a_prob", "Cat A Probability", value = 0.1),
                 numericInput("cat_b_prob", "Cat B Probability", value = 0.6),
                 numericInput("cat_c_prob", "Cat C Probability", value = 0.9)
        )
      ),
      numericInput("num_simulations", "Number of Simulations", value = 100),
      actionButton("run_simulation", "Run Simulation"),
      
      # Explanatory text section
      wellPanel(
        h4("Parameter Explanations"),
        HTML("<p><b>Total Sites:</b> The total number of sites in the landscape.</p>"),
        HTML("<p><b>Total Ash Layers:</b> The total number of ash layers in the landscape.</p>"),
        HTML("<p><b>Cat A Proportion:</b> The proportion of ash layers classified as Category A.</p>"),
        HTML("<p><b>Cat B Proportion:</b> The proportion of ash layers classified as Category B.</p>"),
        HTML("<p><b>Cat C Proportion:</b> The proportion of ash layers classified as Category C.</p>"),
        HTML("<p><b>Cat A Probability:</b> The probability of encountering a Category A layer at a new site.</p>"),
        HTML("<p><b>Cat B Probability:</b> The probability of encountering a Category B layer at a new site.</p>"),
        HTML("<p><b>Cat C Probability:</b> The probability of encountering a Category C layer at a new site.</p>"),
        HTML("<p><b>Number of Simulations:</b> The number of simulations to run.</p>")
      )
    ),
    mainPanel(
      plotOutput("hist_plot"),
      plotOutput("likelihood_plot"),
      verbatimTextOutput("results")
    )
  )
)



server <- function(input, output) {
  observeEvent(input$run_simulation, {
    # Set parameters
    total_sites <- input$total_sites
    total_ash_layers <- input$total_ash_layers
    cat_a_proportion <- input$cat_a_proportion
    cat_b_proportion <- input$cat_b_proportion
    cat_c_proportion <- input$cat_c_proportion
    num_simulations <- input$num_simulations
    cat_a_prob <- input$cat_a_prob
    cat_b_prob <- input$cat_b_prob
    cat_c_prob <- input$cat_c_prob
    
    # Create empty vectors to store results
    site_layers <- matrix(0, nrow = num_simulations, ncol = total_sites)
    average_layers <- numeric(num_simulations)
    likelihood_cat_a_layer <- numeric(num_simulations)
    likelihood_cat_b_layer <- numeric(num_simulations)
    likelihood_cat_c_layer <- numeric(num_simulations)
    
    # Repeat the simulation
    for (i in 1:num_simulations) {
      # Generate the distribution of ash layers
      num_cat_a_layers <- round(cat_a_proportion * total_ash_layers)
      num_cat_b_layers <- round(cat_b_proportion * total_ash_layers)
      num_cat_c_layers <- round(cat_c_proportion * total_ash_layers)
      
      # Calculate the probabilities vector
      probabilities <- numeric(total_ash_layers)
      probabilities[1:num_cat_a_layers] <- cat_a_prob / num_cat_a_layers
      probabilities[(num_cat_a_layers + 1):(num_cat_a_layers + num_cat_b_layers)] <- cat_b_prob / num_cat_b_layers
      probabilities[(num_cat_a_layers + num_cat_b_layers + 1):total_ash_layers] <- cat_c_prob / num_cat_c_layers
      
      # Normalize the probabilities to sum up to 1
      probabilities <- probabilities / sum(probabilities)
      
      # Simulate the sites
      cat_a_site_ash_layers <- sample(1:total_ash_layers, size = total_sites, replace = TRUE, prob = probabilities)
      
      site_layers[i, ] <- cat_a_site_ash_layers
      
      # Calculate the average number of layers in a new site
      average_layers[i] <- length(unique(cat_a_site_ash_layers))
      
      # Calculate the likelihood of encountering cat_a, cat_b, and cat_c layers at a new site
      likelihood_cat_a_layer[i] <- sum(cat_a_site_ash_layers %in% 1:num_cat_a_layers) / total_sites
      likelihood_cat_b_layer[i] <- sum(cat_a_site_ash_layers %in% (num_cat_a_layers + 1):(num_cat_a_layers + num_cat_b_layers)) / total_sites
      likelihood_cat_c_layer[i] <- sum(cat_a_site_ash_layers %in% (num_cat_a_layers + num_cat_b_layers + 1):total_ash_layers) / total_sites
    }
    
    # Calculate overall statistics
    average_layers_mean <- mean(average_layers)
    average_layers_sd <- sd(average_layers)
    confidence_interval <- t.test(average_layers)$conf.int
    
    likelihood_cat_a_layer_mean <- mean(likelihood_cat_a_layer)
    likelihood_cat_a_layer_ci <- t.test(likelihood_cat_a_layer)$conf.int
    
    likelihood_cat_b_layer_mean <- mean(likelihood_cat_b_layer)
    likelihood_cat_b_layer_ci <- t.test(likelihood_cat_b_layer)$conf.int
    
    likelihood_cat_c_layer_mean <- mean(likelihood_cat_c_layer)
    likelihood_cat_c_layer_ci <- t.test(likelihood_cat_c_layer)$conf.int
    
    # Create a data frame for the histogram
    hist_data <- data.frame(Average_Layers = average_layers)
    
    # Create the histogram with confidence interval
    hist_plot <- ggplot(hist_data, aes(x = Average_Layers)) +
      geom_histogram(fill = "darkblue", color = "black", bins = 20) +
      geom_vline(xintercept = average_layers_mean, color = "red", linetype = "dashed", size = 1) +
      geom_text(aes(x = Inf, y = Inf,
                    label = paste("Mean:", round(average_layers_mean, 2))), color = "black", size = 4, hjust = 1, vjust = 1) +
      geom_text(aes(x = Inf, y = Inf,
                    label = paste("Range:", round(confidence_interval[1], 2), "-", round(confidence_interval[2], 2))), color = "black", size = 4, hjust = 1, vjust = 3) +
      labs(title = "Distribution of Average Number of Layers in a New Site",
           x = "Average Number of Layers", y = "Frequency") +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            plot.margin = margin(10, 10, 10, 20, "pt"))
    
    # Create a data frame for the likelihood plot
    likelihood_data <- data.frame(
      Category = c("cat_a", "cat_b", "cat_c"),
      Likelihood = c(likelihood_cat_a_layer_mean, likelihood_cat_b_layer_mean, likelihood_cat_c_layer_mean),
      Lower_CI = c(likelihood_cat_a_layer_ci[1], likelihood_cat_b_layer_ci[1], likelihood_cat_c_layer_ci[1]),
      Upper_CI = c(likelihood_cat_a_layer_ci[2], likelihood_cat_b_layer_ci[2], likelihood_cat_c_layer_ci[2])
    )
    
    # Create the likelihood plot
    likelihood_plot <- ggplot() +
      geom_histogram(data = likelihood_data, aes(x = Category, y = Likelihood, fill = Category),
                     stat = "identity", color = "black", position = "stack", alpha = 0.7) +
      geom_text(data = likelihood_data, aes(x = Category, y = Likelihood,
                                            label = paste(round(Likelihood, 2) * 100, "%")),
                position = position_stack(vjust = 0.5), color = "black") +
      labs(title = "Likelihood of Encountering Tephra Categories",
           x = "Tephra Category", y = "Likelihood") +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            legend.position = "none")
    
    output$hist_plot <- renderPlot({
      hist_plot
    })
    
    output$likelihood_plot <- renderPlot({
      likelihood_plot
    })
    
    output$results <- renderPrint({
      cat("Average number of layers in a new site:", average_layers_mean, "\n")
      cat("95% Confidence Interval for the average number of layers:", confidence_interval[1], "-", confidence_interval[2], "\n")
      cat("Likelihood of encountering a cat_a layer at a new site:", likelihood_cat_a_layer_mean, "\n")
      cat("95% Confidence Interval for the likelihood of encountering a cat_a layer:", likelihood_cat_a_layer_ci[1], "-", likelihood_cat_a_layer_ci[2], "\n")
      cat("Likelihood of encountering a cat_b layer at a new site:", likelihood_cat_b_layer_mean, "\n")
      cat("95% Confidence Interval for the likelihood of encountering a cat_b layer:", likelihood_cat_b_layer_ci[1], "-", likelihood_cat_b_layer_ci[2], "\n")
      cat("Likelihood of encountering a cat_c layer at a new site:", likelihood_cat_c_layer_mean, "\n")
      cat("95% Confidence Interval for the likelihood of encountering a cat_c layer:", likelihood_cat_c_layer_ci[1], "-", likelihood_cat_c_layer_ci[2], "\n")
    })
  })
}

shinyApp(ui, server)
