library(shiny)
library(mixlink)

ui <- shinyUI(fluidPage(
	titlePanel("MixLink Binomial (J=2) vs. Binomial"),
	sidebarLayout(
		sidebarPanel(
			sliderInput("m", "m:", min = 1, max = 50, value = 20),
			sliderInput("mean", "mean:", min = 0.01, max = 0.99, value = 0.5),
			sliderInput("Pi1", "Pi1:", min = 0.01, max = 0.99, value = 0.5),
			sliderInput("kappa", "kappa:", min = 0.01, max = 5, value = 5)
		),
		mainPanel(
			plotOutput("mixlinkPlot", height = 300),
			plotOutput("binPlot", height = 300)
		)
	)
))

server <- shinyServer(function(input, output) {
	output$mixlinkPlot <- renderPlot({
		m <- input$m
		Pi <- c(input$Pi1, 1-input$Pi1)
		ff <- d.mixlink.binom(0:m, m, mean = input$mean, Pi = Pi, kappa = input$kappa)
		barplot(ff, names.arg = 0:m, col = "blue")
		title("MixLinkJ2(m, mean, Pi, kappa) Density")
		box()
	})
	output$binPlot <- renderPlot({
		m <- input$m
		ff <- dbinom(0:m, prob = input$mean, size = m)
		barplot(ff, names.arg = 0:m, col = "green")
		title("Binomial(m, mean) Density")
		box()
	})
})

shinyApp(ui = ui, server = server)
