library(shiny)
library(mixlink)

ui <- shinyUI(fluidPage(
	titlePanel("MixLink Poisson (J=2) vs. Poisson"),
	sidebarLayout(
		sidebarPanel(
			sliderInput("max", "max:", min = 1, max = 100, value = 20),
			sliderInput("mean", "mean:", min = 0.1, max = 20, value = 10),
			sliderInput("Pi1", "Pi1:", min = 0.01, max = 0.99, value = 0.5),
			sliderInput("kappa", "kappa:", min = 0.01, max = 10, value = 5)
		),
		mainPanel(
			plotOutput("mixlinkPlot", height = 300),
			plotOutput("poisPlot", height = 300)
		)
	)
))

server <- shinyServer(function(input, output) {
	output$mixlinkPlot <- renderPlot({
		m <- input$max
		Pi <- c(input$Pi1, 1-input$Pi1)
		ff <- d.mixlink.pois(0:m, mean = input$mean, Pi = Pi, kappa = input$kappa)
		barplot(ff, names.arg = 0:m, col = "blue")
		title("MixLinkJ2(mean, Pi, kappa) Density")
		box()
	})
	output$poisPlot <- renderPlot({
		m <- input$max
		ff <- dpois(0:m, input$mean)
		barplot(ff, names.arg = 0:m, col = "green")
		title("Poisson(mean) Density")
		box()
	})
})

shinyApp(ui = ui, server = server)
