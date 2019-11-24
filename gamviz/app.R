library(shiny)
ui <- fluidPage(
    titlePanel("Old Faithful Geyser Data"),
    sidebarLayout(
        sidebarPanel(
            
        ),
        
        mainPanel(

        )
    )
)

server <- function(input, output) {
    
}

# Run the application 
shinyApp(ui = ui, server = server)

set.seed(123)
x = runif(500)
mu = sin(2*(4*x-2)) + 2*exp(-(16^2)*((x-.5)^2))
y = rnorm(500, mu, .3)
d = data.frame(x,y)

knots = seq(0, 1, by=.1)
knots = knots[-length(knots)]  # don't need the last value
l = 1
bs = sapply(1:length(knots), function(k) ifelse(x >= knots[k], (x-knots[k])^l, 0))

bs = data.frame(int=1, bs)
d2 = data.frame(x, bs) %>% 
    gather(key=bs, value=bsfunc, -x)
# visualizing basis functions
plot_ly(d) %>% 
    add_markers(~x, ~y, marker=list(color='black', opacity=.25), showlegend=F) %>% 
    add_lines(~x, ~bsfunc, color=~bs, colors='Set3', showlegend=F, data=arrange(d2,x)) %>% # colorscale ignored
    layout()
    
# getting regression coefficients
lmMod = lm(y ~ .-1, data=bs) 
bscoefs = coef(lmMod)
bsScaled = sweep(bs, 2, bscoefs,`*`)
colnames(bsScaled) = c('int', paste0('X', 1:10))

d3 = data.frame(x, y, bsScaled) %>% 
    gather(key=bs, value=bsfunc, -x, -y, factor_key = T) %>% 
    dplyr::filter(bsfunc >= min(y) & bsfunc <= max(y))

cs = RColorBrewer::brewer.pal(nlevels(d$x), 'Set3')
# as in fahrmeier
d3 %>%
  group_by(bs) %>%
  plot_ly() %>%
  add_markers(~x, ~y, color=I('rgba(0,0,0,.02)'), colors=cs, showlegend=T) %>% #RColorBrewer::brewer.pal(N, "Set3")
  add_lines(~x, ~bsfunc, color=~bs, colors=cs, showlegend=F) %>%
  layout()

# plot sum of basis functions
plotData = rbind(data.frame(d, sum=fitted(lmMod)),
                 data.frame(x=0, y=bscoefs[1], sum=NA))
filter(plotData) %>% 
    droplevels %>% 
    plot_ly() %>% 
    add_markers(~x, ~y, colors=cs, alpha=.5, showlegend=F) %>% 
    add_lines(~x, ~sum, colors=cs, showlegend=F) %>% 
    add_markers(x=0, y=bscoefs[1], size=~I(20), color=I('salmon'), alpha=.5, showlegend=F) %>% 
    layout()

# change degree of basis to 3
l = 3
bs = sapply(1:length(knots), function(k) ifelse(x >= knots[k], (x-knots[k])^l, 0))
lmModCubicSpline = lm(y ~ poly(x,3) + bs)


d %>%
    mutate(pred=fitted(lmModCubicSpline)) %>% 
    plot_ly() %>% 
    add_markers(~x, ~y, color=I(scales::alpha('black', .5)), colors=cs, showlegend=F) %>% 
    add_lines(~x, ~pred, colors=cs, showlegend=F) %>%
    config(displayModeBar=F) %>%
    layout()

# add other bases & mgcv comparison later