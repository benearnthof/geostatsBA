test <- gam(site ~ s(distance_water), family = binomial(), data = evidence)
plot(test)
library(mgcViz)
p1 <- plot(getViz(test))
p2 <- p1 + 
  ylab("s(Distance to Water)") +
  xlab("Distance to water") +
  theme(text = element_text(size=15))
p2

test <- gam(site ~ s(lon, lat, bs = "gp", m = 2),
            family = binomial(),
            data = evidence)

p3 <- plot(getViz(test))

p4 <- p3 +
  ylab("Latitude") +
  xlab("Longitude") +
  theme(text = element_text(size = 15)) +
  ggtitle("Estimated Spatial Effect", subtitle = "Including no other Predictors")
p4
