# Create project ----------------------------------------------------------

# Install packages --------------------------------------------------------
install.packages("ggplot2")
install.packages("palmerpenguins")
library(ggplot2)
library(palmerpenguins)

# Load data ---------------------------------------------------------------
data()
data(package = "palmerpenguins")
data(penguins)
penguins
str(penguins)

# Base R plotting ---------------------------------------------------------
plot(penguins)
plot(penguins$bill_length_mm, penguins$bill_depth_mm, col = penguins$species)

# Plotting with ggplot2 ---------------------------------------------------
ggplot(penguins, aes(x = bill_length_mm, y = bill_depth_mm, col = species)) +
  geom_point()

?geom_point

ggplot(penguins, aes(x = bill_length_mm, 
                     y = bill_depth_mm, 
                     col = species, 
                     size = body_mass_g)) +
  geom_point()

ggplot(penguins, aes(x = bill_length_mm, 
                     y = bill_depth_mm, 
                     col = species, 
                     shape = sex)) +
  geom_point()

p1 <- ggplot(penguins, aes(x = bill_length_mm, y = bill_depth_mm, col = species))

p1 + geom_point() 
p1 + geom_point() + geom_density_2d()
p1 + geom_point() + stat_ellipse()


# Line plots --------------------------------------------------------------
data(Orange)
Orange

p2 <- ggplot(Orange, aes(age, circumference, group = Tree))
p2 + geom_line() #+ geom_point()
p2 + geom_smooth(se = FALSE)

ggplot(penguins, aes(bill_length_mm, bill_depth_mm, col = species)) +
  geom_point() +
  geom_smooth(method = "lm")

# Bar, boxplot, and violin plots ------------------------------------------
ggplot(penguins, aes(species)) + 
    geom_bar(fill = "lightblue", col = "black", size = 1)

p3 <- ggplot(penguins, aes(species, bill_length_mm))

p3 + geom_bar(stat = "summary", 
              fun = "mean", 
              fill = "lightblue", 
              col = "black", 
              size = 1)

p3 + geom_bar(stat = "summary", 
              fun = "mean", 
              fill = "lightblue", 
              col = "black", 
              size = 1) +
    geom_point(position = position_jitter(0.2)) # with and without jitter

p3 + geom_boxplot(fill = "lightblue", outlier.color = "red") + 
    geom_point(alpha = 0.1)

p3 + geom_violin(fill = "lightblue", draw_quantiles = c(0.25, 0.75)) +
    geom_point(stat = "summary", fun = "median")

# Factorial bar graph -----------------------------------------------------
ggplot(penguins, aes(sex, bill_length_mm, fill = species)) +
    geom_bar(stat = "summary", fun = "median", col = "black", position = "dodge") +
    geom_point(position = position_dodge(0.9)) # jitter_dodge()

# Facetting ---------------------------------------------------------------
p4 <- ggplot(penguins, aes(x = bill_length_mm, y = bill_depth_mm, col = species)) +
    geom_point()

p4 + facet_wrap(~ species)
p4 + facet_wrap(~ species, nrow = 3)
p4 + facet_wrap(~ sex)
p4 + facet_grid(sex ~ year)

# Customizing labels ------------------------------------------------------
p4 + labs(title = "Penguin bill length and depth by species",
          subtitle = "Data collected on Torgersen, Biscoe, and Dream Islands",
          x = "Bill length (mm)", 
          y = "Bill depth (mm)", 
          col = "Species") -> p5

# Customizing the look of your plots with themes --------------------------
p5 + theme_bw()
p5 + theme_classic()
p5 + theme_test()
p5 + theme_dark()
p5 + theme_minimal()
p5 + theme_void()

# Customizing a theme -----------------------------------------------------
?theme

p5 +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_rect(fill = "beige"),
          legend.position = c(0.13, 0.12),
          legend.background = element_rect(colour = "black"))

# Saving plots ------------------------------------------------------------
dir.create("plots")
ggsave("plots/Graph 1.pdf", width = 6, height = 6, units = "in")
