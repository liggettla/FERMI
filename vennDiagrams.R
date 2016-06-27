library(VennDiagram)

# single diagram
grid.newpage()
draw.single.venn(area = 10, category = "Something")

## lty - outline of cirlces
## fill - colour
## alpha - colour transparency
grid.newpage()
draw.single.venn(22, category = "Dog People", lty = "blank", fill = "cornflower blue", 
                 alpha = 0.5)

# Diagram with two circles
grid.newpage()
draw.pairwise.venn(area1 = 22, area2 = 20, cross.area = 11, category = c("Dog People", 
                                                                         "Cat People"))
## cat.pos - position of category titles, represented by degree from the
## middle of the circle
## cat.dist - distance of the category titles from the edge of the circle
grid.newpage()
draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), lty = rep("blank", 
  2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
  0), cat.dist = rep(0.025, 2))

## scaled - TRUE for scaled or FALSE for unscaled cirlces
grid.newpage()
draw.pairwise.venn(22, 20, 11, category = c("Dog People", "Cat People"), lty = rep("blank", 
  2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
  0), cat.dist = rep(0.025, 2), scaled = FALSE)

# non-overlapping circles
## euler.d - TRUE for movable circles; FALSE for unmovable circles. Must be
## TRUE to have space between non-overlapping circles.
## sep.dist - distance between circles
## rotation.degree - degrees the diagram is rotated
grid.newpage()
draw.pairwise.venn(area1 = 22, area2 = 6, cross.area = 0, category = c("Dog People", 
  "Snake People"), lty = rep("blank", 2), fill = c("light blue", "green"), 
  alpha = rep(0.5, 2), cat.pos = c(0, 180), euler.d = TRUE, sep.dist = 0.03, 
  rotation.degree = 45)

# diagram with 3 circles
grid.newpage()
draw.triple.venn(area1 = 22, area2 = 20, area3 = 13, n12 = 11, n23 = 4, n13 = 5, 
                 n123 = 1, category = c("Dog People", "Cat People", "Lizard People"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

# using subset function to pull values for the diagram where d is imported data file
grid.newpage()
draw.triple.venn(area1 = nrow(subset(d, Dog == 1)), area2 = nrow(subset(d, Cat == 
  1)), area3 = nrow(subset(d, Lizard == 1)), n12 = nrow(subset(d, Dog == 1 & 
  Cat == 1)), n23 = nrow(subset(d, Cat == 1 & Lizard == 1)), n13 = nrow(subset(d, 
  Dog == 1 & Lizard == 1)), n123 = nrow(subset(d, Dog == 1 & Cat == 1 & Lizard == 
  1)), category = c("Dog People", "Cat People", "Lizard People"), lty = "blank", 
  fill = c("skyblue", "pink1", "mediumorchid"))