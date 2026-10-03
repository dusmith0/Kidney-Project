#### Kidney Project
#


### Libraries
library('dplyr')
library("readr")
library("ggplot2")
library("patchwork")
library("rlang")
library("car")
library("glmnet")

##### Defining visual colors:
colors = c(
  "#E69F00",  #1 Orange #Trendline 1
  "#56B4E9",  #2 Sky Blue #Fill color
  "#009E73",  #3 Bluish Green #Trendline 3
  "#F5C710",  #4 Darker Amber #Trendline 2
  "#0072B2",  #5 Blue  -- Base color
  "#D55E00",  #6 Vermilion
  "#CC79A7",  #7 Reddish Purple # Background colors?
  "#000080",  #8 Deep Navy
  "#000000"   #9 Black
) #  okabe_ito  -- Color Blind scale



### Loading the data
kidney <- read.csv("ckd-dataset-v2.csv")[-c(1,2),]
unique(kidney['grf'])
  # I have no idea what 'p' could stand for...
kidney <- kidney[which(kidney['grf'] != ' p '),]
str(kidney)
head(kidney)

## Helper loop to allow for quicker ordinal values
## Should return the lower bound value of each interval.
ranking <- function(col){
  unique_values <- unique(col)
  ## subtracting a value for the strictly "<" is needed to keep the smallest value different from the bin above it
  numbers <- parse_number(unique_values) - ifelse(grepl("<", unique_values), 10, 0)
  orders <- unique_values[order(numbers)]
  
  # setting levels
  as.integer(factor(col, levels = orders)) - 1
}

## setting values to make ordinal
to_ordinal = c("bgr", "bu", "sod", "sc", "pot", "pcv", "rbcc", "wbcc", "age", "hemo", "sg", "grf", "al", "su")

ordinal <- kidney %>%
  mutate(across(all_of(to_ordinal), ranking)) %>%
  select(-c(affected))

names(ordinal)
View(ordinal)


## Creating Histograms for each predictor value
dist_check <- function(df){
  med_term <- names(df)[sapply(df, function(col) length(unique(col)) > 2)]

  plots <- lapply(med_term, function(col) { 
    ggplot(data = df, aes(x = !!sym(col))) + 
      geom_bar(fill = colors[5]) +
      labs(title = col, y = "frequency") + 
      theme_minimal() 
  })
  
  grid <- wrap_plots(plots) + 
    plot_annotation(
      title = "Histograms for distribution checking",
    )
  return(grid)
}

dist_check(ordinal)


## Running regression with Kidney Disease as target
ckd <- ordinal %>%
      select(-c(stage)) %>%
      select(-c(bp.limit, pcv, rbcc, sc, su)) %>% 
       ## dropping values Ali dropped due to collinearity
      mutate(class = ifelse(class == "ckd", 1, 0)) %>%
      mutate(across(where(is.character), as.factor))

fit <- glm(class ~ ., 
           family = binomial(link = "logit"),
           data = ckd)

summary(fit)
vif(fit)
## This fails to converge

## Trying this with GLMnet's regurlaization
target <- ckd$class
X <- ckd %>% select(-class)
X <- as.matrix(X)

cv <- cv.glmnet(X, target, family = "binomial", 
              alpha = 1)

coef(cv)













