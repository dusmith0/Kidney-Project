#### Kidney Project
#


##### Libraries
library('dplyr')
library("readr")
library("ggplot2") ## Graphics
library("patchwork") ## Graphics control
library("rlang") ## Macro type of package
library("car") ## General Statistical Package
library("glmnet") ## Regularized GLM package
library("pheatmap")
library("regclass") ## VIF on a Dataframe without specified target
library("ordinal") ## to run frequentest ordinal regression
library("brms") ## Bayesian models

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


## Removes unwanted variables, removing too much here will effect the code below. 
ordinal <- kidney %>%
  mutate(across(all_of(to_ordinal), ranking)) %>%
  mutate(across(where(is.character) & !any_of(c("class","stage")), as.numeric)) %>%
  select(-c(affected))



##### EDA --------------------------------------------------------  #####

## Creating Histograms for each predictor value
dist_check <- function(df){
  ## Removed Binary Predictors
  med_term <- names(df)[sapply(df, function(col) length(unique(col)) > 2)]

  plots <- lapply(med_term, function(col) { 
    ggplot(data = df, aes(x = !!sym(col))) + 
      geom_bar(fill = colors[5]) +
      labs(y = "frequency") + 
      theme_minimal() 
  })
  
  grid <- wrap_plots(plots) + 
    plot_annotation(
      title = "Histograms for distribution checking",
    )
  return(grid)
}

dist_check(ordinal)


## Creating Multi box plot By stage and predictor
## As R does not handle binary the same as Python, I removed the binary predictors. 
box_plot <- function(df, target = "stage") {
  med_term <- setdiff(names(df)[sapply(df, function(col) length(unique(col)) > 2)],
                      target)
  
  plots <- lapply(med_term, function(col) {
    ggplot(data = df, aes(x = factor(!!sym(target)), y = !!sym(col))) + 
      geom_boxplot(fill = colors[5], alpha = .90, outlier.size = 0.8) + 
      theme_minimal(base_size = 9) + 
      labs(y = col, x = NULL) + 
      theme(
        axis.title.y = element_text(size = 14)
      )
  })
  
  grid <- wrap_plots(plots, ncol = 4) + 
    plot_annotation(
      title = "Box-Plot of Predictor by Stage",
      theme = theme(plot.title = element_text(size = 20, face = "bold"))
      )
  return(grid)
}  

grid_plot <- box_plot(ordinal)
ggsave("Descrete_boxplots.png", grid_plot, width = 16, height = 15, dpi = 300)


## Creating Heat Map of before VIF and after VIF corrections
ordinal_corr <- ordinal %>%
  select(-c(class, stage))

heat_map <- pheatmap(cor(ordinal_corr), display_numbers = TRUE)
ggsave("Heatmap.png", heat_map, width = 16, height = 15, dpi = 300)

## Running VIF calculation  -- R usually needs a mock model to run VIF
## the "performance" package maybe able to handle it without the model.
vif_check <- ordinal %>%
  select(where(is.numeric))

# This seems to work, and can recreate the vif, unfortunately it will remove the 
# random variable chosen as the dependent variable. Switching will create the same VIF
# values, and will return the missing one. A loop spit out all, but that seem uncessisary. 
fit <- lm(vif_check[[1]] ~ ., data = vif_check[-1])
vif(fit)



vif_check_new <- ordinal %>%
  select(-c(stage, grf, bp.limit, hemo, class))%>%
  select(where(is.numeric))

fit <- lm(vif_check_new[[2]] ~ ., data = vif_check_new[-2])
vif(fit)


##### Modeling ----------------------------------------------------  #####

## Removing the Unwanted Predictors: 
ordinal_red <- kidney %>%
  mutate(across(all_of(to_ordinal), ranking)) %>%
  mutate(across(where(is.character) & !any_of(c("class","stage")), as.numeric)) %>%
  select(-c(affected, grf, bp.limit, hemo, class))


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


### Ordinal Regression
ordinal_red$stage <- factor(ordinal_red$stage, ordered = TRUE)

model <- clm(stage~., data = ordinal_red,
             link = "logit",
             control = clm.control(
               method = "optim", # This should use the BFGS algorithm for convergence
               maxIter = 5000,
               trace = 0
             ))
summary(model)


### Bayesian Ordinal Logistic Regression
# Define the Horseshoe prior for regression coefficients (class = "b")
# df = 1 specifies the Student-t degrees of freedom for local shrinkage
hs_prior <- set_prior(horseshoe(df = 1), class = "b")

# Fit the ordinal logistic model (cumulative logit family)
fit_horseshoe <- brm(
  formula = stage ~ .,
  data = ordinal_red,
  family = cumulative(link = "logit"),
  prior = hs_prior,
  chains = 4,
  cores = 4,             # Parallel processing across CPU cores
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.99) # Prevents divergent transitions common with Horseshoe priors
)

# Model Summary and Visualization
summary(fit_horseshoe)
plot(fit_horseshoe)

unique(ordinal_red$sod)
### Running Cumulative Ordinal Models
df_bayes <- ordinal_red %>%
  select(bp..Diastolic., bu, sod, pcv, age, stage) %>%
  mutate(
    
    stage = ordered(stage), 

    sod_factor = factor(case_when(
      sod <= 3 ~ "1",
      sod == 4 ~ "2",
      sod == 5 ~ "3",
      sod >= 6 ~ "4"),
      levels = c("2", "1", "3", "4") # To set "2" as the baseline
      ),
    ## Selecting BU measures
    ## Daka BU is 2.14 times the western equivalence, and common literature
    ## uses Serum Urea instead of blood urea in mmol/L instead of mg/L
    ## the conversion is Serum Urea * (2.8 * 2.14) = Blood Urea
    ## Literature suggests T1 <10.5 SU == 62.916
    ##                     T2 < 15.1 SU == 90.4792
    ##                     T3 > 15.1 SU == 90.4792
    bu_factor = factor(case_when(
      bu == 0 ~ "0_baseline",
      bu == 1 ~ "1_mild",
      bu == 2 ~ "2_mod",
      bu == 3 ~ "3_mod",
      bu == 4 ~ "4_sev",
      bu == 5 ~ "5_sev",
      bu == 6 ~ "6_sev",
      bu == 7 ~ "7_extreme"
    ))
  ) %>%
  select(-sod, - bu)


priors <- c(
  # Flat Prior defined by Alistar
  set_prior("normal(0, 50)", class = "b"),
  
  # Custom Priors based on Literature
  # Sod
  set_prior("normal(0.37, 0.2)", class = "b", coef = "sod_factor1"),
  set_prior("normal(0.16, 0.2)", class = "b", coef = "sod_factor3"),
  set_prior("normal(0.20, 0.2)", class = "b", coef = "sod_factor4"),
  
  # Bu
  set_prior("normal(0.25, 0.20)", class = "b", coef = "bu_factor1_mild"),
  set_prior("normal(0.55, 0.20)", class = "b", coef = "bu_factor2_mod"),
  set_prior("normal(0.55, 0.20)", class = "b", coef = "bu_factor3_mod"),
  set_prior("normal(0.55, 0.20)", class = "b", coef = "bu_factor4_sev"),
  set_prior("normal(0.55, 0.20)", class = "b", coef = "bu_factor5_sev"),
  set_prior("normal(0.55, 0.20)", class = "b", coef = "bu_factor6_sev"),
  set_prior("normal(0.55, 0.20)", class = "b", coef = "bu_factor7_extreme"),
  
  # pcv and age
  set_prior("normal(8, 2)", class = "b", coef = "pcv"),
  set_prior("normal(0, 0.5)", class = "b", coef = "age"),
  set_prior("normal(0, 0.5)", class = "b", coef = "bp..Diastolic.")
)

# fitting the model
cumulative_model <- brm(
  formula = stage ~ .,
  data = df_bayes,
  family = cumulative(link = "logit"),
  prior = priors,
  chains = 4, 
  iter = 2000, 
  warmup = 1000,
  cores = 4,
  seed = 123
)


summary(cumulative_model)







