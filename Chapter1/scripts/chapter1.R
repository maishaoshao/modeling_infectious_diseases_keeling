#Note------
#script for this paper: https://www.science.org/doi/full/10.1126/science.1061020?versioned=true 
#The Foot-and-Mouth Epidemic in Great Britain: Pattern of Spread and Impact of Interventions
#date: 2025-03-01

#load the library--------
library(readxl)
library(sparr)
library(dplyr)
library(lubridate)
library(geosphere)
library(ggplot2)
library(RColorBrewer)  # For color palettes
library(haven)
library(nlme)
library(mgcv)
library(hrbrthemes)
library(sparr)
library(purrr)  # For the map function
library(knitr)
library(rcompanion)



#set wd-------
setwd("/Users/...")
#load the data----
#uk data
d_2001 <- read_excel("UK_FAM_2001.xlsx")
d_1967 <- read.csv2("UK_FAM_1967-68.csv", 
                    sep = ",")
#stimulate data
d2 <- read.csv2("epi.tab.run1.csv",
                sep = ",")

#uk Cumbria data
data(fmd) #from sparr package

#1.UK data 2001 and 1967-----
#clean the data
#add cases since the start of the pandemic
#add days since the start of the pandemic
d_2001 %>%
  dplyr::select(Date, `New cases`) %>%
  rename(date = Date,
         new_cases = `New cases`) %>% 
  mutate(date = as.Date(date)) %>%
  as.data.frame() %>%
  arrange(date) %>%
  mutate(cum_cases = cumsum(new_cases),
         days = as.numeric(date - min(date))) -> d_2001


head(d_2001, 10)
summary(d_2001)
str(d_2001) #'data.frame':	41 obs. of  4 variables:

#add cases since the start of the pandemic
#add days since the start of the pandemic

d_1967 %>%
  rename(date = day.of.the.month,
         new_cases = daily_total_of_outbreaks) %>%
  mutate(date = as.Date(dmy(date))) %>%
  arrange(date) %>%
  mutate(cum_cases = cumsum(new_cases),
         days = as.numeric(date - min(date))) -> d_1967

str(d_1967) #'data.frame':	224 obs. of  4 variables

d_1967$year <- "1967"
d_2001$year <- "2001"
rbind(d_1967, d_2001) -> d


### plot the data -----
d %>%
  filter(days <= 60) %>%
  ggplot(aes(x = days, y = new_cases, colour = year, fill = year)) +
  geom_bar(stat = "identity", position = "dodge") +  # Changed from geom_histogram()
  theme_bw() +
  theme(legend.position = "right") +
  labs(title = "Foot-and-Mouth Epidemic in Great Britain",
       x = "Date",
       y = "New cases")

#save the plot
ggsave("foot_and_mouth_epidemic.png", width = 10, height = 6, dpi = 300)


#2.stimualate data----
head(d2)
str(d2)
d2 %>%
  mutate(latitude = as.numeric(latitude),
         longitude = as.numeric(longitude)) -> d2

distances_to_first_km <- round(sapply(2:nrow(d2), function(i) {
  distHaversine(d2[1, c("longitude", "latitude")], d2[i, c("longitude", "latitude")])
}) / 1000,2)

distances_to_first_km <- c(0, distances_to_first_km)
d2$distances_to_first_km <- distances_to_first_km
range(d2$distances_to_first_km) #0.00 185.76

d2 %>%
  mutate(distances_to_first_km_cat = cut(distances_to_first_km, 
                                         breaks = seq(0, 200, by = 20))) -> d2  

###plot the data----
head(d2)
# Ensure the distance category is a factor and set the levels in order
d2$distances_to_first_km_cat <- factor(d2$distances_to_first_km_cat, levels = unique(d2$distances_to_first_km_cat))

# Creating a bar plot of distance categories
ggplot(d2, aes(x = distances_to_first_km_cat)) +
  geom_bar(fill = "steelblue") +  # You can adjust colors here
  labs(title = "Frequency of Cases by Distance Category to First Case",
       x = "Distance to First Case (km)",
       y = "Count of Cases") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability

#save the plot
ggsave("distance_to_first_case.png", width = 10, height = 6, dpi = 300)

colors <- brewer.pal(n = length(unique(d2$herd.type)), name = "Set2")
ggplot(d2, aes(x = distances_to_first_km_cat, fill = herd.type)) +
  geom_bar(position = "stack") +  # Stack is default but included for clarity
  scale_fill_manual(values = colors, name = "Herd Type") +  # Manually set colors with legend title
  labs(title = "Frequency of Cases by Distance Category to First Case",
       x = "Distance to First Case (km)",
       y = "Count of Cases") +
  theme_bw() +  # Minimal theme with some tweaks
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Improve x-axis labels readability
        plot.title = element_text(face = "bold", hjust = 0.5),  # Center and bold the title
        legend.position = "right",  # Adjust legend position
        legend.title.align = 0.5) +  # Center align legend title
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) + 
  facet_wrap(~ herd.type) # Adjust legend guide

#save the plot
ggsave("distance_to_first_case_herd_type.png", width = 10, height = 6, dpi = 300)

#3. UK Cumbria data----
f <- spattemp.density(fmd$cases,h=6,lambda=8)
g <- bivariate.density(fmd$controls,h0=6)
fmdrr <- spattemp.risk(f,g,tolerate=TRUE)
plot(fmdrr,sleep=0.1,fix.range=TRUE)
plot(fmdrr,type="conditional",sleep=0.1,tol.type="two.sided",
     tol.args=list(levels=0.05,drawlabels=FALSE)) 

#check if the pattern is clustered, random, or regular
k_test <- Kest(fmd$cases)  # Ripley's K function
plot(k_test)

#envelope test for significance
envelope_test <- envelope(fmd$cases, Kest, nsim = 10)
plot(envelope_test, main = "Pointwise Envelope Test for fmd$cases")

#extract data from fmd
if(is.ppp(fmd$cases)) {
  # Extract coordinates
  coords <- data.frame(x = fmd$cases$x, y = fmd$cases$y)
  
  # Check if there are marks and include them
  if(!is.null(marks(fmd$cases))) {
    coords$marks <- marks(fmd$cases)
  }
  
  # Now, coords is a data frame with the x, y coordinates and possibly marks
  print(coords)
}

summary(coords)  # Get a summary of your data
head(coords)
ggplot(coords, aes(x = x, y = y, color = marks)) +
  geom_point(alpha = 0.6) +  # Adjust transparency with alpha
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient from low to high marks
  labs(title = "Spatial Distribution of Cases", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()

#
library(deSolve)

# Model definition
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS1 <- -beta * S1 * I1 + m21 * S2 - m12 * S1
    dI1 <- beta * S1 * I1 - gamma * I1 + m21 * I2 - m12 * I1
    dR1 <- gamma * I1 + m21 * R2 - m12 * R1
    dS2 <- -beta * S2 * I2 + m12 * S1 - m21 * S2
    dI2 <- beta * S2 * I2 - gamma * I2 + m12 * I1 - m21 * I2
    dR2 <- gamma * I2 + m12 * R1 - m21 * R2
    
    return(list(c(dS1, dI1, dR1, dS2, dI2, dR2)))
  })
}

# Initial state values
initial_state <- c(S1 = 990, I1 = 10, R1 = 0, S2 = 1000, I2 = 0, R2 = 0)

# Parameters
parameters <- c(beta = 0.4, gamma = 0.1, m12 = 0.05, m21 = 0.05)

# Time
times <- seq(0, 50, by = 0.1)  # More frequent time points for stability

# Solve ODE
output <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

# Convert output to a dataframe
output_df <- as.data.frame(output)

# Plot results
library(ggplot2)
output_long <- reshape2::melt(output_df, id.vars = "time")
ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(title = "SIR Model Across Two Patches", x = "Time (days)", y = "Population count",
       color = "Compartment") +
  theme_minimal()

