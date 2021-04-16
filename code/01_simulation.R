library(ggplot2)
library(dplyr)

# Setup key variables
n_cells <- 1000000 # In theory, we'd want the number of T cells in a mouse but any sufficiently large number is fine
n_SLO <- 400
n_chemotaxis <- 5 # number of SLOs affected by chemotaxis
n_count_SLOs_as_hits <- 1 # how many SLOs do we want to see as "seen" ? 5? 1?
weight_chemotaxis <- 0.1 # number between 0 and 1 that ascribes the amount of probability to be placed on the SLOs
non_chemotaxis_SLOs <- n_SLO-n_chemotaxis
n_rounds <- 16 # One round = 12 hours; so n_rounds/2 is number of days

scenario_1_mat <- matrix(rep(0, n_rounds*n_cells), nrow = n_cells)
scenario_2_mat <- scenario_1_mat

set.seed(14651)

# Loop is nasty in R but most readable
for(i in 1:n_rounds){
  
  # assign a SLO for each cell under a uniform probability model
  slo_assign_scenario1 <- sample(1:n_SLO, size = n_cells, replace = TRUE, prob = rep(1/n_SLO, n_SLO) )
  
  # assign a SLO for each cell under a model accounting for chemotaxis variation
  
  prob_vec_scenario2 <- c(rep(weight_chemotaxis/n_chemotaxis, n_chemotaxis), 
                rep((1-weight_chemotaxis)/(non_chemotaxis_SLOs), non_chemotaxis_SLOs))
  slo_assign_scenario2 <- sample(1:n_SLO, size = n_cells, replace = TRUE, prob = prob_vec_scenario2 )
  
  if(i > 1){
    scenario_1_mat[,i] <- (slo_assign_scenario1 <= n_count_SLOs_as_hits) | scenario_1_mat[,i-1]  # for cumulative counting
    scenario_2_mat[,i] <- (slo_assign_scenario2 <= n_count_SLOs_as_hits) | scenario_2_mat[,i-1]  # for cumulative counting
  } else {
    scenario_1_mat[,i] <- slo_assign_scenario1 <= n_count_SLOs_as_hits
    scenario_2_mat[,i] <- slo_assign_scenario2 <= n_count_SLOs_as_hits
  }
}

data.frame(
  DPI = 1:n_rounds/2,  # divide by 2 b/c we're saying there's turnover every 
  scenario1 = colMeans(scenario_1_mat), 
  scenario2 = colMeans(scenario_2_mat)
) %>% reshape2::melt(id.vars = "DPI") -> melted_df

ggplot(melted_df, aes(x = DPI, y = value *100, color = variable)) + 
  geom_point() + geom_line() +
  labs(x = "DPI", y = "% of cells that have encountered SLOs of interest", color = "") +
  theme_classic()+ scale_color_manual(values = c("firebrick", "dodgerblue3"))
  