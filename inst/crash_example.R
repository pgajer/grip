rds.file <- "~/current_projects/ZB/analysis_output/mg/mg_spp_k5_graph.rds"
g <- readRDS(rds.file)

names(g)
## [1] "adj_list"    "weight_list"

x <- grip.layout(
+     adj_list = g$adj_list,
+     weight_list = g$weight_list,
+     dim = 3,
+     final_rounds = 20,
+     num_init = 16,
+     num_nbrs = 10,
+     r = 0.15,
+     s = 3.0,
+     tinit_factor = 6,
+     seed = 6
+ )

head(x)
##      [,1] [,2] [,3]
## [1,]  NaN  NaN  NaN
## [2,]  NaN  NaN  NaN
## [3,]  NaN  NaN  NaN
## [4,]  NaN  NaN  NaN
## [5,]  NaN  NaN  NaN
## [6,]  NaN  NaN  NaN
>
