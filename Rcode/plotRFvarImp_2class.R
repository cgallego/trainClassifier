plotRFvarImp_2class <- function(allselvarImp) {
  library("lattice")
  library("latticeExtra")
  require(ggplot2)
  
  # plot
  print(
  ggplot(allselvarImp, aes(x=reorder(allselvarImp$features, allselvarImp$varImp), y=varImp, fill=SelectedFeatureGroup )) +
  geom_bar() + coord_flip() + 
    theme_bw(base_size = 14) + 
    labs(title = "Comparison of selected features for Mass and Non-mass groups",
         x="", y="variable importance ranking")
  )

  return()  
}