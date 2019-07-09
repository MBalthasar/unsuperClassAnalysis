#' Unsupervised Classification Saturation Function
#'
#' This function aims to identify the number of classes required for an unsupervised
#' classification in order to get the highest accuracy results when mapping the
#' distribution of a single class.
#'
#' @param raster A raster file.
#'
#' @param nSamples Integer. Number of random samples to draw to fit cluster map.
#'
#' @param nClasses Vector of integers, containing the number of classes.
#'
#' @param nStarts Integer. Number of random starts for kmeans algorithm.
#'
#' @param valData SpatialPolygonsDataFrame or SpatialPointsDataFrame with validation data.
#'
#' @param responseCol Character. Column containing the validation data in attribute table of valData.
#'
#' @param valSamples Integer. Number of pixels to sample for validation (only applies to polygons).
#'
#' @return A list containing:
#' \itemize{
#'   \item Raster. The output image of the best unsupervised classification containing the presence
#'   probability of the desired class' presence.
#'   \item Data.Frame. Overall accuracies for each class number.
#'   \item Classificaation metrics for the best unsupervised classification.
#' }
#'
#' @details The output images contains values ranging from 1 - 4. These pixel values represent:
#' \itemize{
#'   \item 1 - High probability
#'   \item 2 - Medium probability
#'   \item 3 - Low probability
#'   \item 4 - Least probability
#' }
#'
#' @examples
#' library(sp)
#' library(raster)
#' library(RStoolbox)
#' library(ggplot2)
#'
#' # Load sample raster file
#' my_raster <- raster::brick(system.file(package = "unsuperClassAnalysis",
#'                            "extdata", "landsat_sample.tif"))
#'
#' # Load sample validation data
#' my_val <- raster::shapefile(system.file(package = "unsuperClassAnalysis", "extdata",
#'                             "validation_sample.shp"))
#'
#' # Get vector of class name for each polygon and check if it matches with "deforest"
#' # If there's a match it will return 1, else it will return NA
#' deforest_vec <- match(my_val@data[,2], "deforest")
#'
#' # Replace all NAs with 0
#' deforest_vec[is.na(deforest_vec)] <- 0
#'
#' # Add new binary vector as new column in shapefile
#' # Now every polygon has a binary code whether it's deforested area (1) or not (0)
#' my_val[["Unsup_Id"]] <- deforest_vec
#'
#' # Plot raster data with validation data on top
#' raster::plotRGB(my_raster, 3,2,1, stretch="lin")
#' plot(my_val[my_val$Unsup_Id == 0,], col="red", add=T)
#' plot(my_val[my_val$Unsup_Id == 1,], col="blue", add=T)
#'
#' # Execute function
#' x <- unsupSaturation(raster = my_raster,
#'                     nSamples = 1000,
#'                     nClasses = c(2,3,5,10),
#'                     nStarts = 25,
#'                     valData = my_val,
#'                     responseCol = "Unsup_Id",
#'                     valSamples = 1000)
#'
#' # Get best unsupervised classification
#' x1 <- x[[1]]
#'
#' # Get data.frame with accuracy statistics
#' x2 <- x[[2]]
#'
#' # Get accuracy matrix of best classification
#' x3 <- x[[3]]
#'
#' # Get unique values of classification
#' unique_vals <- unique(x1)
#'
#' # Change numeric (continuous) values to factor (discrete) for discrete legend
#' x1[] <- as.factor(x1@data@values)
#'
#' # Plot with ggplot
#' RStoolbox::ggR(x1, geom_raster = T)+
#'  ggtitle("Unsupervised Classification of Deforested Areas\n")+
#'  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14))+
#'  theme(legend.title = element_text(size = 11, face = "bold"))+
#'  theme(legend.text = element_text(size = 10))+
#'  scale_fill_manual(values = c("red", "khaki", "beige"),
#'                    labels = c("High", "Medium", "Least"),
#'                    name = "Deforestation\nProbability")+
#'  xlab("")+
#'  ylab("")
#'
#' @export
unsupSaturation <- function(raster, nSamples, nClasses, nStarts, valData, responseCol, valSamples){
  # Define presence data ####
  presence_data <- valData
  presence_values <- valData@data[[responseCol]]
  presence_data[["new_column"]] <- presence_values
  presence_data <- presence_data[presence_data$new_column == 1,]
  presence_data <- raster::aggregate(presence_data)

  # Define max class number and new probaility pixel vals
  max_class_number <- max(nClasses)
  high_prob_value <- max_class_number +1
  medium_prob_value <- max_class_number +2
  low_prob_value <- max_class_number +3
  least_prob_value <- max_class_number +4

  # Create df which will be filled with accuracy information
  my_df <- data.frame(matrix(data = NA, nrow = length(nClasses), ncol = 2))
  names(my_df) <- c("classes", "accuracy")
  my_df[,1] <- nClasses

  # start for loop for each class number
  for (i in 1:length(nClasses)){
    print(paste0("Starting with classification ", i, " of ", length(nClasses)))
    # Define current class numbers
    current_class_number <- nClasses[i]
    # execute unsupervised classification
    set.seed(7)
    current_unsup <- RStoolbox::unsuperClass(raster, nSamples = nSamples,
                                             nClasses = current_class_number,
                                             nStarts = nStarts)
    # Get current map from unsupervised classification
    current_map <- current_unsup$map
    # Extract values from current classification based on presence data and ignore NA values
    current_values <- raster::extract(current_map, presence_data)
    current_values <- current_values[[1]][!is.na(current_values[[1]])]
    # List frequency of all values and reorder them with in decreasing manner
    current_values_table <- table(current_values)
    current_values_table_sort <- sort(current_values_table, decreasing = T)

    # Check how many different pixel values are available
    if (length(current_values_table_sort) >= 3){
      current_map[current_map == as.numeric(names(current_values_table_sort[1]))] <- high_prob_value
      current_map[current_map == as.numeric(names(current_values_table_sort[2]))] <- medium_prob_value
      current_map[current_map == as.numeric(names(current_values_table_sort[3]))] <- low_prob_value
      current_map[current_map < high_prob_value] <- least_prob_value
    } else if (length(current_values_table_sort) == 2){
      current_map[current_map == as.numeric(names(current_values_table_sort[1]))] <- high_prob_value
      current_map[current_map == as.numeric(names(current_values_table_sort[2]))] <- medium_prob_value
      current_map[current_map < high_prob_value] <- least_prob_value
    } else {
      current_map[current_map == as.numeric(names(current_values_table_sort[1]))] <- high_prob_value
      current_map[current_map < high_prob_value] <- least_prob_value
    }

    # create map for validation
    val_map <- current_map
    val_map[val_map == high_prob_value] <- 1
    val_map[val_map != 1] <- 0
    # start validation
    current_validatation <- RStoolbox::validateMap(val_map, valData = valData,
                                                   responseCol = responseCol,
                                                   nSamples = valSamples,
                                                   mode = "classification")
    current_accuracy <- current_validatation$performance$overall[1]
    # fill df with accuracy information
    my_df[i,2] <- current_accuracy

    # check highest accuracy
    if (i == 1){
      final_accuracy <- current_accuracy
      final_validation <- current_validatation
      final_map <- current_map
      final_class_number <- current_class_number
    } else {
      # Check if current accuracy is higher than the one before
      if (current_accuracy > final_accuracy){
        final_accuracy <- current_accuracy
        final_validation <- current_validatation
        final_map <- current_map
        final_class_number <- current_class_number
      }
    }
  }

  # Change pixel probability values
  # get number of classes from final map
  final_classes <- unique(final_map@data@values)
  # remove NA values
  final_classes <- final_classes[!is.na(final_classes)]

  if (length(final_classes) == 4){
    final_map[final_map == high_prob_value] <- 1
    final_map[final_map == medium_prob_value] <- 2
    final_map[final_map == low_prob_value] <- 3
    final_map[final_map > 3] <- 4
  } else if (length(final_classes) == 3){
    final_map[final_map == high_prob_value] <- 1
    final_map[final_map == medium_prob_value] <- 2
    final_map[final_map > 2] <- 4
  } else {
    final_map[final_map == high_prob_value] <- 1
    final_map[final_map > 1] <- 4
  }

  # Plot my_df
  print(ggplot2::ggplot(data=my_df, aes(x=classes, y=accuracy, group=1)) +
          geom_line()+
          geom_point()+
          labs(title="Overall accuracy per class number",
               x="class number", y = "Overall accuracy")+
          theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
                axis.text=element_text(size=10),
                axis.title=element_text(size=12),
                legend.title=element_text(size=13, face="bold"),
                legend.text=element_text(size=12)))

  # Change names of data.frame
  names(my_df) <- c("No_of_Classes", "Overall_Accuracy")

  # Print accuracy statistics and ideal number of classes
  print(final_validation)
  print(paste0("Ideal number of classes is: ", final_class_number))

  # Create output list
  my_output <- list(final_map, my_df, final_validation)
  return(my_output)
}
