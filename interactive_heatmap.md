# Interactive ComplexHeatmap (R Code)
**Author:** Shadi Mohammadabadi  
**Purpose:** This markdown file contains R code for an RNA-seq analysis pipeline originally written in Jupyter Notebook.

### 1. Load Libraries
```r
# Load required libraries
library(shiny)
library(InteractiveComplexHeatmap)
```
### 2. Read a heatmap object
```r
# Read a heatmap object in .rds format
ht <- readRDS("heatmap_obj.rds")
```
### 3. Load the Interactive heatmap
```r
ui = fluidPage(
    InteractiveComplexHeatmapOutput(width1 = 600, height1 = 600)  
)

server = function(input, output, session) {
    makeInteractiveComplexHeatmap(input, output, session, ht)
}

shiny::shinyApp(ui, server)
```



