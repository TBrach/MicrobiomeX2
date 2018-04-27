# MicrobiomeX Pipeline (MXP)

The MXP illustrates a streamlined analysis of microbiome data starting with a phyloseq object. 
The analysis includes exploratory analysis, alpha-diversity analysis, beta-diversity (ordination) analysis, and comparison of phyla ratios.

- Where to start?: 
    - The pipeline makes use of "R Projects" and the "bookdown package". 
        - Start a new project in R Studio using this folder as "from existing directory". Then open the index.Rmd file in R studio.
        - When you later build the final document with bookdown, all Rmd files (the Chapters) in the folder are considered in alphabetic order. For more info: <https://www.rstudio.com/resources/webinars/introducing-bookdown/>
- What you need?:
    - make sure that all packages under **Load packages (and color schemes)** are installed and up to date.
    - you need a prepared phyloseq object, saved as **.rds** file on your computer (saveRDS command) (best in the project folder). Just like the one in Example_phyloseq_objects. Info: <https://joey711.github.io/phyloseq/import-data.html>
- What then?
    - Adjust the input parameters for the analysis under **Load phyloseq object and set input**. I hope most input is self-explanatory in particular with the annotation given in the section. In this section: all lines where user input is required or possible are marked with **NB: USER INPUT** and **NB: user input**, respectively
       - NB: In the section: **Load phyloseq object and set input**, I intentionally put the input tests directly after the input is given. That way it takes longer to write the input. I did so, because I think it is a good idea to go once through the entire section for each analysis. It helps in understanding the input and discovering potential problems.
- How to run it?
    - When you are certain about the input, you build the document in Rstudio with bookdown. You do so under Build >> Build All. Or by using the Build tab >> Build Book. This will run the entire analysis and produce the output document.
        - By default I use bookdown::html_document2, so an html file **_main.html (look at it in your browser)** is created.
        - You can change the output to pdf or gitbook by setting the output: on top of the document to bookdown::pdf_book or bookdown::gitbook, respectively. However, wide tables that are wider than the page size are currently only fully visible in html.


   