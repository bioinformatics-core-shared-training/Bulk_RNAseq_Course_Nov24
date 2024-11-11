# Introduction to Bulk RNA-seq data analysis
### 2nd, 3rd, 4th Oct 2024
#### In person
#### Bioinformatics Training Facility, Craik-Marshall Building, Downing Site, University of Cambridge

![](Bulk_RNAseq_Course_Base/images/CRUK_Cambridge_Major_Centre_logo.jpg)

## Instructors

* Ashley Sawle (CRUK Cambridge Institute)
* Abbi Edwards (CRUK Cambridge Institute)
* Betty Wang (Dpt Clinical Neurosciences, University of Cambridge)
* Ian Tsang (NIAB)
* Jiayin Hong (Dpt Biochemistry, University of Cambridge)
* Ramy Ragheb (Stem Cell Institute, University of Cambridge)

## Outline

In this workshop, you will be learning how to analyse RNA-seq data. This will
include read alignment, quality control, quantification against a reference,
reading the count data into R, performing differential expression analysis, and
gene set testing, with a focus on the DESeq2 analysis workflow. You will learn
how to generate common plots for analysis and visualisation of gene expression
data, such as boxplots and heatmaps.

This workshop is aimed at biologists interested in learning how to perform
differential expression analysis of RNA-seq data.

> ## Prerequisites
>
> __**Some basic experience of using a UNIX/LINUX command line is assumed**__
>
> __**Some R knowledge is assumed and essential. Without it, you will struggle on this course.**__
> If you are not familiar with the R statistical programming language we
> strongly encourage you to work through an introductory R course before
> attempting these materials.
> We recommend our [Introduction to R course](https://bioinformatics-core-shared-training.github.io/r-intro/)

## Timetable

### Day 1

**Trainers:** Abbi, Betty, Ian, Ramy, Hugo

9:30 - 9:45 - [Welcome](https://docs.google.com/presentation/d/13mwQfPCpYzj0iz5EyTQHXRd1KNwfds2ILdY8EiK2020/edit?usp=sharing)! - Bajuna

9:45 - 10:15 - [Introduction to RNAseq Methods](Bulk_RNAseq_Course_Base/Markdowns/01_Introduction_to_RNAseq_Methods.html) - Betty

10:15 - 11:00 [Raw read file format and QC](Bulk_RNAseq_Course_Base/Markdowns/02_FastQC_introduction.html) - Betty

  - [Practical](Bulk_RNAseq_Course_Base/Markdowns/02_FastQC_practical.html)
  - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/02_FastQC_solutions.html)

11:00 - 13:30 [Alignment and Quantification of Gene Expression with Salmon](Bulk_RNAseq_Course_Base/Markdowns/03_Quantification_with_Salmon_introduction.html) - Hugo

  - [Practical](Bulk_RNAseq_Course_Base/Markdowns/03_Quantification_with_Salmon_practical.html)  
  - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/03_Quantification_with_Salmon_solutions.html)

13:30 - 14:30 Lunch

14:30 - 15:30 [QC of alignment](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_introduction.html) - Abbi

  - [Practical](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_practical.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_practical.pdf))  
   - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_solutions.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_solutions.pdf))

15.30 - 17.30 [Data Exploration in R](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration_Slides.html) - Ian

  - [Practical](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.pdf))
  <!-- - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration_solutions.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration_solutions.pdf)) -->
  <!-- - [Live Script](live_scripts/20240621_initial_exploration.R) -->


### Day 2

**Trainers:** Ash, Abbi (AM), Ian, Ramy

9:30 - 10:15  [Introduction to RNAseq Analysis in R](Bulk_RNAseq_Course_Base/Markdowns/06_Introduction_to_RNAseq_Analysis_in_R.html) - Abbi

10:15 - 11:30 [Statistical Analysis of Bulk RNAseq Data](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/07_Statistics_1.pdf) - Abbi

11:30 - 17:30 [Linear Models in R and DESeq2 (Slides)](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/07_Statistics_2.pdf) - Abbi/Ash

  - [Practical - Differential Expression for RNA-seq](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2.pdf)) - Ash/Abbi
  <!-- - [practical solutions](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2_solutions.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2_solutions.pdf)) -->
  <!-- - [live script](live_scripts/DESeq2_script.R)  -->
  
  -  [Linear Models in R and DESeq2 (Worksheet)](Bulk_RNAseq_Course_Base/Markdowns/07_Linear_Models.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/07_Linear_Models.pdf))  
  - [DESeq2 results extraction cheatsheet](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/DESeq2_results_cheatsheet.pdf)


13:00 - 14:00 Lunch


### Day 3

**Trainers:** Ash, Ian, Jiayin, Ramy

9:30 - 9:45 - [Recap of Day 1 and 2](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/Analysis_of_RNA-seq_data_day3recap.pdf) - Ash

9:45 - 11:00 [Annotation of RNA-seq results](Bulk_RNAseq_Course_Base/Markdowns/09_Annotation_Slides.html) - Ash. 
  - [Practical](Bulk_RNAseq_Course_Base/Markdowns/09_Annotation.html)
  <!-- - [live script](live_scripts/visualisation.R) -->

11:00 - 12.30 [Visualisation differential expression results](Bulk_RNAseq_Course_Base/Markdowns/10_Data_Visualisation_Slides.html) - Jiayin

   - [Practical](Bulk_RNAseq_Course_Base/Markdowns/10_Data_Visualisation.html)
  <!-- - [practical solutions](Bulk_RNAseq_Course_Base/Markdowns/10_Data_Visualisation_solutions.html) -->
  <!-- - [live script](live_scripts/visualisation.R) -->

12.30 - 13.30 Lunch

13.30 - 16:30  [Gene-set testing](Bulk_RNAseq_Course_Base/Markdowns/11_Gene_set_testing_introduction.html) - Ian
   - [Practical (html)](Bulk_RNAseq_Course_Base/Markdowns/11_Gene_set_testing.html) [(pdf)](Bulk_RNAseq_Course_Base/Markdowns/11_Gene_set_testing.pdf)
  <!-- - [Practical solutions (html)](Bulk_RNAseq_Course_Base/Markdowns/11_Gene_set_testing_solutions.html) (pdf)](Bulk_RNAseq_Course_Base/Markdowns/11_Gene_set_testing_solutions.pdf) -->
  <!-- - [Live Script](live_scripts/gene_set_testing.R) -->


## Data, software and materials

The lecture slides and other source materials, including R code and
practical solutions, can be found in the course's [Github
repository](https://github.com/bioinformatics-core-shared-training/Bulk_RNASeq_Course_March23)

The full data used in the course can be [downloaded from dropbox](https://www.dropbox.com/sh/sz44que2vha44xw/AABISE1DdBSS6s_zLoW1vuCGa?st=z95zfjyg&dl=0). 

Instructions to install software are available from the "[Software installation instructions](Bulk_RNAseq_Course_Base/Markdowns/setup.html)" page.

## Extended materials

The [Extended Materials](Extended_index.md) contain extensions to some of the
sessions and additional materials, including instruction on downloading and
processing the raw data for this course, a link to an excellent R course, and
where to get further help after the course.

## Additional Resources

* [Bioconductor for relevant R packages](https://bioconductor.org/)
* [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)  
* [RNAseq Workflow](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)  
* [RStudio CheatSheets](https://rstudio.com/resources/cheatsheets/)

## Acknowledgements

This course is based on the course [RNAseq analysis in
R](http://combine-australia.github.io/2016-05-11-RNAseq/) prepared by [Combine
Australia](https://combine.org.au/) and delivered on May 11/12th 2016 in
Carlton. We are extremely grateful to the authors for making their materials
available; Maria Doyle, Belinda Phipson, Matt Ritchie, Anna Trigos, Harriet
Dashnow, Charity Law.

![](Bulk_RNAseq_Course_Base/images/combine_banner_small.png)

The materials have been rewritten/modified/corrected/updated by various
contributors over the past 5 years including:

Abigail Edwards
Ashley D Sawle
Chandra Chilamakuri
Dominique-Laurent Couturier
Guillermo Parada González
Hugo Tavares
Jon Price
Mark Dunning
Mark Fernandes
Oscar Rueda
Sankari Nagarajan
Stephane Ballereau
Tom Smith
Zeynep Kalender Atak

Apologies if we have missed anyone!
