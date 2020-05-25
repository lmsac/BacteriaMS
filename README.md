# BacteriaMS
Identification and assessment for MALDI-MS based bacterial analysis. It combines new measures for spectra similarity and a novel bootstrapping assessment that could lead to more reliable bacterial typing.


## Requirement
* [R] (https://www.r-project.org/). As an alternative, the latest version of [Microsoft R Open] (https://mran.microsoft.com/open/) should be fine.

* [RStudio] (https://www.rstudio.com/) is recommended but optional.


## Database
BacteriaMS uses peak lists for bacteria identification. A peak list is stored in a space-separated text (.txt) file including the *m/z* and *intensity* of all the peaks from a mass spectrum of a strain of bacteria with the first column as *m/z* and the second column as *intensity*, e.g.:

```
4376.517617 6369.468
4545.36698 14923.287008
4616.504007 3381.736378
4778.424564 6757.753407
4830.992604 10117.37512
……
9899.339424 4964.88939
10159.44407 1130.018489
10323.81304 5016.852179
10707.91308 1477.481354
11240.44102 1693.14288
```

Some example peak lists are provided in this repo.

In addition, below are some accessible databases where users can download MALDI-TOF spectra of bacteria, under their terms of usage.

* European Consortium of Microbial Resources Centres (EM-baRC): [http://www.embarc.eu/deliverables/EMbaRC_D.JRA2.1.4_D15.33_MALDI-TOF_DB.pdf]

* Robert Koch-Institute: [http://www.microbe-ms.com/]

* The Public Health Agency of Sweden: [http://spectra.folkhalsomyndigheten.se/spectra/]

The spectra are only compatible with the analysis software Biotyper (Bruker), but can be converted to peak list files using Bruker's software.


## Usage
If you have everything installed, you can run identification for a sample spectrum as follows:

1. Run [main.R] (main.R) to load the functions.

2. Load the reference peak lists.
```R
setwd(REF_PATH) 
# REF_PATH is the path of the directory that contains the reference peak list files.

reference.peaklists = local({
  filenames = list.files(pattern = '.txt')
  peaklists = lapply(filenames, read.table)
  names(peaklists) = sub('.txt', '', filenames)
  lapply(peaklists, normalize.intensity)
})
```

3. Load the sample peak list.
```R
sample.peaklist = read.table(SAMPLE_PATH)
# SAMPLE_PATH is the path of the sample peak list file.
sample.peaklist = normalize.intensity(sample.peaklist)
```


4. Search database.
```R
result = search.datebase.bootstrap(
  sample.peaklist,
  reference.peaklists,
  tolerance = 2000, # Unit: ppm
  method = c('eu', 'ieu', 'cosine'), # Similarity scoring method(s)
  bootstrap.times = 1000
)

result = lapply(result, function(res) {
  list(
    res$sample,
    genus = get.bootstrap.score(
      res$bootstrap, 
      function(x) sapply(strsplit(x, ' '), function(x) paste(x[1]))
    ),
    species = get.bootstrap.score(
      res$bootstrap, 
      function(x) sapply(strsplit(x, ' '), function(x) paste(x[1], x[2]))
    )      
  )
})

# Print the results (cosine as an example)
result$cosine[[1]]    # Print the matches ranked by the cosine scores.
result$cosine$species # Print the confidence scores at the species level.
result$cosine$genus   # Print the confidence scores at the species level.
```

The script [bootstrapping.R] (scripts/bootstrapping.R) runs a cross-validation test on a set of spectra. Each spectrum in the dataset was selected as the sample spectrum, and the other spectra were used to construct a reference database, respectively. Plots receiver operating characteristic (ROC) curves for identification.


## Publications
Yang, Y., Lin, Y., Chen, Z., Gong, T., Yang, P., Girault, H., Liu, B., Qiao, L. Bacterial whole cell typing by mass spectra pattern matching with bootstrapping assessment. *Anal Chem* **89**, 12556–12561 (2017).

## License
BacteriaMS is distributed under a BSD license. See the LICENSE file for details.


## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn  
