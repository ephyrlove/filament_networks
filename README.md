# filament_networks

## /data
This directory contains a lot of data and so is not in the git. Talk to me or Dr. Maroulas for access to the data.

## /fig
This folder contains the output figures, mostly from R. These figures were used in ERL's dissertation and likely to be used in manuscript with supergroup.

## /doc
These are some old doccumentation files. All recent work is is done in the latex directory under src/latex/.

## /src

### /cpp
This contains the ripser src code. You will need to compile for your machine from https://github.com/Ripser/ripser.

### /python
This contains a brief project looking at applying a CNN to the filament images. Very little productive work in here.

### /bash
There are some bash scripts that make it easy to parallely compute rips diagrams with ripser on a linux machine.

### /latex
This contains the original manuscript code. The complete history on just the manuscript is available in branch **main**

### /R
This is the meat of this repository, so I will outline all files:

* classify_cells_from_skeletons.R
  * This is the most recent and most important file. This file contains code to go from Andreas's skeletonized (preprocessed) images of cells to a subsampled rips diagram, summary stats, SVM classifier, and visualization of clustering results.
* classify_cells.R
   * This is the essentially the same code as above, but inteded for a raw cell with no preprocessing. Therefore, this file includes preprocessing code to enhance filaments. This file contains code to go from Andreas's skeletonized (preprocessed) images of cells to a subsampled rips diagram, summary stats, SVM classifier, and visualization of clustering results.
* classify_was.R
* compare_old_new.R
* dpc.R
* filament_images.R
* filament_images_ripser.R
* ko_2ko_wt_micro.R
* mutant_vs_wildtype_ANsum.R
  * This looks specifically at Andreas's summary data, though this data became important to many analyses.
* mutant_vs_wildtype_micro_noThresh_6f_newData.R
  * The files with \*6f\* were just my rewriting of the subsampling code to be able to take the number of subsamples as an argument.
* mutant_vs_wildtype_micro_noThresh_6f.R
* mutant_vs_wildtype_micro_noThresh.R
* mutant_vs_wildtype_micro.R
* persistance_images.R
* plot_several_cells.R
* prepare_cells_for_quiz.R
* simplicial_complex demo.R
* synth_eda_pd_alternatives.R
* synth_eda.R
* synth_eda_subsampled.R
