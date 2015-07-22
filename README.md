# Dynamic-Binary-Choice-Model
---
This repository contains the R code that implements the dynamic discrete choice model we estimate the paper [Too Proud to Stop](http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2465840). 

**Note**: _We do not publish the data here to which we fit the model in our paper. You may request the data from us directly._

## Content

There are three essential files we provide here

+ functions.R: This files contains all the necessary functions we need. This is where you will find how we implemented our model.
+ estimate_models.R: This file sources functions.R and fits a model to the **pooled** data.
+ estimate_individual_models.R: Does the same as estimate_models.R, but fits the model to the per-subject data.
