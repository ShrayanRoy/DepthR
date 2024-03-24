# DepthR

Collection of functions to estimate depth from Natural Images.

Install the package using the following R code.

```
   devtools::install_github("ShrayanRoy/DepthR")
```

## Instructions

* Package images are stored in `imgDepthR` object. Use `imgDepthR[[i]]` with $1 \leq i \leq 9$ to use these images.

* Segmentation data produced by **SAM** for these images are stored as `segmeniData` and `bboxiData` object with $1 \leq i \leq 9$
