# Bundled BSDS500 examples

This folder contains five small examples used by the single-image demos:

- 3063
- 22090
- 104055
- 124084
- 228076

Each example uses the following naming convention:

```text
<id>.jpg       input image
<id>_gt.mat    BSDS500 groundTruth cell array
<id>_rcf.png   precomputed RCF edge map
```

`demo_eval.m` runs HSSS on image 3063 at K=400, evaluates BR and IBR against
all annotations in `3063_gt.mat`, and visualizes the result using annotation 4.
`22090_labels.csv` is retained as an additional saved example result.

`selected_image_ids.txt` lists the 27 BSDS500 image IDs used for an additional
selection. The corresponding dataset images and annotations are not included
in this repository; obtain BSDS500 separately and match files by image ID.
BSDS500 remains subject to its own distribution terms.
