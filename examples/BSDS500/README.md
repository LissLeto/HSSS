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
