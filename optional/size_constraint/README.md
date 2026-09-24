# Optional size-constrained variant

This folder contains an experimental extension of HSSS. The original
implementation in `core/` is unchanged and remains the primary method.

The extension adds a soft penalty to the candidate-merge score:

```text
lambda_size * median_D * max(0, merged_area / target_area - size_tau)^2
```

The default values are `LambdaSize = 2` and `SizeTau = 1.25`. They were
selected as a balanced setting for `K = 400`, rather than as universal
parameters for every image or target count.

Run `demo_size_constraint.m` for a single-image example. Use
`hsss_segment_size` in the same way as `hsss_segment`, with the two extra
option fields shown above.

For robustness with arbitrary edge maps, the public wrapper detects a
stalled size-constrained hierarchy and uses the original merge rule to finish
the remaining requested count. The returned `info.UsedOriginalFallback` flag
reports whether this happened. The underlying experimental merge function is
unchanged.

In our local 500-image BSDS500 check at `K = 400` (strict mean over all
available ground-truth segmentations), the size constraint changed mean ASA
from 0.952411 to 0.959441 and mean BR from 0.907289 to 0.910883. Exact values
depend on the edge maps and evaluation protocol.
