# PET

`PET.cpp` is a small C++/ROOT demonstration of tomographic reconstruction inspired by positron emission tomography (PET). It studies how the angular sampling of a detector affects a sinogram and the images reconstructed from it. It is an educational numerical example, not a clinical PET simulator or a tool for interpreting medical images.

## What it does

The `PET()` macro builds a synthetic 2D activity profile from simple shapes (circles, a rectangle, and ellipses) on a 256 × 256 grid covering x and y from -1 to 1. It then:

1. Computes a sinogram by numerically integrating the profile along projection lines. The projection angle spans 0 to π, and each projection has 256 position bins.
2. Reconstructs an image by directly backprojecting the sinogram.
3. For each projection angle, computes a 1D FFT along the detector-position axis of the sinogram, applies a selected frequency-domain weight, and performs the corresponding 1D inverse FFT along that same axis. No 2D FFT is used.
4. Backprojects the filtered projections to produce a filtered reconstruction.
5. Saves plots of the source profile, sinogram, reconstructions, real and imaginary FFT components, inverse FFT, and filter.

The macro can optionally replace one sinogram channel with a chosen value or add non-negative Gaussian noise to every sinogram bin. Its filters are custom frequency weights implemented in `TSinogram.cpp`; filter 4, the default, is a triangular weight over the frequency-bin index.

## Projection and backprojection

For each angle `φ = iπ/Nφ`, the forward projector samples parallel lines at detector offsets `s`. It parameterizes each line with a coordinate `z` along the line:

```text
x = s cos(φ) - z sin(φ)
y = s sin(φ) + z cos(φ)
```

The sinogram value at `(s, φ)` is a numerical line integral: the code samples the source profile along `z` and accumulates `profile(x, y) × Δs × Δz`. The histogram's first axis is intended to represent detector offset and its second axis represents angle `φ`; the angles are sampled from 0 up to, but not including, π.

The backprojector then visits each output image-grid point `(x, y)`, looks up a sinogram value for each sampled angle, and sums those values. It uses the histogram bin containing the lookup coordinate (no interpolation) and divides the sum by the number of detector-position bins.

**Implementation caveat:** the forward-projection loop fills successive detector bins while its computed offset `s` decreases from positive to negative, opposite the histogram's increasing first-axis labels. Also, the backprojector looks up `z = -x sin(φ) + y cos(φ)` as the detector offset, whereas the forward projector computes the perpendicular offset as `s = x cos(φ) + y sin(φ)`. Since `z` is along the projection line, not perpendicular to it, the backprojector does not select the forward projector's matching line coordinate. This describes the current implementation.

## Requirements and running

Run from the repository directory with [ROOT](https://root.cern/) installed. ROOT's ACLiC compiles the macro and its included source files:

```sh
root -l -q 'PET.cpp+(1, 4, false, false, 0.0005)'
```

The first argument selects how many angular samplings to run. The available projection-bin counts, in processing order, are `64, 2, 4, 8, 16, 32, 128`. For example, run all seven:

```sh
root -l -q 'PET.cpp+(7)'
```

To run with a sinogram channel set to zero and noise enabled:

```sh
root -l -q 'PET.cpp+(1, 4, true, true, 0.0005, 0.0)'
```

## `PET()` arguments

Arguments are positional; omitted arguments use the defaults shown:

| Argument | Default | Description |
| --- | --- | --- |
| `nPhiBins` | `1` | Number of angular samplings to process, from 1 to 7. The macro uses the first `nPhiBins` entries in the sequence above. |
| `ifilter` | `4` | Frequency weight to apply before the inverse FFT. Values 0–4 select the custom weights in `TSinogram::EvalFilter`; other values leave the spectrum unchanged. |
| `killSinoChannels` | `false` | If true, replace sinogram channel 121 at every angle with `killValue`. |
| `smearSino` | `false` | If true, add the absolute value of a Gaussian random draw to each sinogram bin. |
| `sigma` | `0.0005` | Standard deviation of the added Gaussian noise. Used only when `smearSino` is true. |
| `killValue` | `0.0` | Replacement value for channel 121. Used only when `killSinoChannels` is true. |
| `noTitle` | `true` | Hide ROOT plot titles and adjust plot margins. |

For example, to run the first four angular samplings with filter 3:

```sh
root -l -q 'PET.cpp+(4, 3)'
```

## Output

PNG plots are written to the current directory. For each angular sampling, the macro saves canvases named `PET_<N>iphi.png`, `PET_FFT_<N>iphi.png`, and `PET_FFT_proj_<N>iphi.png`, where `<N>` is the projection-bin count. It also writes individual diagnostic plots whose names start with `can_`, plus a filter plot such as `filter_Filter4.png`.

The reconstruction is a simple demonstration of projection, FFT filtering, and backprojection. It does not model PET event detection or provide a validated physical detector model.
