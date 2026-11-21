# Layout optimization: questions and requests for the simulation group

Context: we re-score arbitrary detector layouts against the existing CORSIKA surface-particle output of `official_simulations/datasets/v1/volcano_valley` (exact re-mapping reproduces the official hits with zero mismatches), then reconstruct with LlhamaFit (`TRIGGER_5x1`, `PlaneTiming`). Signal is ν_τ weighted with TamboSim v1.0.1 OneWeights and the TamboMakie reference flux. The background gate uses the proton sets with the all-particle cosmic-ray flux. The optimizer is limited by four things in the inputs.

## 1. The upstream ν_τ channel is undersampled (request: enriched production)

The ν_τ sets mix two injection channels: `ForcedNeutrinoInteractionPoint` (interaction forced in the detector region; OneWeight includes the interaction probability) and `UpstreamNeutrinoInteractionPoint` (ν→τ conversion during Earth transit via TauRunner; OneWeight = 1/surface pdf, about 10⁷ times larger). Upstream events carry most of the predicted rate at 0.1–10 PeV but are very rare in the MC:

| set | forced events | upstream events | upstream share of predicted rate |
|---|---|---|---|
| nu_tau/00010 | 32,028 | 2 | ~90% |
| nu_tau/00011 | 31,805 | 4 | ~50% |
| nu_tau/00012 | 31,849 | 186 | ~58% |
| nu_tau/00013 | 3,097 | 52 | ~12% |
| nu_tau/00014 | 3,040 | 130 | ~12% |
| nu_tau/00019 | 160,943 | 0 | 0% |

A single upstream event in nu_tau/00010 (run 1, event 2198, 0.87 PeV, θ = 99°) carries 88% of that set's rate. For now we optimize on the forced channel only. The bootstrap relative error on the upstream-channel rate is 40–130%.

- **Request:** a production that oversamples the upstream channel, or a forced-conversion option for Earth transit, in the 0.1–10 PeV range.
- **Question:** is the upstream OneWeight (1/surface pdf, with no survival or conversion factor) the intended normalization for TamboSim v1.0.1?

## 2. Proton statistics limit the background gate

After the quality cut (fitted, θ̂ ≤ 90°, u3 > 0.1), well-performing layouts leak 0–2 proton MC events in total. Whether a layout leaks 0, 1 or 2 is mostly chance at this sample size, but each leak changes the ranking.

- **Request:** more proton (and Fe) showers, especially at high energy, so a leak rate can be quoted with a useful upper limit.

## 3. What is the deployable slope limit?

`official125` has 258 of its 917 units on 35–48° terrain. Its per-zone counts on terrain of 35° or less are identical to `TamboSim.place_detector_units(spacing = 125 m)`, which uses a 35° default. Those steep units raise the training-split signal rate by about 70% (0.0215 vs 0.0125 per year for the forced channel), so the slope limit strongly affects any optimized layout.

- **Question:** is ~48–50° deployable, or should optimized layouts respect 35°?

## 4. Conventions to confirm

- Local frame: +x = east and +y = north? Zone names such as "east_facing" assume this.
- Flux: the TamboMakie reference flux (`flux_nu` ∝ E^-2.52, `flux_CR` = 1.8e4 E^-2.7). Is this the agreed benchmark?
- Are ν̄_τ included in the ν_τ sets, or do rates need a factor of 2?
