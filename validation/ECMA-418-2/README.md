# ECMA-418-2:2022 - Sottek Hearing Model tonality, loudness and roughness sound quality metrics: Validation of MATLAB implementation
Validation of the MATLAB algorithms has been undertaken by comparison with outputs calculated using HEAD Acoustics ArtemiS v15.6 software.

The audio signals used for the validation comprise:
1. reference calibration signal for tonality and loudness: 1 kHz sinusoid at 40 dB sound pressure level (5 seconds, mono)
1. reference calibration signal for roughness: 1 kHz sinusoid modulated at 70 Hz (modulation factor 1), at 60 dB sound pressure level (5 seconds, mono)
1. binaural audio recording of a 'busy city street' environment (30 seconds, 2-channel binaural)

The reference calibration signals 1 and 2 were generated using [acousticSHMGenerateRefSignals.m](mlab/acousticSHMGenerateRefSignals.m).

The binaural audio recording was extracted from the [EigenScape](https://zenodo.org/doi/10.5281/zenodo.1012808) database (Green et al., [2017](https://doi.org/10.3390/app7111204),  licenced under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0)).

The reference ArtemiS results are included in the [reference folder](reference).

Calculated sound quality values and reference comparison figures were generated using the [validation script](mlab/acousticSHMValidation.m), which calls the algorithms from the [ECMA-418-2](../../src/mlab/ECMA_418-2) folder. The full set of comparison plots can be displayed by running the validation script; the selection presented below is sufficient to validate all capabilities of the algorithms.

Signal 1 is unmodulated, which yields 0 asper roughness, so time-dependent and specific roughnesses are not displayed for this signal.

# Tonality
## Time-dependent tonality

<img src='results/tonalSHMTDepSine1kHz40dB.png' width=400>

![Signal 3 time-dependent tonality](results/tonalSHMTDepBusySt.png)

## Time-dependent specific tonality

![Signal 1 time-dependent specific tonality](results/tonalSHMSpecTDepSine1kHz40dB.png)

![Signal 3 time-dependent specific tonality](results/tonalSHMSpecTDepBusySt.png)

## Time-aggregated specific tonality

![Signal 3 time-aggregated specific tonality](results/tonalSHMSpecTAggBusySt.png)

## Overall tonality

<img src='results/tonalSHMsingles.png' width=400>

# Loudness
## Time-dependent loudness

<img src='results/loudSHMTDepSine1kHz40dB.png' width=400>

![Signal 3 time-dependent loudness](results/loudSHMTDepBusySt.png)

## Time-dependent specific loudness

![Signal 1 time-dependent specific loudness](results/loudSHMSpecTDepSine1kHz40dB.png)

![Signal 3 time-dependent specific loudness](results/loudSHMSpecTDepBusySt.png)

## Time-aggregated specific loudness

![Signal 3 time-aggregated specific loudness](results/loudSHMSpecTAggBusySt.png)

## Time-dependent specific binaural loudness

![Signal 3 time-dependent specific binaural loudness](results/loudSHMSpecTDepBinBusySt.png)

## Overall loudness

<img src='results/loudSHMsingles.png' width=400>

# Roughness
## Time-dependent roughness

<img src='results/roughSHMTDepSine1kHz70Hz60dB.png' width=400>

![Signal 3 time-dependent roughness](results/roughSHMTDepBusySt.png)

## Time-dependent specific roughness

![Signal 2 time-dependent specific roughness](results/roughSHMSpecTDepSine1kHz70Hz60dB.png)

![Signal 3 time-dependent specific roughness](results/roughSHMSpecTDepBusySt.png)

## Time-aggregated specific roughness

![Signal 3 time-aggregated specific roughness](results/roughSHMSpecTAggBusySt.png)

## Time-dependent specific binaural roughness

![Signal 3 time-dependent specific binaural roughness](results/roughSHMSpecTDepBinBusySt.png)

## Overall roughness

<img src='results/roughSHMsingles.png' width=400>
