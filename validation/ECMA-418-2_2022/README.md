# ECMA-418-2:2022 - Hearing Model of Sottek tonality, loudness and roughness sound quality metrics: Validation of MATLAB implementation
Validation of the MATLAB algorithms has been undertaken by comparison with outputs calculated using HEAD Acoustics ArtemiS v15.6 software.

The audio signals used for the validation comprise:
1. reference calibration signal for tonality and loudness: 1 kHz sinusoid at 40 dB sound pressure level (5 seconds, mono)
1. reference calibration signal for roughness: 1 kHz sinusoid modulated at 70 Hz (modulation factor 1), at 60 dB sound pressure level (5 seconds, mono)
1. binaural audio recording of a 'busy city street' environment, extracted from [EigenScape](https://zenodo.org/doi/10.5281/zenodo.1012808) database and reused under Creative Commons Attribution 4.0 licence (30 seconds, 2-channel binaural)

The reference calibration signals 1 and 2 were generated using [acousticHMSGenerateRefSignals.m](mlab/acousticHMSGenerateRefSignals.m).

The reference ArtemiS results are included in [reference folder](reference).

Calculated sound quality values and reference comparison figures were generated using the [validation script](mlab/acousticHMSValidation.m). The full set of comparison plots can be displayed by running the validation script; the selection presented below is sufficient to validate all capabilities of the algorithms.

Signal 1 is unmodulated, which yields 0 asper roughness, so time-dependent and specific roughnesses are not displayed.

# Tonality
## Time-dependent tonality

![Signal 1 time-dependent tonality](results/tonalHMSTDepSine1kHz40dB.png | width=150)

![Signal 3 time-dependent tonality](results/tonalHMSTDepBusySt.png)

## Time-dependent specific tonality

![Signal 1 time-dependent specific tonality](results/tonalHMSSpecTDepSine1kHz40dB.png)

![Signal 3 time-dependent specific tonality](results/tonalHMSSpecTDepBusySt.png)

## Time-aggregated specific tonality

![Signal 3 time-aggregated specific tonality](results/tonalHMSSpecTAggBusySt.png)

## Overall tonality

![All signals overall tonality](results/tonalHMSsingles.png)

# Loudness
## Time-dependent loudness

![Signal 1 time-dependent loudness](results/loudHMSTDepSine1kHz40dB.png | width=150)

![Signal 3 time-dependent loudness](results/loudHMSTDepBusySt.png)

## Time-dependent specific loudness

![Signal 1 time-dependent specific loudness](results/loudHMSSpecTDepSine1kHz40dB.png)

![Signal 3 time-dependent specific loudness](results/loudHMSSpecTDepBusySt.png)

## Time-aggregated specific loudness

![Signal 3 time-aggregated specific loudness](results/loudHMSSpecTAggBusySt.png)

## Time-dependent specific binaural loudness

![Signal 3 time-dependent specific binaural loudness](results/loudHMSSpecTDepBinBusySt.png)

## Overall loudness

![All signals overall loudness](results/loudHMSsingles.png)

# Roughness
## Time-dependent roughness

![Signal 2 time-dependent roughness](results/roughHMSTDepSine1kHz70Hz60dB.png | width=150)

![Signal 3 time-dependent roughness](results/roughHMSTDepBusySt.png)

## Time-dependent specific roughness

![Signal 2 time-dependent specific roughness](results/roughHMSSpecTDepSine1kHz70Hz60dB.png)

![Signal 3 time-dependent specific roughness](results/roughHMSSpecTDepBusySt.png)

## Time-aggregated specific roughness

![Signal 3 time-aggregated specific roughness](results/roughHMSSpecTAggBusySt.png)

## Time-dependent specific binaural roughness

![Signal 3 time-dependent specific binaural roughness](results/roughHMSSpecTDepBinBusySt.png)

## Overall roughness

![All signals overall roughness](results/roughHMSsingles.png)