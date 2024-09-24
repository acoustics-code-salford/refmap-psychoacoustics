# ECMA-418-2:2022 - Hearing Model of Sottek tonality, loudness and roughness sound quality metrics validation of MATLAB implementation
Validation of the MATLAB algorithms has been undertaken by comparison with outputs calculated using HEAD Acoustics ArtemiS v15.6 software.

The audio signals used for the validation comprise:
1. reference calibration signal for tonality and loudness: 1 kHz sinusoid at 40 dB sound pressure level (5 seconds, mono)
1. reference calibration signal for roughness: 1 kHz sinusoid modulated at 70 Hz (modulation factor 1), at 60 dB sound pressure level (5 seconds, mono)
1. binaural audio recording of a 'busy city street' environment, extracted from [EigenScape](https://zenodo.org/doi/10.5281/zenodo.1012808) database and reused under Creative Commons Attribution 4.0 licence

The reference calibration signals 1 and 2 were generated using [acousticHMSGenerateRefSignals.m](validation/ECMA-418-2_2022/mlab/acousticHMSGenerateRefSignals.m).

The reference ArtemiS results are included in [reference folder](validation/ECMA-418-2_2022/reference).

Calculated sound quality values and reference comparison figures were generated using the [validation script](validation/ECMA-418-2_2022/mlab/acousticHMSValidation.m). The full set of comparison plots can be displayed by running the validation script; the selection presented below is sufficient to validate all capabilities of the algorithms.

Signal 1 is unmodulated, which yields 0 asper roughness, so time-dependent and specific roughnesses are not displayed.

# Tonality
## Time-dependent tonality

[![Signal 1 time-dependent tonality](validation/ECMA-418-2_2022/results/tonalHMSTDepSine1kHz40dB.pdf)]

[![Signal 3 time-dependent tonality](validation/ECMA-418-2_2022/results/tonalHMSTDepBusySt.pdf)]

## Time-dependent specific tonality

[![Signal 1 time-dependent specific tonality](validation/ECMA-418-2_2022/results/tonalHMSSpecTDepSine1kHz40dB.png)]

[![Signal 3 time-dependent specific tonality](validation/ECMA-418-2_2022/results/tonalHMSSpecTDepBusySt.png)]

## Time-aggregated specific tonality

[![Signal 3 time-aggregated specific tonality](validation/ECMA-418-2_2022/results/tonalHMSSpecTAggBusySt.pdf)]

## Overall tonality

[![All signals overall tonality](validation/ECMA-418-2_2022/results/tonalHMSsingles.pdf)]

# Loudness
## Time-dependent loudness

[![Signal 1 time-dependent loudness](validation/ECMA-418-2_2022/results/loudHMSTDepSine1kHz40dB.pdf)]

[![Signal 3 time-dependent loudness](validation/ECMA-418-2_2022/results/loudHMSTDepBusySt.pdf)]

## Time-dependent specific loudness

[![Signal 1 time-dependent specific loudness](validation/ECMA-418-2_2022/results/loudHMSSpecTDepSine1kHz40dB.png)]

[![Signal 3 time-dependent specific loudness](validation/ECMA-418-2_2022/results/loudHMSSpecTDepBusySt.png)]

## Time-aggregated specific loudness

[![Signal 3 time-aggregated specific loudness](validation/ECMA-418-2_2022/results/loudHMSSpecTAggBusySt.pdf)]

## Time-dependent specific binaural loudness

[![Signal 3 time-dependent specific binaural loudness](validation/ECMA-418-2_2022/results/loudHMSSpecTDepBinBusySt.png)]

# Roughness
## Time-dependent roughness

[![Signal 2 time-dependent roughness](validation/ECMA-418-2_2022/results/roughHMSTDepSine1kHz70Hz60dB.pdf)]

[![Signal 3 time-dependent roughness](validation/ECMA-418-2_2022/results/roughHMSTDepBusySt.pdf)]

## Time-dependent specific roughness

[![Signal 2 time-dependent specific roughness](validation/ECMA-418-2_2022/results/roughHMSSpecTDepSine1kHz70Hz60dB.png)]

[![Signal 3 time-dependent specific roughness](validation/ECMA-418-2_2022/results/roughHMSSpecTDepBusySt.png)]

## Time-aggregated specific roughness

[![Signal 3 time-aggregated specific roughness](validation/ECMA-418-2_2022/results/roughHMSSpecTAggBusySt.pdf)]

## Time-dependent specific binaural roughness

[![Signal 3 time-dependent specific binaural roughness](validation/ECMA-418-2_2022/results/roughHMSSpecTDepBinBusySt.png)]
