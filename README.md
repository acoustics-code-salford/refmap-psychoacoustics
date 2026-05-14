[![REFMAP logo](assets/proj/horizontal_REFMAP_FINAL_LOGO-black-bg.png)](https://www.refmap.eu/)
# refmap-psychoacoustics
RefMap project task 4.3: Psychoacoustics of sound from unmanned aircraft systems (UAS) / unmanned aerial vehicles (UAVs), and urban air mobility (UAM)

## Sound quality metrics
### MATLAB
Verified implementations of the standardised Sottek Hearing Model sound quality metrics (SQMs) for loudness, tonality and roughness [1, 2] are provided in the [ECMA 418-2](src/mlab/ECMA_418-2) folder. Versions of these implementations have also been incorporated into the [SQAT](https://github.com/ggrecow/SQAT) (Sound Quality Analysis Toolbox) repository.

An implementation of the aural detectability metric [3, 4] is available in the [Metrics](src/mlab/Metrics) folder.

### Python
The MATLAB implementations of the Sottek Hearing Model SQMs have been translated into a python package, [sottek-hearing-model](https://github.com/mlotinga/sottek-hearing-model).

A python version of the aural detectability metric is also available in the [metrics](src/refmap_psychoacoustics/metrics) folder.

# References
[1] Ecma International. (2025). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 4th Edition/June 2025). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf).

[2] Lotinga, M. J. B., Torjussen, M, & Felix Greco, G. (2025). Verified implementations of the Sottek psychoacoustic Hearing Model standardised sound quality metrics (ECMA-418-2 loudness, tonality and roughness). 11th Convention of the European Acoustics Association (Forum Acusticum / Euronoise), 23-26 June 2025, Malaga, Spain. [https://www.researchgate.net/publication/392904348](https://www.researchgate.net/publication/392904348).

[3] Fidell, S., Pearsons, K. S., & Bennett, R. (1974). Prediction of aural detectability of noise signals. Human Factors, 16(4), 373–383. [DOI: 10.1177/001872087401600405](https://doi.org/10.1177/001872087401600405).

[4] Rizzi, S. A., Christian, A., Letica, S. J., & Lympany, S. V. (2025). Annoyance model assessments of urban air mobility vehicle operations. Journal of Aircraft, 0(0), 1–16. [DOI: 10.2514/1.C038188](https://doi.org/10.2514/1.C038188).

# Ownership and licensing
This repository is maintained by [mlotinga](https://github.com/mlotinga)

The code developed under this project is licensed for open-access use as defined in [the licensing page](https://github.com/acoustics-code-salford/refmap-psychoacoustics/blob/main/LICENSE). Please be aware that sub-components used in the code may have alternative licensing applicable - the relevant licence information is included alongside such components where required to comply with the applicable licence conditions.

# How to cite code from this repository
If you use, reuse or adapt any original code from this repository, please include an appropriate citation, including the following information:
> Lotinga, M.J.B. (2026). refmap-psychoacoustics. GitHub. [https://github.com/acoustics-code-salford/refmap-psychoacoustics](https://github.com/acoustics-code-salford/refmap-psychoacoustics)

# Disclaimer
This repository contains 'work in progress' code for a research project. It is not a completed package or code library. As per the licensing information, the code is published WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
