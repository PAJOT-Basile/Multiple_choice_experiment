# Vcf2Colony

This folder is used to transform a vcf format to a colony input format (see colony documentation available [here](https://www.zsl.org/about-zsl/resources/software/colony) to see the format). This folder contains two directories: `Data` and `Scripts`.

The `Data` directory contains the details about the conda environment and the input file to be filled up for the transformations to work.
The `Input.yaml` is the file to modify. Each line is a variable that is used to run the format transformation. They can be modified in accordance with the colony input format.

To execute the pipeline, run:
```
sbatch ./Scripts.Vcf2Colony.sh
```