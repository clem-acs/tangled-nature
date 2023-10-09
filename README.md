# Tangled Nature Model in Rust

The Tangled Nature model is a theoretical framework used to understand the evolution and dynamics of species in complex ecosystems [Jensen 2018](https://iopscience.iop.org/article/10.1088/1361-6404/aaee8f/meta). This Rust implementation offers a fast and efficient simulation of the model, capturing the intricate relationships between individuals and species.

## Table of Contents

- [Overview](#overview)
- [Installation and Usage](#installation-and-usage)
- [Model Parameters](#model-parameters)
- [File Outputs](#file-outputs)

## Overview

This codebase provides an implementation of the Tangled Nature model, which seeks to understand:

- How species emerge, evolve, and go extinct.
- The interactions between species and their environment.
- The dynamics of population sizes and species distributions over time.

## Installation and Usage

1. Ensure you have the Rust programming language and Cargo (the Rust package manager) installed on your machine.
2. Clone this repository.
3. To run the simulation:

```bash
cargo run -- [seed] [j_file_path] [mutation_rate]
```

- `seed`: (Optional) Seed for the random number generator. Defaults to the current system time.
- `j_file_path`: (Optional) Path to an existing J matrix file. If not provided, a new J matrix will be generated and saved.
- `mutation_rate`: (Optional) Rate of mutation. Defaults to `0.001`.

## Model Parameters

The model uses several parameters to define its behavior:

- `L`: Length of the genome.
- `GENOMES`: Total number of different species.
- `N_INIT`: Initial total population size.
- `GENERATIONS`: Length of the simulation in generations.
- `THETA`: For J matrix initialization.
- `P_KILL`: Probability of killing an individual.
- `R` and `W`: Parameters for calculating reproduction probability.

## File Outputs

The simulation produces several output files detailing the state and progress of the simulation:

- `j_[seed].txt`: Contains the interaction matrix J.
- `simlog/species_genome_[seed]_[mutation_rate].txt`: Logs the genomes of species at specific intervals.
- `simlog/species_population_[seed]_[mutation_rate].txt`: Logs the population sizes of species at specific intervals.
