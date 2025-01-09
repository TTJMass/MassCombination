# TTJetMassCombination

This repo centralizes all inputs, the combination setup, plotting and fitting scripts, and the results.

## How to install
To clone the repo with all submodules do
```bash
git clone --recursive https://github.com/TTJMass/MassCombination.git
```
or
```bash
git clone https://github.com/TTJMass/MassCombination.git
cd MassCombination
git submodule update --init --recursive
```

## Inputs
### Theory
The theoretical calculations are stored in `powheg_generations/` and `powheg_generations_8tev/` respectively. They are synched from `/eos` and can be downloaded via `downloadTheory.sh`.
The files are managed via [Git-LFS](https://git-lfs.com/).

### Experiment
