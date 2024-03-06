# Chemical Reaction Compound Atom-to-atom Mapping

A chemical reaction can be defined as a transformation of one set of chemical compounds to another. The correct mapping
of the rearrangement of chemical compound atoms during the transformation is essential for capturing the essence of a
chemical reaction. This task, widely known as atom-to-atom mapping, has proven quite challenging. Nevertheless, novel
approaches are being published frequently. The main goal of this repository is to curate existing open-source chemical
reaction compound atom-to-atom mapping libraries.


## Execution Environment

To use the [***crc_a2am***](/crc_a2am) package, please ensure that the
[***chytorch-rxnmap***](https://github.com/chython/chytorch-rxnmap),
[***epam.indigo***](https://github.com/epam/Indigo), and [***rxnmapper***](https://github.com/rxn4chemistry/rxnmapper)
libraries are available. The execution environment can be set up using [***conda***](https://docs.conda.io/en/latest)
and [***pip***](https://pip.pypa.io/en/stable) as follows:

```shell
conda env create -f environment.yaml

conda activate crc-a2am-env

pip install --no-build-isolation -e . --user
```


## Scripts
The [***scripts***](/scripts) directory is primarily meant to illustrate how to utilize the [***crc_a2am***](/crc_a2am)
package to run the chemical reaction atom-to-atom mapping on chemical reaction data.

```shell
# Example 1: Map a chemical reaction SMILES string.
python -m scripts.map_reaction_smiles --library "chytorch_rxnmap" --reaction_smiles "BrCCBr.COC(=O)c1cc(n[nH]1)C(F)(F)F>>COC(=O)c1cc(nn1CCBr)C(F)(F)F"
python -m scripts.map_reaction_smiles --library "epam_indigo" --reaction_smiles "CCN(CC)CCCl.[O-][N+](=O)c1ccc(S)cc1>>CCN(CC)CCSc1ccc(cc1)[N+]([O-])=O"
python -m scripts.map_reaction_smiles --library "rxnmapper" --reaction_smiles "ClC(=O)c1ccc(cc1)C#N.CCCCCCCc1ccccc1>>CCCCCCCc1ccc(cc1)C(=O)c1ccc(cc1)C#N"
```


## Supported Libraries

Currently, the [***crc_a2am***](/crc_a2am) package supports the following open-source chemical reaction compound
atom-to-atom mapping libraries:

1. The [**Chytorch RxnMap**](https://github.com/chython/chytorch-rxnmap) library utilizes a Transformer model adapted
   for processing chemical compound graphs. [**[1]**](#References)
2. The [**EPAM Indigo**](https://github.com/epam/Indigo) library utilizes a chemical compound graph-matching algorithm.
   [**[2]**](#References)
3. The [**RXNMapper**](https://github.com/rxn4chemistry/rxnmapper) library utilizes a chemically agnostic
   attention-guided Transformer model. [**[3]**](#References)


## License Information

This repository is published under the [**MIT**](/LICENSE) license. Please refer to the original publications for the
license information of individual open-source chemical reaction compound atom-to-atom mapping libraries.


## Contact

If you are interested in contributing to this repository by reporting bugs, submitting feedback or anything else that
might be beneficial, please feel free to do so via **GitHub Issues** or **e-mail**.


## References

1. Nugmanov, R., Dyubankova, N., Gedich, A., and Wegner, J.K. **Bidirectional Graphormer for Reactivity Understanding:
   Neural Network Trained to Reaction Atom-to-Atom Mapping Task**. *J. Chem. Inf. Model., 2022, 62, 14, 3307–3315*.
   DOI: https://doi.org/10.1021/acs.jcim.2c00344.
2. **EPAM Indigo**: https://lifescience.opensource.epam.com/indigo/index.html. Accessed on: March 6th, 2024.
3. Schwaller, P., Hoover, B., Reymond, J., Strobelt, H., and Laino, T. **Extraction of Organic Chemistry Grammar from
   Unsupervised Learning of Chemical Reactions**. *Sci. Adv., 2021, 7, 15*. DOI:
   https://doi.org/10.1126/sciadv.abe4166.
