""" The ``crc_a2am.chytorch_rxnmap`` package ``atom_to_atom_mapping`` module. """

from functools import partial
from logging import Logger
from typing import Collection, List, Optional, Tuple

from pqdm.processes import pqdm

from chython.files.daylight.smiles import smiles


class ChytorchRxnMapAtomToAtomMappingUtilities:
    """
    The `Chytorch RxnMap <https://pubs.acs.org/doi/10.1021/acs.jcim.2c00344>`_ library chemical reaction compound
    atom-to-atom mapping utilities class.
    """

    @staticmethod
    def map_reaction_smiles(
            reaction_smiles: str,
            custom_logger: Logger = None,
            **kwargs
    ) -> Tuple[Optional[str], Optional[float]]:
        """
        Map a chemical reaction SMILES string using the `Chytorch RxnMap <https://github.com/chython/chytorch-rxnmap>`_ library.

        :parameter reaction_smiles: The chemical reaction SMILES string.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions and methods: {
            ``chython.files.daylight.smiles.smiles``,
            ``chython.algorithms.mapping.attention.Attention.reset_mapping``
        }.
        :parameter custom_logger: The custom logger.

        :returns: The mapped chemical reaction SMILES string, and the chemical reaction compound atom-to-atom mapping score.
        """

        try:
            chytorch_rxnmap_reaction = smiles(
                reaction_smiles,
                ignore=kwargs.get("ignore", True),
                remap=kwargs.get("remap", False),
                ignore_stereo=kwargs.get("ignore_stereo", False),
                ignore_bad_isotopes=kwargs.get("ignore_bad_isotopes", False),
                keep_implicit=kwargs.get("keep_implicit", False),
                ignore_carbon_radicals=kwargs.get("ignore_carbon_radicals", False),
                ignore_aromatic_radicals=kwargs.get("ignore_aromatic_radicals", True)
            )

            chytorch_rxnmap_score = chytorch_rxnmap_reaction.reset_mapping(
                return_score=True,
                multiplier=kwargs.get("multiplier", 1.75),
                keep_reactants_numbering=kwargs.get("keep_reactants_numbering", False)
            )

            return format(chytorch_rxnmap_reaction, "m"), chytorch_rxnmap_score

        except Exception as exception_handle:
            if custom_logger is not None:
                custom_logger.exception(
                    msg=exception_handle
                )

            return None, None

    @staticmethod
    def map_reaction_smiles_strings(
            reaction_smiles_strings: Collection[str],
            **kwargs
    ) -> List[Tuple[Optional[str], Optional[float]]]:
        """
        Map the chemical reaction SMILES strings using the `Chytorch RxnMap <https://github.com/chython/chytorch-rxnmap>`_ library.

        :parameter reaction_smiles_strings: The chemical reaction SMILES strings.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying methods: {
            ``crc_a2am.chytorch_rxnmap.atom_to_atom_mapping.ChytorchRxnMapAtomToAtomMappingUtilities.map_reaction_smiles``
        }.

        :returns: The mapped chemical reaction SMILES strings, and the chemical reaction compound atom-to-atom mapping scores.
        """

        return pqdm(
            array=reaction_smiles_strings,
            function=partial(
                ChytorchRxnMapAtomToAtomMappingUtilities.map_reaction_smiles,
                **kwargs
            ),
            n_jobs=1,
            total=len(reaction_smiles_strings),
            ncols=150,
            ascii=True,
            desc="Mapping the chemical reaction SMILES strings (Library: Chytorch RxnMap | CPU: ???)"
        )
