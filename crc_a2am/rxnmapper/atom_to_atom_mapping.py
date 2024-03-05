""" The ``crc_a2am.rxnmapper`` package ``atom_to_atom_mapping`` module. """

from logging import Logger
from math import ceil
from typing import Collection, List, Optional, Tuple

from rxnmapper.core import RXNMapper
from tqdm.auto import tqdm


class RXNMapperAtomToAtomMappingUtilities:
    """
    The `RXNMapper <https://www.science.org/doi/10.1126/sciadv.abe4166>`_ library chemical reaction compound
    atom-to-atom mapping utilities class.
    """

    @staticmethod
    def map_reaction_smiles(
            reaction_smiles: str,
            custom_logger: Logger = None,
            **kwargs
    ) -> Tuple[Optional[str], Optional[float]]:
        """
        Map a chemical reaction SMILES string using the `RXNMapper <https://github.com/rxn4chemistry/rxnmapper>`_
        library.

        :parameter reaction_smiles: The chemical reaction SMILES string.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying methods: {
            ``rxnmapper.core.RXNMapper.get_attention_guided_atom_maps`` }.
        :parameter custom_logger: The custom logger.

        :returns: The mapped chemical reaction SMILES string, and the chemical reaction compound atom-to-atom mapping
            score.
        """

        try:
            rxnmapper_model = RXNMapper()

            rxnmapper_model_output = rxnmapper_model.get_attention_guided_atom_maps(
                rxns=[reaction_smiles],
                zero_set_p=kwargs.get("zero_set_p", True),
                zero_set_r=kwargs.get("zero_set_r", True),
                canonicalize_rxns=kwargs.get("canonicalize_rxns", True),
                detailed_output=kwargs.get("detailed_output", False),
                absolute_product_inds=kwargs.get("absolute_product_inds", False),
                force_layer=kwargs.get("force_layer", None),
                force_head=kwargs.get("force_head", None)
            )

            return rxnmapper_model_output[0].get("mapped_rxn", None), rxnmapper_model_output[0].get("confidence", None)

        except Exception as exception_handle:
            if custom_logger is not None:
                custom_logger.exception(
                    msg=exception_handle
                )

            return None, None

    @staticmethod
    def map_reaction_smiles_strings(
            reaction_smiles_strings: Collection[str],
            batch_size: int = 10,
            custom_logger: Logger = None,
            **kwargs
    ) -> Optional[List[Tuple[Optional[str], Optional[float]]]]:
        """
        Map the chemical reaction SMILES strings using the `RXNMapper <https://github.com/rxn4chemistry/rxnmapper>`_
        library.

        :parameter reaction_smiles_strings: The chemical reaction SMILES strings.
        :parameter batch_size: The batch size.
        :parameter custom_logger: The custom logger.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying methods: {
            ``rxnmapper.core.RXNMapper.get_attention_guided_atom_maps`` }.

        :returns: The mapped chemical reaction SMILES strings, and the chemical reaction compound atom-to-atom mapping
            scores.
        """

        try:
            rxnmapper_model = RXNMapper()

            mapped_reaction_smiles_strings_and_scores = list()

            for reaction_smiles_index in tqdm(
                iterable=range(0, len(reaction_smiles_strings), batch_size),
                total=ceil(len(reaction_smiles_strings) / batch_size),
                ncols=150,
                ascii=True,
                desc="Map the chemical reaction SMILES strings (Library: RXNMapper | Batch Size: {0:d} | GPU: ???)".format(
                    batch_size
                )
            ):
                rxnmapper_model_outputs = rxnmapper_model.get_attention_guided_atom_maps(
                    rxns=list(reaction_smiles_strings[reaction_smiles_index: reaction_smiles_index + batch_size]),
                    zero_set_p=kwargs.get("zero_set_p", True),
                    zero_set_r=kwargs.get("zero_set_r", True),
                    canonicalize_rxns=kwargs.get("canonicalize_rxns", True),
                    detailed_output=kwargs.get("detailed_output", False),
                    absolute_product_inds=kwargs.get("absolute_product_inds", False),
                    force_layer=kwargs.get("force_layer", None),
                    force_head=kwargs.get("force_head", None)
                )

                for rxnmapper_model_output in rxnmapper_model_outputs:
                    mapped_reaction_smiles_strings_and_scores.append((
                        rxnmapper_model_output.get("mapped_rxn", None),
                        rxnmapper_model_output.get("confidence", None),
                    ))

            return mapped_reaction_smiles_strings_and_scores

        except Exception as exception_handle:
            if custom_logger is not None:
                custom_logger.exception(
                    msg=exception_handle
                )

            return None
