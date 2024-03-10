""" The ``crc_a2am.epam_indigo`` package ``atom_to_atom_mapping`` module. """

from functools import partial
from indigo import Indigo
from logging import Logger
from typing import Collection, List, Optional, Tuple

from pqdm.processes import pqdm


class EpamIndigoAtomToAtomMappingUtilities:
    """
    The `EPAM Indigo <https://lifescience.opensource.epam.com/indigo>`_ library chemical reaction compound
    atom-to-atom mapping utilities class.
    """

    @staticmethod
    def map_reaction_smiles(
            reaction_smiles: str,
            timeout_period_ms: int = 10000,
            handle_existing_atom_to_atom_mapping: str = "discard",
            ignore_charges: bool = False,
            ignore_isotopes: bool = False,
            ignore_valences: bool = False,
            ignore_radicals: bool = False,
            canonicalize_reaction_smiles: bool = True,
            custom_logger: Logger = None
    ) -> Tuple[Optional[str], Optional[bool]]:
        """
        Map a chemical reaction SMILES string using the
        `EPAM Indigo <https://lifescience.opensource.epam.com/indigo/api/index.html#reaction-atom-to-atom-mapping>`_ library.

        :parameter reaction_smiles: The chemical reaction SMILES string.
        :parameter timeout_period_ms: The maximum amount of time in milliseconds that may be spent on the chemical
            reaction compound atom-to-atom mapping.
        :parameter handle_existing_atom_to_atom_mapping: The indicator of how the existing chemical reaction compound
            atom-to-atom mapping should be handled. The exclusive value choices are: {
                ``"alter"``, ``"clear"``, ``"discard"``, ``"keep"``
            }.
        :parameter ignore_charges: The indicator of whether the chemical reaction compound atom charges should be ignored.
        :parameter ignore_isotopes: The indicator of whether the chemical reaction compound atom isotopes should be ignored.
        :parameter ignore_valences: The indicator of whether the chemical reaction compound atom valences should be ignored.
        :parameter ignore_radicals: The indicator of whether the chemical reaction compound atom radicals should be ignored.
        :parameter canonicalize_reaction_smiles: The indicator of whether the chemical reaction SMILES string should be canonicalized.
        :parameter custom_logger: The custom logger.

        :returns: The mapped chemical reaction SMILES string, and the indicator of whether the chemical reaction
            compound atom-to-atom mapping was completed without errors.
        """

        try:
            epam_indigo = Indigo()

            epam_indigo.setOption(
                option="aam-timeout",
                value1=timeout_period_ms
            )

            epam_indigo_reaction = epam_indigo.loadReactionSmarts(
                string=reaction_smiles
            )

            epam_indigo_status_code = epam_indigo_reaction.automap(
                mode="".join([
                    handle_existing_atom_to_atom_mapping
                    if handle_existing_atom_to_atom_mapping in ["alter", "clear", "discard", "keep", ] else "discard",
                    " ignore_charges" if ignore_charges else "",
                    " ignore_isotopes" if ignore_isotopes else "",
                    " ignore_valence" if ignore_valences else "",
                    " ignore_radicals" if ignore_radicals else ""
                ])
            )

            if canonicalize_reaction_smiles:
                return epam_indigo_reaction.canonicalSmiles(), True if epam_indigo_status_code == 1 else False

            else:
                return epam_indigo_reaction.smiles(), True if epam_indigo_status_code == 1 else False

        except Exception as exception_handle:
            if custom_logger is not None:
                custom_logger.exception(
                    msg=exception_handle
                )

            return None, None

    @staticmethod
    def map_reaction_smiles_strings(
            reaction_smiles_strings: Collection[str],
            number_of_cpu_cores: int = 1,
            **kwargs
    ) -> List[Tuple[Optional[str], Optional[bool]]]:
        """
        Map the chemical reaction SMILES strings using the
        `EPAM Indigo <https://lifescience.opensource.epam.com/indigo/api/index.html#reaction-atom-to-atom-mapping>`_ library.

        :parameter reaction_smiles_strings: The chemical reaction SMILES strings.
        :parameter number_of_cpu_cores: The number of CPU cores that should be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying methods: {
            ``crc_a2am.epam_indigo.atom_to_atom_mapping.EpamIndigoAtomToAtomMappingUtilities.map_reaction_smiles``
        }.

        :returns: The mapped chemical reaction SMILES strings, and the indicators of whether the chemical reaction
            compound atom-to-atom mappings were completed without errors.
        """

        return pqdm(
            array=reaction_smiles_strings,
            function=partial(
                EpamIndigoAtomToAtomMappingUtilities.map_reaction_smiles,
                **kwargs
            ),
            n_jobs=number_of_cpu_cores,
            total=len(reaction_smiles_strings),
            ncols=150,
            ascii=True,
            desc="Mapping the chemical reaction SMILES strings (Library: EPAM Indigo | CPU: {0:d})".format(
                number_of_cpu_cores
            )
        )
