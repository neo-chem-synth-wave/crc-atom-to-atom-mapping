""" The ``scripts`` directory ``map_reaction_smiles`` script. """

from argparse import ArgumentParser, Namespace
from warnings import filterwarnings

filterwarnings(
    action="ignore"
)

from crc_a2am.chytorch_rxnmap import ChytorchRxnMapAtomToAtomMappingUtilities
from crc_a2am.epam_indigo import EpamIndigoAtomToAtomMappingUtilities
from crc_a2am.rxnmapper import RXNMapperAtomToAtomMappingUtilities


def get_script_arguments(
) -> Namespace:
    """
    Get the script arguments.

    :returns: The script arguments.
    """

    argument_parser = ArgumentParser()

    argument_parser.add_argument(
        "-l",
        "--library",
        type=str,
        choices=[
            "chytorch_rxnmap",
            "epam_indigo",
            "rxnmapper"
        ],
        required=True,
        help="The indicator of which chemical reaction compound atom-to-atom mapping library should be utilized."
    )

    argument_parser.add_argument(
        "-s",
        "--reaction_smiles",
        type=str,
        required=True,
        help="The chemical reaction SMILES string."
    )

    return argument_parser.parse_args()


if __name__ == "__main__":
    script_arguments = get_script_arguments()

    if script_arguments.library == "chytorch_rxnmap":
        chytorch_rxnmap_mapped_reaction_smiles, chytorch_rxnmap_score = \
            ChytorchRxnMapAtomToAtomMappingUtilities.map_reaction_smiles(
                reaction_smiles=script_arguments.reaction_smiles
            )

        print("Mapped Reaction SMILES: '{0:s}'".format(chytorch_rxnmap_mapped_reaction_smiles))
        print("Score (Chytorch RxnMap): {0:.5f}".format(chytorch_rxnmap_score))

    elif script_arguments.library == "epam_indigo":
        epam_indigo_mapped_reaction_smiles, epam_indigo_indicator = \
            EpamIndigoAtomToAtomMappingUtilities.map_reaction_smiles(
                reaction_smiles=script_arguments.reaction_smiles
            )

        print("Mapped Reaction SMILES: '{0:s}'".format(epam_indigo_mapped_reaction_smiles))
        print("Completed without Errors (EPAM Indigo): {0:}".format(epam_indigo_indicator))

    elif script_arguments.library == "rxnmapper":
        rxnmapper_mapped_reaction_smiles, rxnmapper_score = \
            RXNMapperAtomToAtomMappingUtilities.map_reaction_smiles(
                reaction_smiles=script_arguments.reaction_smiles
            )

        print("Mapped Reaction SMILES: '{0:s}'".format(rxnmapper_mapped_reaction_smiles))
        print("Score (RXNMapper): {0:.5f}".format(rxnmapper_score))
