#!/usr/bin/env python
import argparse
from ete3 import Tree
import pandas as pd
import re
from collections import Counter, defaultdict


def search_mat(
    newickfile, annotationfile, background, include_tips=False, skip_transitions=True
):
    """Function to search a mutation annotated tree from UShER and matUtils.

    Parameters
    ----------
    newickfile: str
        A path to a Newick format tree from [UShER](https://github.com/yatisht/usher)

    annotationfile: str
        A path to a tsv file from UShER's matUtils summary --translate function.

    substitution: str
        A background to annotate the mutations on in the format of an amino acid substitution, i.e. N501Y

    include_tips: bool
        Should mutations on tip branches be included in the analysis - default is False

    skip_transitions: bool
        Should mutations on branches with a substitution of interest be included on that background - default is False


    Returns
    -------
    pd.Dataframe
        Tidy dataframe with the columns `background`, `substitution`, and `count`.

    """

    #### START, PROCESS USER INPUTS ####
    print(f"\nStarting analysis for background: {background}\n")
    if include_tips:
        print("Including tip branches.\n")
    if skip_transitions:
        print(
            "Mutations on the branch of a new background are assigned to the parental background.\n"
        )

    # Regex to format amino acid substitutions
    aaregex = re.compile("([a-zA-Z*]+)([0-9]+)([a-zA-Z*]+)")

    # Format the background of interest
    search_substitution = aaregex.match(background).groups()

    # Load the tree from Newick file - takes the most time.
    print(f"Loading in the tree for {newickfile}\n")
    tree = Tree(newickfile, format=1)
    print("Finished loading tree.\n")

    print(f"Processing the annotation data for {annotationfile}\n")
    # Load and proces the annotations file
    annotation_df = pd.read_csv(annotationfile, sep="\t", low_memory=False)

    # parse the annotations in **Spike** into a dictionary object
    annotation_dict = {
        k: [sub.split(":")[1] for sub in v.split(";") if sub.split(":")[0] == "S"]
        for k, v in annotation_df.set_index("node_id")["aa_mutations"].to_dict().items()
    }

    # populate the dataframe with the missing nodes - presumably these have no mutations.
    existing_mutations = set(annotation_dict.keys())
    for node in [
        node.name for node in tree.traverse() if node.name not in existing_mutations
    ]:
        annotation_dict[node] = []

    # get the annotations and positions in dictionaries with the node names as keys
    annotations = {
        k: {
            aaregex.match(sub).groups()
            for sub in v
            if aaregex.match(sub).group(1) != aaregex.match(sub).group(3)
        }  # TODO: Remove "Substitutions" that are actually synonymous
        if v
        else []
        for k, v in annotation_dict.items()
    }  # full mutations
    positions = {
        k: {
            aaregex.match(sub).group(2)
            for sub in v
            if aaregex.match(sub).group(1) != aaregex.match(sub).group(3)
        }  # TODO: Remove "Substitutions" that are actually synonymous
        if v
        else []
        for k, v in annotation_dict.items()
    }  # just positions

    print("Finished processing annotation data.\n")
    #### END, PROCESS USER INPUTS ####

    #### START, TREE TRAVERSAL AND ANNOTATION ####
    print("Starting the tree traversal...\n")
    # Dictionary to hold the state of each node in the tree
    background_state = dict()

    # Dictionary to hold the mutations associated with each background
    substitution_dict = defaultdict(list)

    # Populate the background dictionary with the background of the root
    # TODO: This assumes that the root state is the first amino acid in the search sub. Not always true.
    background_state[tree.name] = (
        search_substitution[0],
        search_substitution[1],
        search_substitution[0],
    )

    # Traverse the tree, skipping the root (already annotated)
    for node in tree.iter_descendants("preorder"):

        # Skip tips (leaf nodes) if include_tips is True
        if node.is_leaf() and not include_tips:
            continue

        # There are no new Spike substitutions on the branch preceding this node.
        if not annotations[node.name]:
            background_state[node.name] = background_state[
                node.up.name
            ]  # inherit state from parent

        # None of the Spike substitutions are at the search position
        elif search_substitution[1] not in positions[node.name]:
            background_state[node.name] = background_state[
                node.up.name
            ]  # inherit state from parent

            # add these spike mutations to the background that's been assigned
            substitution_dict[background_state[node.name]].extend(
                annotations[node.name]
            )

        # The position of interest is present on the branch preceding this node
        else:
            # TODO: assumes that there is only one position of interest on the branch, handle errors where annotation is wrong
            for substitution in annotations[node.name]:
                if substitution[1] == search_substitution[1]:
                    background_state[
                        node.name
                    ] = substitution  # change state of this node

            # Mutation on these branches are considered to belong to the background of parent (not the new background)
            if skip_transitions:
                substitution_dict[background_state[node.up.name]].extend(
                    annotations[node.name]
                )
            # Mutations on these branches are appended to current background
            else:
                substitution_dict[background_state[node.name]].extend(
                    annotations[node.name]
                )

    print("Finished the tree traversal.\n")

    #### END, TREE TRAVERSAL AND ANNOTATION ####

    #### START, PROCESS SUBSTITUTIONS AND RETURN ####

    # List to store substitutions
    substitution_count_list = list()

    # Convert the dictionary to a list of lists for pandas
    for background in substitution_dict.keys():
        # Count the substitutions in a given background
        substitution_count = Counter(substitution_dict[background])

        # List of columns as tuples
        substitution_count_list.extend(
            [
                ("".join(background), "".join(k), k[0], k[1], k[2], v)
                for k, v in substitution_count.items()
            ]
        )

    # Return a dataframe of counts of substitutions on each background
    return pd.DataFrame(
        substitution_count_list,
        columns=["Background", "Substitution", "WT", "POS", "MUT", "Count"],
    )

    #### END, PROCESS SUBSTITUTIONS AND RETURN ####


if __name__ == "__main__":

    # Command line interface
    parser = argparse.ArgumentParser(
        description="Search a mutation annotated tree for substitutions on a background of interest."
    )

    parser.add_argument(
        "--tree",
        "-t",
        type=str,
        required=True,
        help="Path to a newick format tree generated by UShER.",
    )
    parser.add_argument(
        "--annotations",
        "-a",
        type=str,
        required=True,
        help="Path to a *.tsv file from `matUtils summary --translate`.",
    )
    parser.add_argument(
        "--substitution",
        "-s",
        type=str,
        required=True,
        help="Background substitution to search the tree with, i.e. N501Y",
    )
    parser.add_argument(
        "--outpath",
        "-o",
        type=str,
        required=True,
        help="Desired path for output *.csv file containing counts of substitutions on each background.",
    )
    parser.add_argument(
        "--include_tips",
        type=bool,
        required=False,
        default=False,
        help="Should mutations on tips be included in the output? Default: False",
    )
    parser.add_argument(
        "--skip_transitions",
        type=bool,
        required=False,
        default=True,
        help="Should substitutions on branches leading up to a change in background be included on the previous background? Default: True",
    )

    args = parser.parse_args()

    # Main function
    substitution_df = search_mat(
        args.tree,
        args.annotations,
        args.substitution,
        args.include_tips,
        args.skip_transitions,
    )

    substitution_df.to_csv(args.outpath, index=False)

    print(f"Done! A table of substitution counts is saved at {args.outpath}.")