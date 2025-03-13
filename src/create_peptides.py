import pandas as pd
import os 
import argparse
import dask.bag as db
from dask.bag import Bag
from dask.diagnostics import ProgressBar
from typing import Union, List, Tuple

def format_time(duration_seconds: Union[int, float]) -> str:
    """
    Convert a seconds to a readable string in hours, minutes, and seconds.
    
    Args:
        duration_seconds (Union[int, float]): The time duration in seconds.
    
    Returns:
        str: A formatted string
    """
    hours, remainder = divmod(duration_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds"


def init_bag(building_blocks: Union[List[str], str]) -> Bag:
    """
    Create a Dask Bag from a list of building blocks or from a CSV file containing them.
    
    If a string is provided, it is treated as the file path to a CSV file. The CSV file must have a column
    named 'SMILES'. If a list is provided, it is directly used as the source for the bag.
    
    Args:
        building_blocks (Union[List[str], str]): Either a list of SMILES strings or a file path (str) to a CSV file.
    
    Returns:
        Bag: A Dask Bag containing the building block SMILES strings.
    """
    if isinstance(building_blocks, str):
        # Read the CSV file using the provided file path.
        df = pd.read_csv(building_blocks)  # Expected CSV format: columns include 'SMILES'
        smiles_list = df["SMILES"].tolist()  # Call the method to obtain the list
        bag = db.from_sequence(smiles_list, npartitions=16)
    elif isinstance(building_blocks, list):
        bag = db.from_sequence(building_blocks, npartitions=16)
    else:
        raise ValueError("building_blocks must be either a file path (str) or a list of strings.")
    return bag


def cartesian_combination(bag: Bag) -> Bag:
    """
    Generate all possible 5-length Cartesian combinations of building blocks from the input bag.
    
    The combination is built in three steps:
      1. Generate 2-length combinations (pairs) using the Cartesian product of the bag with itself.
      2. Form 4-length combinations by taking the Cartesian product of the pairs with themselves.
      3. Append one additional building block to each 4-length combination to form a 5-length combination.
    
    Args:
        bag (Bag): A Dask Bag containing building block SMILES strings.
    
    Returns:
        Bag: A Dask Bag where each element is a tuple of 5 SMILES strings.
    """
    # Create all possible pairs (2-length combinations)
    pairs = bag.product(bag)  # Each element is a tuple: (block_i, block_j)

    # Create 4-length combinations by taking the product (cartesian) of (x,y)(y,z)
    quadruple_nested = pairs.product(pairs)  # Each element is ((x1, x2), (x3, x4))
    def flatten_quadruple(nested: Tuple[Tuple[str, str], Tuple[str, str]]) -> Tuple[str, str, str, str]:
        (x1, x2), (x3, x4) = nested
        return (x1, x2, x3, x4)
    quadruples = quadruple_nested.map(flatten_quadruple)
    
    #Create 5-length combinations by appending one more building block to each 4-length tuple
    quintuple_nested = quadruples.product(bag)  # Each element is ((x1, x2, x3, x4), x5)
    def flatten_quintuple(nested: Tuple[Tuple[str, str, str, str], str]) -> Tuple[str, str, str, str, str]:
        quadruple, x5 = nested
        return (*quadruple, x5)
    quintuples = quintuple_nested.map(flatten_quintuple)
    
    return quintuples


def combine_smiles(building_blocks: Union[List[str], Tuple[str, ...]]) -> str:
    """
    Combine a sequence of building block SMILES strings into a single peptide SMILES string.
    
    The combination is achieved by concatenating the SMILES strings and appending an 'O' at the end.
    
    Args:
        building_blocks (Union[List[str], Tuple[str, ...]]): A sequence (e.g., tuple) of building block SMILES strings.
    
    Returns:
        str: The combined peptide SMILES string.
    """
    return "".join(building_blocks) + 'O'


def save_bag(bag: Bag, npartitions: int, output_path: str) -> None:
    """
    Save the contents of a Dask Bag to text files, with one file per partition.
    
    The bag is first repartitioned to the specified number of partitions. Each partition is then written
    to a separate text file with a name based on the given output prefix.
    
    Args:
        bag (Bag): A Dask Bag containing peptide sequences.
        npartitions (int): The number of partitions; one file will be generated per partition.
        output_prefix (str): The output file path prefix. Files will be saved as '<output_prefix>-part-*.txt'.
    """
    import os
    # Ensure the directory for the output prefix exists
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    
    bag = bag.repartition(npartitions=npartitions)
    with ProgressBar():
        bag.to_textfiles(f"{output_path}-part-*.txt")
    print("Exporting completed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Peptide forge tool for generating peptide SMILES strings from building blocks."
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to a CSV file with building blocks. The CSV file must include a column named 'SMILES'."
    )
    parser.add_argument(
        "--building_blocks",
        type=str,
        nargs="*",
        help="List of building block SMILES strings to use if no input CSV file is provided."
    )
    parser.add_argument(
        "--npartitions",
        type=int,
        default=250,
        help="Number of partitions to use for saving the output. One file will be generated per partition."
    )
    parser.add_argument(
        "--output_path",
        type=str,
        help="Output file path. Files will be saved as '<output_path>-part-*.txt'. The number of parts is the same as"
        "the number of partitions"
    )
    args = parser.parse_args()

    import time
    start_time = time.time()

    # If an input CSV file is provided, use that; otherwise, use the provided list.
    if args.input:
        building_blocks_bag = init_bag(args.input)
    elif args.building_blocks:
        building_blocks_bag = init_bag(args.building_blocks)

    else:
        raise ValueError("An input csv path for SMILES or building block as argumnets must be passed")

    # Generate 5-length combinations of building blocks.
    combined_combinations = cartesian_combination(bag=building_blocks_bag)
    
    # Convert each 5-length tuple of building blocks into a combined peptide SMILES string.
    combined_smiles_bag = combined_combinations.map(combine_smiles)
    
    # Save the resulting bag to text files.
    save_bag(combined_smiles_bag, npartitions=args.npartitions, output_path=args.output_path)
    
    elapsed_time = time.time() - start_time
    print(f"Peptide forge took: {format_time(elapsed_time)} to complete")