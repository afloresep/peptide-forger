import sys
import os
from typing import List

# Append the parent directory to sys.path to import modules from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.create_peptides import * 

def test_cartesian():
    """
    Test the cartesian products of a single bag returns the 
    correct amount of combinations
    """
    bb_bag_size2 = init_bag(["BB1", "BB2"])
    bb_bag_size4 = init_bag(["BB1", "BB2", "BB3", "BB4"])

    #Do cartesian combination 
    bb_combined2 = cartesian_combination(bag=bb_bag_size2)
    bb_combined4 = cartesian_combination(bag=bb_bag_size4)

    bb_combined2 = bb_combined2.map(combine_smiles)
    bb_combined4 = bb_combined4.map(combine_smiles)

    assert bb_combined4.count().compute() == (4**5), "Cartesian product returned wrong number of combinations"
    assert bb_combined2.count().compute() == (2**5), "Cartesian product returned wrong number of combinations"


