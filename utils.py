import numpy as np
import os
from itertools import chain, combinations

def row_to_block(row):
    if row < 13:
        return 1
    elif row < 25:
        return 2
    elif row < 37:
        return 3
    else:
        raise Exception('Too many rows in array. Add Blocks to Django Model (Gal File) and Use them!')

def all_same(items):
    return all(x == items[0] for x in items)

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def all_subsets(ss):
  return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))