

def row_to_block(row):
    if row < 13:
        return 1
    elif row < 25:
        return 2
    elif row < 37:
        return 3
    else:
        raise Exception('Too many rows in array. Add Blocks to Django Model (Gal File) and Use them!')
