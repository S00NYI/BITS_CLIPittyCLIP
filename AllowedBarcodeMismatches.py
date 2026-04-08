import argparse
import itertools
import sys

def barcodes_N(barcode, n):
    positions = range(len(barcode))
    for pos in itertools.combinations(positions, n):
        new_barcode = list(barcode)
        for p in pos:
            new_barcode[p] = 'N'
        yield ''.join(new_barcode)

def main():
    parser = argparse.ArgumentParser(description='Generate barcodes with N substitutions')
    parser.add_argument('-barcodes', required=True, help='List of input barcodes separated by spaces (e.g. "ATCACG CGATGT TTAGGC")')
    parser.add_argument('-mismatches', type=int, required=True, help='Number of N substitutions allowed')

    args = parser.parse_args()

    if not args.barcodes and not args.mismatches:
        print("Use -h to see options.")
        sys.exit(0)

    barcodes = args.barcodes.split()
    Ns_allowed = args.mismatches

    new_barcodes = []

    for barcode in barcodes:
        for combination in barcodes_N(barcode, Ns_allowed):
            new_barcodes.append(combination)

    # Check whether all the barcodes with Ns are unique or not.
    # If True, there is no overlap.
    # If False, there is overlap (we lose uniqueness).
    result = len(new_barcodes) == len(list(set(new_barcodes)))

    if result == True:
        print(f"Uniqueness of barcodes are maintained with {Ns_allowed} N substitutions.")
    else:
        print(f"Uniqueness of barcodes are not maintained with {Ns_allowed} N substitutions.")

if __name__ == "__main__":
    main()
