import math

def get_bin(x, bin_size, offset=1088):

    '''
    Get the bin number based on the base pair number in a sequence (ex: 2500th bp is in the 2nd bin for bins of size 2048).
    Note: x is the distance to the start of the sequence, not the distance to mat_start !!!

    '''
    x_bin = math.ceil(x/bin_size) - offset

    return x_bin
