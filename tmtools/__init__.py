from ._bindings import tm_align as tm_align_function  # noqa

SEQ_LENGTH_MIN = 3

def tm_align(*args, **kwargs):
    """
    Wrapper function for `tm_align` for error handling. 
    """
    # Check seq1 and seq2
    coords1, coords2, seq1, seq2 = args 
    if len(seq1) < SEQ_LENGTH_MIN or len(seq2) < SEQ_LENGTH_MIN:
        raise ValueError(f"Sequences must be at least {SEQ_LENGTH_MIN} characters long.")
    try:
        return tm_align_function(*args, **kwargs)
    except Exception as e:
        print(f"Error: {e}")
        
    
