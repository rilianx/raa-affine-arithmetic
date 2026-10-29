"""Global input box shared by the RAA forms (one interval [lb, ub] per variable)."""

_BOX = []


def set_box(box):
    """Set the input box: a sequence of (lb, ub) pairs, one per variable."""
    global _BOX
    new = [(float(a), float(b)) for a, b in box]
    for a, b in new:
        if a > b:
            raise ValueError(f"empty interval [{a}, {b}]")
    _BOX = new


def get_box():
    return _BOX
