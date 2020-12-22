
def limit_title(title, max_title_length):
    """
    Limit the length of each line of a tile to the given number of characters,
    adding an ellipsis if a given line is too long.

    Parameters
    ----------
    title : str
        The title to limit

    max_title_length : int
        The maximum length of each line of the title

    Returns
    -------
    title : str
        The title after limiting each line's length
    """
    parts = list()
    for part in title.split('\n'):
        if len(part) > max_title_length:
            part = '{}...'.format(part[0:max_title_length])
        parts.append(part)
    title = '\n'.join(parts)
    return title
