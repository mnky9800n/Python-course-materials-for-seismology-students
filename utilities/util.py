def replace_unique_items(iterable, replace_with=None):
    """
    replaces items after the first unique item in a list
    with replace_with.
    
    ::default behavior::
    ----------------------------------
    In  : x = [1, 1, 2, 2]
    Out : replace_unique_items(x)
    [1, None, 2, None]
    """
    result = []
    for item in iterable:
        if item not in result:
            result.append(item)
        else:
            result.append(replace_with)
    return result
