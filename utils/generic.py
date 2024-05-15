def flatten_list_to_set(li):
    flat_set = set()
    for xs in li:
        for x in xs:
            flat_set.add(x)
    return flat_set


def flatten_list_to_list(li):
    flat_list = list()
    for xs in li:
        for x in xs:
            flat_list.append(x)
    return flat_list
