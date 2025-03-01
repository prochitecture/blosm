
def printStreetContent(street, processBundle=False):
    if street.bundle and not processBundle:
        return
    
    print(
        "{0}:".format(street.id) +\
        " - ".join(
            (
                "(Int)" if item.__class__.__name__ == "Intersection" else\
                "{0}({1})({2})".format(item.__class__.__name__, item.id, "aeroway" if "aeroway" in item.tags else "highway")
            )\
            for item in street.iterItems()
        )
    )


def printBundleContent(bundle):
    print("Bundle:{0}".format(bundle.id))
    for street in bundle.streetsHead:
        printStreetContent(street, processBundle=True)
    print("\n")