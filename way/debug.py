
def printStreetContent(street, processBundle=False):
    if street.bundle and not processBundle:
        return
    
    for item in street.iterItems():
        if item.__class__.__name__ == "Intersection":
            t = 0
    
    print(
        "{0}:".format(street.id) +\
        " - ".join(
            (
                "({0}/{1})".format(
                    item.leftHead.item.head.tags["highway"] if item.leftHead else '-',
                    item.rightHead.item.head.tags["highway"] if item.rightHead else '-'
                )\
                if item.__class__.__name__ == "Intersection" else "{0}({1})".format(item.__class__.__name__, item.id)
            )\
            for item in street.iterItems()
        )
    )


def printBundleConten(bundle):
    print("Bundle:{0}".format(bundle.id))
    for street in bundle.streetsHead:
        printStreetContent(street, processBundle=True)
    print("\n")