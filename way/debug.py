
def printStreetContent(street, processBundle=False):
    if street.bundle and not processBundle:
        return
    
    print(
        "{0}:".format(street.id) +\
        " - ".join(
            (
                "({0}({1})/{2}({3})".format(
                    (item.leftHead.item.head.tags["highway"] if item.leftHead.leaving else item.leftHead.item.tail.tags["highway"]) if item.leftHead else '-',
                    (item.leftHead.item.head.id if item.leftHead.leaving else item.leftHead.item.tail.id) if item.leftHead else '',
                    (item.rightHead.item.head.tags["highway"] if item.rightHead.leaving else item.rightHead.item.tail.tags["highway"]) if item.rightHead else '-',
                    (item.rightHead.item.head.id if item.rightHead.leaving else item.rightHead.item.tail.id) if item.rightHead else ''
                )\
                if item.__class__.__name__ == "Intersection" else "{0}({1})".format(item.__class__.__name__, item.id)
            )\
            for item in street.iterItems()
        )
    )


def printBundleContent(bundle):
    print("Bundle:{0}".format(bundle.id))
    for street in bundle.streetsHead:
        printStreetContent(street, processBundle=True)
    print("\n")