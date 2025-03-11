
def printStreetContent(street, processBundle=False):
    if street.bundle and not processBundle:
        return
    
    print(
        street.debugInfo() +\
        " - ".join( item.debugInfo() for item in street.iterItems() )
    )


def printBundleContent(bundle):
    print("Bundle:{0}".format(bundle.id))
    for street in bundle.streetsHead:
        printStreetContent(street, processBundle=True)
    print("\n")