
def getStreetDebugInfo(street):
    return street.debugInfo() + " - ".join( item.debugInfo() for item in street.iterItems() )


def printStreetContent(street, includeEnds):
    if includeEnds:
        pred = "None" if street.pred is None else "I({})".format(street.pred.intersection.id)
        succ = "None" if street.succ is None else "I({})".format(street.succ.intersection.id)
        print("{0}{1}{2} -> {3} -> {4}{5}".format(
            street.debugInfo(),
            pred,
            '' if street.pred is None else "[{0:.2f},{1:.2f}]".format(street.pred.intersection.location[0], street.pred.intersection.location[1]),
            " - ".join( item.debugInfo() for item in street.iterItems() ),
            succ,
            '' if street.succ is None else "[{0:.2f},{1:.2f}]".format(street.succ.intersection.location[0], street.succ.intersection.location[1])
        ))
    else:
        print( getStreetDebugInfo(street) )


def printBundleContent(bundle):
    # check if we got a simple case: <bundle.streetsHead> and <bundle.streetsTail> have the same Streets
    streetsHead, streetsTail = bundle.streetsHead, bundle.streetsTail
    
    print("###########")
    print("Bundle:{}".format(bundle.id))
    print("##########")
    
    if len(streetsHead) == len(streetsTail) == sum(1 for streetHead, streetTail in zip(streetsHead, streetsTail) if streetHead is streetTail):
        # a simple case: <bundle.streetsHead> and <bundle.streetsTail> have the same Streets
        for street in streetsHead:
            print( getStreetDebugInfo(street) )
    else:
        streets = dict()
        intersections = dict()
        unprocessedIntersections = set()
        
        for street in streetsHead:
            if (street.pred is None or street.pred.item is bundle) and (street.succ is None or street.succ.item is bundle):
                # still a simple case: <street> spans the whole <bundle>
                printStreetContent(street, includeEnds=False)
            else:
                streets[street.id] = street
                leaving = street.pred is None or street.pred.item is bundle
                intersection = (street.succ if leaving else street.pred).intersection
                if not intersection.id in unprocessedIntersections:
                    intersections[intersection.id] = intersection
                    unprocessedIntersections.add(intersection.id)
                if leaving:
                    print("{0} -> I({1})".format( getStreetDebugInfo(street), intersection.id ))
        
        for street in streetsTail:
            if not ( (street.pred is None or street.pred.item is bundle) and (street.succ is None or street.succ.item is bundle) ):
                streets[street.id] = street
                leaving = street.pred is None or street.pred.item is bundle
                intersection = (street.succ if leaving else street.pred).intersection
                if not intersection.id in intersections:
                    intersections[intersection.id] = intersection
                    unprocessedIntersections.add(intersection.id)
                if leaving:
                    print("{0} -> I({1})".format( getStreetDebugInfo(street), intersection.id ))
        
        while unprocessedIntersections:
            _id = next(iter(unprocessedIntersections))
            processBundleIntersection(intersections[_id], intersections, unprocessedIntersections, streets)
            unprocessedIntersections.remove(_id)
        
        for _id in intersections:
            intersection = intersections[_id]
            print("I({0}) [{1:.2f},{2:.2f}]:".format(_id, intersection.location[0], intersection.location[1]))
            for connector in intersection.iterConnectors():
                if connector.leaving:
                    street = connector.item
                    endStr = "None" if street.succ is None else\
                        ("B({})".format(street.succ.item.id) if street.succ.item is bundle else "I({})".format(street.succ.intersection.id))
                    print("  {0} -> {1}".format( getStreetDebugInfo(street), endStr) )
    
    print("\n")


def processBundleIntersection(intersection, intersections, unprocessedIntersections, streets):
    connector = intersection.startConnector
    while True:
        street = connector.item
        if not street.id in streets:
            streets[street.id] = street
            oppositeConnector = street.pred if street.succ is connector  else street.succ
            if oppositeConnector and not oppositeConnector.intersection.id in intersections:
                intersections[oppositeConnector.intersection.id] = oppositeConnector.intersection
                unprocessedIntersections.add(oppositeConnector.intersection.id)
        connector = connector.succ
        if connector is intersection.startConnector:
            break