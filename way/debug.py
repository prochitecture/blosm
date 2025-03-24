
def getStreetDebugInfo(street):
    return street.debugInfo() + " - ".join( item.debugInfo() for item in street.iterItems() )


def printStreetContent(street):
    print( getStreetDebugInfo(street) )


def printBundleContent(bundle):
    # check if we got a simple case: <bundle.streetsHead> and <bundle.streetsTail> have the same Streets
    streetsHead, streetsTail = bundle.streetsHead, bundle.streetsTail
    
    if len(streetsHead) == len(streetsTail) == sum(1 for streetHead, streetTail in zip(streetsHead, streetsTail) if streetHead is streetTail):
        # a simple case: <bundle.streetsHead> and <bundle.streetsTail> have the same Streets
        print("Bundle:{}".format(bundle.id))
        for street in streetsHead:
            print( getStreetDebugInfo(street) )
    else:
        streets = dict()
        intersections = dict()
        unprocessedIntersections = set()
        
        for street in streetsHead:
            if (street.pred is None or street.pred.item is bundle) and (street.succ is None or street.succ.item is bundle):
                # still a simple case: <street> spans the whole <bundle>
                printStreetContent(street)
            else:
                streets[street.id] = street
                #unvisitedStreets[street.id] = street
                intersection = (street.succ if street.pred.item is bundle else street.pred).intersection
                if not intersection.id in unprocessedIntersections:
                    intersections[intersection.id] = intersection
                    unprocessedIntersections.add(intersection.id)
        
        for street in streetsTail:
            if not ( (street.pred is None or street.pred.item is bundle) and (street.succ is None or street.succ.item is bundle) ):
                streets[street.id] = street
                #unvisitedStreets[street.id] = street
                intersection = (street.succ if (street.pred is None or street.pred.item is bundle) else street.pred).intersection
                if not intersection.id in intersections:
                    intersections[intersection.id] = intersection
                    unprocessedIntersections.add(intersection.id)
        
        while unprocessedIntersections:
            _id = next(iter(unprocessedIntersections))
            processBundleIntersection(intersections[_id], intersections, unprocessedIntersections, streets)
            unprocessedIntersections.remove(_id)
    
    print("\n")


def processBundleIntersection(intersection, intersections, unprocessedIntersections, streets):
    connector = intersection.startConnector
    while True:
        street = connector.item
        if not street.id in streets:
            streets[street.id] = street
            oppositeConnector = street.pred if street.succ is connector  else street.succ
            if not oppositeConnector.id in intersections:
                intersections[oppositeConnector.id] = oppositeConnector.intersection
                unprocessedIntersections.add(oppositeConnector.intersection.id)
        connector = connector.succ
        if connector is intersection.startConnector:
            break