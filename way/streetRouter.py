from .item.section import Section

def lengthOfStreet(street):
    length = 0.0
    for item in street.iterItems():
        if isinstance(item, Section):
            length += item.polyline.length()
    return length


def categoryOfStreet(street):
    section = street.head
    if isinstance(section, Section):
        return section.category
    return None

# Unit vector pointing from start of the street inwards
def startVec(street):
    firstSection = street.head
    if isinstance(firstSection, Section):
        v1,v2 = firstSection.polyline[0], firstSection.polyline[1]
        return (v2-v1)/(v2-v1).length
    return None

# Unit vector pointing from end of the street outwards
def endVec(street):
    lastSection = street.tail
    if isinstance(lastSection, Section):
        v1,v2 = lastSection.polyline[-2], lastSection.polyline[-1]
        return (v2-v1)/(v2-v1).length
    return None

# Tries to find the street segments, that fill tha gab between <fromStreet> 
# and <toStreet>. The angle between the vector at the end of <fromStreet>
# the one at the start of <toStreet> needs to be less than 60°.
def streetsOnRoute(fromStreet, toStreet, dist_max, max_angle=0.9):
    fromVec = endVec(fromStreet)
    toVec = startVec(toStreet)
    cosAngle = fromVec.dot(toVec)
    if cosAngle < 0.5:  # corresponds to 60°
        return []
    
    category = fromStreet.head.category
    goalIntersection = toStreet.pred.intersection

    currFromVec = fromVec
    currIntersection = fromStreet.succ.intersection
    currStreet = fromStreet
    routedStreets = []
    while True:
        # Find the leaving street of the current intersection, that
        # has the same category and that has an angle that is not very sharp.
        found = False
        for way in currIntersection.leaveWays:
            if way.leaving:
                leavingStreet = way.street
                if leavingStreet == currStreet:
                    continue
                if categoryOfStreet(leavingStreet) == category:
                    currToVec = startVec(leavingStreet)
                    cosAngle = currFromVec.dot(currToVec)
                    if cosAngle > max_angle:
                        found = True
                        break

        # Only short streets are allowed to fill a gap
        if lengthOfStreet(leavingStreet) > dist_max:
            found = False
            return []

        # if theere is a leaving way, continue with the intersection
        # at its end.
        if found and leavingStreet.succ:
            routedStreets.append(leavingStreet)
            currIntersection = leavingStreet.succ.intersection
            if currIntersection == goalIntersection:
                break
            currStreet = leavingStreet
            currVector = endVec(leavingStreet)
        else:
            return []

    return routedStreets

    
