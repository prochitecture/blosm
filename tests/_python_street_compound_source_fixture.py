"""Original metre-space controls for a Beijing compound surface regression.

This reduced input deliberately tests the physical compound, not whole-network
Bundle generation. Its missing remote junction context is unsuitable as a full
city import fixture. No generated artifact file is needed to run the test.
"""
from dataclasses import replace
import math
import numpy as np
from tested_python_streets import model


def compound_source_fixture(angle=0., translation=(0., 0.)):
    sources = (
        (397, 397, 10.8, ((833.64599609375, 532.2640380859375),
            (830.8951416015625, 544.7649536132812))),
        (401, 401, 10.8, ((830.8951416015625, 544.7649536132812),
            (831.1597290039062, 545.4774169921875))),
        (402, 402, 10.8, ((830.8951416015625, 544.7649536132812),
            (831.457763671875, 553.30322265625))),
        (406, 406, 5.4, ((831.457763671875, 553.30322265625),
            (817.9664306640625, 552.3109741210938),
            (771.8734741210938, 552.0167846679688),
            (760.4146118164062, 548.8541259765625))),
        (1003, 1003, 2., ((770.3237915039062, 508.3905029296875),
            (777.0947875976562, 512.4097900390625),
            (781.5173950195312, 518.076416015625),
            (783.0279541015625, 526.1027221679688),
            (783.752685546875, 536.1215209960938),
            (783.6492919921875, 545.4166870117188))),
        (25, 2267, 10.8, ((825.8921508789062, 300.2511291503906),
            (825.7032470703125, 309.91363525390625))),
        (29, 2267, 10.8, ((825.7032470703125, 309.91363525390625),
            (823.9221801757812, 355.5766906738281),
            (821.952880859375, 404.82421875),
            (820.5657958984375, 439.94537353515625),
            (820.283203125, 447.203369140625),
            (820.1602172851562, 479.3301696777344),
            (820.0282592773438, 515.3976440429688),
            (830.8951416015625, 544.7649536132812))),
        (398, 2309, 5., ((833.64599609375, 532.2640380859375),
            (843.1497802734375, 532.5433959960938))),
        (405, 2347, 10.8, ((831.457763671875, 553.30322265625),
            (832.1964111328125, 591.93115234375),
            (830.2948608398438, 646.1991577148438))),
        (439, 2349, 5.4, ((760.4146118164062, 548.8541259765625),
            (771.9595336914062, 545.4266357421875),
            (783.6492919921875, 545.4166870117188))),
        (998, 2349, 5.4, ((783.6492919921875, 545.4166870117188),
            (818.2745361328125, 545.3646850585938),
            (831.1597290039062, 545.4774169921875))),
        (32, 2417, 10.8, ((842.3983764648438, 300.74273681640625),
            (841.2329711914062, 337.3667297363281),
            (841.0104370117188, 342.2536315917969),
            (839.7092895507812, 371.4971008300781),
            (839.281005859375, 383.2857971191406),
            (838.9041748046875, 393.5827941894531),
            (836.8914184570312, 448.5966796875))),
        (190, 2417, 10.8, ((836.8914184570312, 448.5966796875),
            (836.1804809570312, 468.67864990234375),
            (834.613037109375, 512.9613647460938),
            (833.64599609375, 532.2640380859375))),
    )
    rotation = np.array(((math.cos(angle), -math.sin(angle)),
                         (math.sin(angle), math.cos(angle))))
    sections = []
    for identity, street, width, points in sources:
        moved = np.asarray(points) @ rotation.T + translation
        sections.append(model.Section(identity, street,
            tuple((float(x), float(y), 0.) for x, y in moved), width=width,
            rise=.07999999821186066))

    def junction(identity, ends, radii):
        return model.Intersection(identity, tuple(model.End(*end) for end in ends),
            radii[0], radii, ('STYLE',) * len(ends))

    def member(*visits):
        return model.BundleMember(tuple(model.Traversal(*visit) for visit in visits))

    return model.Network('original_compound_controls', tuple(sections), (
        junction(802, ((401, 0), (402, 0), (29, 1), (397, 1)), (4., 4., 4., 4.)),
        junction(794, ((398, 0), (397, 0), (190, 1)), (2., 4., 4.)),
        junction(810, ((405, 0), (406, 0), (402, 1)), (4., 2., 4.)),
        junction(2005, ((439, 1), (998, 0), (1003, 1)), (2., 2., 1.))),
        bundles=(model.Bundle(11, (member((25, True), (29, True)),
                                  member((32, True), (190, True)))),
                 model.Bundle(35, (member((406, False)),
                                   member((439, True), (998, True))))),
        connections=(model.Connection(16, (model.End(25, 1), model.End(29, 0))),
                     model.Connection(162, (model.End(998, 1), model.End(401, 1)))),
        settings=replace(model.Settings(), cell_size=.699999988079071, exaggeration=1.))
