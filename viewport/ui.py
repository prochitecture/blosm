"""Provider-independent layout and controls. Coordinates are logical pixels.

No Blender imports: layout, focus, clipping and drag behavior are testable in
ordinary Python. Containers own children; the session owns the single router.
"""
from dataclasses import dataclass

HANDLED = 'HANDLED'
PASS_THROUGH = 'PASS_THROUGH'


@dataclass(frozen=True)
class Rect:
    x: float
    y: float
    width: float
    height: float

    def contains(self, x, y):
        return self.x <= x < self.x+self.width and self.y <= y < self.y+self.height


@dataclass(frozen=True)
class Frame:
    region: Rect
    available: Rect
    scale: float
    revision: int = 0


@dataclass(frozen=True)
class InputEvent:
    type: str
    value: str
    x: float = 0.
    y: float = 0.
    shift: bool = False
    ctrl: bool = False
    alt: bool = False


def shell_layout(rect):
    """Responsive canvas/toolbar/card example shared by current consumers."""
    pad = min(12., rect.width/10, rect.height/10)
    width = max(0., rect.width-2*pad)
    toolbar = Rect(rect.x+pad, rect.y+pad, width, min(78., max(0., rect.height-2*pad)))
    card = (Rect(rect.x+rect.width-232-pad, rect.y+112, 232., max(0., rect.height-140))
            if rect.width >= 850 and rect.height >= 360 else None)
    return toolbar, card


class Widget:
    focusable = False
    def __init__(self, key, rect=None):
        self.key, self.rect = key, rect or Rect(0,0,0,0)
        self.visible, self.enabled = True, True

    def handle(self, event, router):
        return PASS_THROUGH

    def draw(self, painter, router):
        pass

    def dispose(self):
        pass


class Button(Widget):
    focusable = True
    def __init__(self, key, label, action):
        super().__init__(key)
        self.label, self.action = label, action

    def handle(self, event, router):
        if event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            router.capture = self
            return HANDLED
        if event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            captured = router.capture is self
            router.capture = None
            if captured and self.rect.contains(event.x,event.y):
                self.action()
            return HANDLED
        if event.type in {'RET','SPACE'} and event.value == 'PRESS':
            self.action()
            return HANDLED
        return PASS_THROUGH

    def draw(self, painter, router):
        color = (.22,.32,.42,1) if router.focus is self or router.hover is self else (.13,.17,.22,1)
        painter.rect(self.rect, color if self.enabled else (.1,.1,.1,1))
        painter.text(self.label,self.rect.x+8,self.rect.y+8,max_width=max(0,self.rect.width-16))


class Slider(Widget):
    focusable = True
    def __init__(self,key,label,get_value,set_value):
        super().__init__(key)
        self.label,self.get_value,self.set_value = label,get_value,set_value

    def handle(self,event,router):
        if event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            router.capture = self
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            router.capture = None
            return HANDLED
        elif event.type == 'MOUSEMOVE' and router.capture is self:
            pass
        elif event.type in {'LEFT_ARROW','RIGHT_ARROW'} and event.value == 'PRESS':
            self.set_value(max(0.,min(1.,self.get_value()+(.05 if event.type=='RIGHT_ARROW' else -.05))))
            return HANDLED
        else:
            return PASS_THROUGH
        self.set_value(max(0.,min(1.,(event.x-self.rect.x)/max(1.,self.rect.width))))
        return HANDLED

    def draw(self,painter,router):
        painter.rect(self.rect,(.10,.13,.17,.95))
        painter.rect(Rect(self.rect.x,self.rect.y,self.rect.width*self.get_value(),4),(.25,.65,.95,1))
        painter.text('%s %d%%' % (self.label,round(100*self.get_value())),self.rect.x+8,self.rect.y+9,
                     max_width=max(0,self.rect.width-16))


class Container(Widget):
    """Ordered composition. Later children draw/hit-test above earlier ones."""
    def __init__(self,key,children=()):
        super().__init__(key)
        self.children = list(children)

    def leaves(self):
        for child in self.children:
            if child.visible:
                if isinstance(child,Container):
                    yield from child.leaves()
                else:
                    yield child

    def hit(self,x,y):
        if not self.visible or not self.rect.contains(x,y):
            return None
        for child in reversed(self.children):
            if not child.visible or not child.enabled:
                continue
            if isinstance(child,Container):
                hit = child.hit(x,y)
                if hit is not None: return hit
            elif child.rect.contains(x,y):
                return child
        return None

    def draw(self,painter,router):
        with painter.clip(self.rect):
            for child in self.children:
                if child.visible: child.draw(painter,router)

    def dispose(self):
        for child in self.children: child.dispose()
        self.children.clear()


class Router:
    def __init__(self,root):
        self.root = root
        self.capture = self.focus = self.hover = None

    def reset_capture(self):
        self.capture = self.hover = None
        if self.focus not in tuple(self.root.leaves()): self.focus = None

    def handle(self,event):
        if event.type == 'TAB' and event.value == 'PRESS':
            controls = [w for w in self.root.leaves() if w.focusable and w.enabled]
            if controls:
                index = controls.index(self.focus) if self.focus in controls else (-1 if not event.shift else 0)
                self.focus = controls[(index+(-1 if event.shift else 1)) % len(controls)]
            return HANDLED
        self.hover = self.root.hit(event.x,event.y)
        target = self.capture
        if target is None and event.type in {'RET','SPACE','LEFT_ARROW','RIGHT_ARROW'}:
            target = self.focus
        if target is None: target = self.hover
        if target is not None:
            if event.type=='LEFTMOUSE' and event.value=='PRESS': self.focus = target
            return target.handle(event,self)
        return PASS_THROUGH
