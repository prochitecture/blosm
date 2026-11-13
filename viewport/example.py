"""Provider-free API example. Invoke manually; never registered at startup.

    from blosm.viewport.example import show
    show(bpy.context)  # from a 3D View context

An optional existing image is borrowed; closing does not delete it. No network,
scene IDs or import transactions are needed by a custom surface.
"""
from . import open
from .session import Surface
from .ui import Container,Button,Slider,Router,Rect,shell_layout
from .assets import ImageRef


class Example(Surface):
    def __init__(self,image=None):
        self.image=image
        self.reference=None
        self.opacity=.8

    def attach(self,session):
        super().attach(session)
        if self.image is not None:self.reference=ImageRef(self.image)
        self.close=Button('close','Close example',session.request_close)
        self.slider=Slider('opacity','Opacity',lambda:self.opacity,self.set_opacity)
        self.root=Container('example',[self.close,self.slider])
        self.router=Router(self.root)

    @property
    def capturing(self):return self.router.capture is not None
    def set_opacity(self,value):self.opacity=value
    def cancel_gesture(self):
        if hasattr(self,'router'):self.router.reset_capture()

    def layout(self,frame):
        self.root.rect=frame.available
        self.toolbar,self.card=shell_layout(frame.available)
        r=self.toolbar
        self.close.rect=Rect(r.x+6,r.y+6,min(130,max(0,r.width-12)),30)
        self.slider.visible=r.width>350
        self.slider.rect=Rect(r.x+150,r.y+6,min(220,max(0,r.width-162)),30)

    def handle(self,event):return self.router.handle(event)

    def draw(self,frame,painter):
        if self.reference:painter.image(self.reference,frame.available,opacity=self.opacity)
        else:painter.rect(frame.available,(.04,.08,.13,self.opacity))
        painter.rect(self.toolbar,(.04,.04,.05,.95))
        self.root.draw(painter,self.router)
        if self.card:
            painter.rect(self.card,(.13,.18,.22,.95))
            painter.text('A reusable side card',self.card.x+12,self.card.y+self.card.height-25,max_width=self.card.width-24)
            painter.text('Collapses in a narrow viewport',self.card.x+12,self.card.y+self.card.height-50,size=11,max_width=self.card.width-24)

    def detach(self,reason):
        if self.reference:self.reference.release()
        if hasattr(self,'root'):self.root.dispose()


def show(context,image=None):return open(context,Example(image))
