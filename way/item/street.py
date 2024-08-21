from .item import Item
from way.item.section import Section


class Street(Item):
    
    ID = 1  # Must not start with zero to get unambiguous connector indices!
    
    def __init__(self, src, dst):
        super().__init__()
        self.id = Street.ID
        Street.ID += 1
        
        self.style = None

        self.src = src
        self.dst = dst
        
        self.head = None
        self.tail = None

        self.pred = None
        self.succ = None

    def insertEnd(self, item):
        item.street = self
        if self.head is None:
            self.head = item 
            self.tail = item
            item.pred = None
            item.succ = None
            return

        # If list is not empty 
        item.pred = self.tail
        item.succ = None
        self.tail.succ = item
        self.tail = item

    def insertStreetEnd(self,newstreet):
        for item in newstreet.iterItems():
            item.street = self

        if self.head is None:
            self.head = newstreet.head
            self.tail = newstreet.tail
            self.dst = newstreet.dst
            self.src = newstreet.src
            self.head.pred = None
            self.tail.succ = None
            return
        
        # If list is not empty 
        newstreet.head.pred = self.tail
        newstreet.tail.succ = None
        self.tail.succ = newstreet.head
        self.tail = newstreet.tail
        self.dst = newstreet.dst

    def insertFront(self, item):
        item.street = self
        if self.head is None:
            self.head = item 
            self.tail = item
            item.pred = None
            item.succ = None
            return

        # If list is not empty 
        item.succ = self.head
        item.pred = None
        self.head.pred = item
        self.head = item

    def insertStreetFront(self,newstreet):
        for item in newstreet.iterItems():
            item.street = self

        if self.head is None:
            self.head = newstreet.head
            self.tail = newstreet.tail
            self.src = newstreet.src
            self.dst = newstreet.dst
            self.head.pred = None
            self.tail.succ = None
            return

        # If list is not empty 
        newstreet.tail.succ = self.head
        newstreet.head.pred = None
        self.head.pred = newstreet.head
        self.head = newstreet.head
        self.src = newstreet.src

    def iterItems(self):
        if self.head is None:
            yield None
        else:
            current = self.head
            while current is not None:
                yield current
                current = current.succ
            
    def getCategory(self):
        return self.head.getCategory()
    
    def setStyleForItems(self):
        style = self.style
        if self.head is self.tail:
            self.head.setStyleBlockFromTop(style)
        else:
            for item in self.iterItems():
                item.setStyleBlockFromTop(style)
    
    def getName(self):
        return self.head.getName()
    
    def print(self):
        for item in self.iterItems():
            print(item.id, type(item))

    def plot(self,color,width,fwd):
        for item in self.iterItems():
            if isinstance(item,Section):
                
                if fwd:
                    t = item.polyline.d2t(10)
                    item.polyline.trimmed(0,t).plotWithArrows(color,width)
                else:
                    t = item.polyline.d2t(item.polyline.length()-10)
                    item.polyline.trimmed(t,len(item.polyline)).plotWithArrows(color,width)