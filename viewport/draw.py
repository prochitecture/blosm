"""Small GPU painter. All methods require Blender's main graphics context."""
from contextlib import contextmanager
import blf
import gpu
from gpu_extras.batch import batch_for_shader
from .ui import Rect

DISPLAY_RGB_GLSL = '''
vec3 display_rgb(vec3 c) {
    c=max(c,vec3(0.0));
    return mix(12.92*c,1.055*pow(c,vec3(1.0/2.4))-0.055,step(vec3(0.0031308),c));
}
'''


class Painter:
    def __init__(self):
        self.shader = gpu.shader.from_builtin('UNIFORM_COLOR')
        self.batch = batch_for_shader(self.shader,'TRI_FAN',{'pos':((0,0),(1,0),(1,1),(0,1))})
        self.frame = None
        self.clip_rect = None
        info=gpu.types.GPUShaderCreateInfo()
        info.vertex_in(0,'VEC2','pos')
        info.push_constant('VEC4','bounds')
        info.push_constant('FLOAT','opacity')
        info.push_constant('INT','premultiplied')
        info.push_constant('VEC4','clip_rect')
        info.sampler(0,'FLOAT_2D','image_texture')
        interface=gpu.types.GPUStageInterfaceInfo('blosm_image_interface')
        interface.smooth('VEC2','uv')
        info.vertex_out(interface)
        info.fragment_out(0,'VEC4','fragColor')
        info.vertex_source('void main(){uv=pos;gl_Position=vec4((bounds.xy+pos*bounds.zw)*2.0-1.0,0.0,1.0);}')
        info.fragment_source(DISPLAY_RGB_GLSL+'''
        void main(){
            if(any(lessThan(gl_FragCoord.xy,clip_rect.xy)) || any(greaterThanEqual(gl_FragCoord.xy,clip_rect.zw))) discard;
            vec4 c=texture(image_texture,uv);
            if(premultiplied==1 && c.a>0.00001)c.rgb/=c.a;
            fragColor=vec4(display_rgb(c.rgb),c.a*opacity);}
        ''')
        self.image_shader=gpu.shader.create_from_info(info)
        self.image_batch=batch_for_shader(self.image_shader,'TRI_FAN',{'pos':((0,0),(1,0),(1,1),(0,1))})

    @contextmanager
    def state(self,frame):
        self.frame = frame
        blend,depth,mask = gpu.state.blend_get(),gpu.state.depth_test_get(),gpu.state.depth_mask_get()
        self.clip_rect=frame.region
        try:
            gpu.state.blend_set('ALPHA')
            gpu.state.depth_test_set('NONE')
            gpu.state.depth_mask_set(False)
            with self.clip(frame.available): yield self
        finally:
            gpu.state.depth_mask_set(mask)
            gpu.state.depth_test_set(depth)
            gpu.state.blend_set(blend)

    @contextmanager
    def clip(self,rect):
        # Blender exposes no scissor-enable getter. Clip geometry/shader pixels
        # instead of changing an enable bit that cannot safely be restored.
        old=self.clip_rect
        left,bottom=max(old.x,rect.x),max(old.y,rect.y)
        right,top=min(old.x+old.width,rect.x+rect.width),min(old.y+old.height,rect.y+rect.height)
        self.clip_rect=Rect(left,bottom,max(0,right-left),max(0,top-bottom))
        try: yield
        finally:self.clip_rect=old

    def pixel_clip(self):
        x,y,_,_=gpu.state.viewport_get()
        r,s=self.clip_rect,self.frame.scale
        return (x+r.x*s,y+r.y*s,x+(r.x+r.width)*s,y+(r.y+r.height)*s)

    def rect(self,rect,color):
        s = self.frame.scale
        r=self.clip_rect
        left,bottom=max(r.x,rect.x),max(r.y,rect.y)
        right,top=min(r.x+r.width,rect.x+rect.width),min(r.y+r.height,rect.y+rect.height)
        rect=Rect(left,bottom,max(0,right-left),max(0,top-bottom))
        with gpu.matrix.push_pop():
            gpu.matrix.translate((rect.x*s,rect.y*s,0))
            gpu.matrix.scale((rect.width*s,rect.height*s,1))
            self.shader.bind()
            self.shader.uniform_float('color',color)
            self.batch.draw(self.shader)

    def text(self,value,x,y,*,size=13,color=(.94,.96,1,1),max_width=None):
        s = self.frame.scale
        r=self.clip_rect
        if x<r.x or y<r.y+2 or y+size>r.y+r.height:return
        available=max(0,r.x+r.width-x)
        max_width=available if max_width is None else min(available,max_width)
        value = str(value).replace('\n',' ')
        blf.size(0,size*s)
        if max_width is not None:
            # Text fields are deliberately bounded (metadata may be long).
            value = value[:512]
            if blf.dimensions(0,value)[0] > max_width*s:
                while value and blf.dimensions(0,value+'…')[0] > max_width*s: value=value[:-1]
                value += '…' if value else ''
        blf.position(0,x*s,y*s,0)
        blf.color(0,*color)
        blf.draw(0,value)

    def image(self,reference,rect,*,fit='contain',opacity=1.):
        """Draw a borrowed ImageRef; preserve aspect unless stretch is explicit."""
        if fit not in {'contain','cover','stretch'}: raise ValueError('Unknown image fit')
        target=rect
        if fit!='stretch':
            width,height=reference.size
            ratio=(min if fit=='contain' else max)(rect.width/width,rect.height/height)
            target=Rect(rect.x+(rect.width-width*ratio)/2,rect.y+(rect.height-height*ratio)/2,width*ratio,height*ratio)
        region=self.frame.region
        with self.clip(rect):
            self.image_shader.bind()
            self.image_shader.uniform_float('bounds',(target.x/region.width,target.y/region.height,target.width/region.width,target.height/region.height))
            self.image_shader.uniform_float('opacity',opacity)
            self.image_shader.uniform_int('premultiplied',int(reference.premultiplied))
            self.image_shader.uniform_float('clip_rect',self.pixel_clip())
            self.image_shader.uniform_sampler('image_texture',reference.texture)
            self.image_batch.draw(self.image_shader)

    def release(self):
        self.batch = self.shader = None
        self.image_batch = self.image_shader = None
