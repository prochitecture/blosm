"""Finite intersections of two leads on the same original physical curb."""
import numpy as np


def same_bank(first,last):
    """Explicit native banks require their exact authored offset as well."""
    if first['path'] is not last['path']:return False
    if first.get('explicit_bank') or last.get('explicit_bank'):
        return first['offset']==last['offset']
    return abs(first['offset']-last['offset'])<=1e-8


def _intervals(curve,lower,upper,initial):
    from .source_curb_coordinates import original_bank_intervals
    native=original_bank_intervals(curve,lower,upper,initial)
    if native:return native
    leads=getattr(curve,'leads',None)
    lengths=getattr(curve,'lengths',None)
    if leads is None or lengths is None or len(lengths)!=4:return ()
    result=[]
    for part,index in ((0,0),(2,1)):
        start=max(lower*curve.length,float(lengths[part]))
        end=min(upper*curve.length,float(lengths[part+1]))
        if end<=start+1e-12:continue
        approach,parameters,cumulative,distance=leads[index]
        if cumulative[-1]<=1e-12:continue
        a,b=np.interp(np.array((start,end))-lengths[part],cumulative,parameters)
        lo,hi=sorted((float(a),float(b)))
        seed=float(np.interp(np.clip(initial*curve.length-lengths[part],0.,cumulative[-1]),cumulative,parameters))
        path=approach.path
        if not getattr(path,'outward',False):continue
        for j,piece in enumerate(path.pieces):
            base=float(path.stations[j]);length=piece.path.length
            first,last=max(lo,base),min(hi,base+length)
            if last<=first+1e-12:continue
            global_station=lambda s: s-base if piece.forward else length-(s-base)
            low,high=sorted((global_station(first),global_station(last)))
            ordered_parameters,ordered_cumulative=(parameters,cumulative) if parameters[-1]>=parameters[0] else (parameters[::-1],cumulative[::-1])
            result.append(dict(path=piece.path,offset=distance*(1. if piece.forward else -1.),
                low=low,high=high,seed=global_station(np.clip(seed,first,last)),
                forward=piece.forward,base=base,length=length,parameters=ordered_parameters,
                cumulative=ordered_cumulative,curve_base=float(lengths[part]),curve_length=curve.length))
    return result


def contact(curves,initial,lower,upper):
    """Choose a common source station, never a proximity-based connection.

    A Newton intersection of coincident finite leads has a singular Jacobian.
    The identical original Path and oriented offset identify the common curve;
    only their bounded overlapping station domains can supply a contact.
    """
    intervals=[_intervals(c,float(lo),float(hi),float(seed))
               for c,lo,hi,seed in zip(curves,lower,upper,initial)]
    candidates=[]
    for a in intervals[0]:
        for b in intervals[1]:
            if not same_bank(a,b):continue
            lo,hi=max(a['low'],b['low']),min(a['high'],b['high'])
            if hi<lo:continue
            station=float(np.clip(.5*(a['seed']+b['seed']),lo,hi))
            values=[]
            for row in (a,b):
                outward=row['base']+(station if row['forward'] else row['length']-station)
                arc=float(np.interp(outward,row['parameters'],row['cumulative']))
                values.append((row['curve_base']+arc)/row['curve_length'])
            values=np.asarray(values)
            if np.any(values<lower-1e-9) or np.any(values>upper+1e-9):continue
            if np.linalg.norm(curves[0](values[0])-curves[1](values[1]))>1e-7:continue
            candidates.append((float(np.sum((values-initial)**2)),values))
    return min(candidates,key=lambda item:item[0])[1] if candidates else None


def source_station(row,parameter):
    """Map one proven lead parameter back to its original Section station."""
    arc=parameter*row['curve_length']-row['curve_base']
    xp,fp=row['cumulative'],row['parameters']
    if xp[-1]<xp[0]:xp,fp=xp[::-1],fp[::-1]
    outward=float(np.interp(arc,xp,fp))-row['base']
    return outward if row['forward'] else row['length']-outward


def curve_parameter(row,station):
    """Invert source_station using the same original finite lead table."""
    outward=row['base']+(station if row['forward'] else row['length']-station)
    arc=float(np.interp(outward,row['parameters'],row['cumulative']))
    return (row['curve_base']+arc)/row['curve_length']
