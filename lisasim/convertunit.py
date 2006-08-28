# $Id$
# $Date$
# $Author$
# $Revision$

# so far we can convert
#
# - Degree to Radian
# - Degree/Minute/Second (DMS), space-separated, to Radian
# - Hour/Minute/Second (HMS), space-separated, to Radian

import math

def convertUnit(param,unitin,unitout,paramname=''):
    if unitout == unitin:
        return (param,unitin)

    if unitout == "Radian":
        if unitin == "Degree":
            return ((math.pi/180.0) * float(param),'Radian')
        elif unitin == "DMS":
            # space separated
            dms = map(float, param.split())

            return ((math.pi/180.0) * (dms[0] + dms[1]/60.0 + dms[2]/3600.0),'Radian')
        elif unitin == "HMS":
            # space separated
            hms = map(float, param.split())
            
            return ((15.0*math.pi/180.0) * (hms[0] + hms[1]/60.0 + hms[2]/3600.0),'Radian')

    if unitout == 'String':
        if paramname == 'SpectralType':
            if unitin == '1':
                if float(param) == 0.0:
                    return ('WhiteFrequency','String')
                elif float(param) == 2.0:
                    return ('WhitePhase','String')
                elif float(param) == -2.0:
                    return ('WhiteAcceleration','String')
                elif float(param) == -4.0:
                    return ('RedAcceleration','String')

    if unitout == '1':
        if paramname == 'SpectralType':
            if unitin == 'String':
                if param == 'WhiteFrequency':
                    return ('0.0','1')
                elif param == 'WhitePhase':
                    return ('2.0','1')
                elif param == 'WhiteAcceleration':
                    return ('-2.0','1')                
                elif param == 'RedAcceleration':
                    return ('-4.0','1')

    raise NotImplementedError, "convertUnit(): cannot convert %s from %s to %s (parameter %s)" % (param,unitin,unitout,paramname)

# conversion rules (a tuple consisting of a list of required parameters with their
# units, and a function returning the transformed parameter with its unit

conversionRules = {}

conversionRules['TimeOffset'] =     ( [('Prebuffer','Second')],
                                      lambda x: ( str(-float(x)), 'Second' ) ) 

conversionRules['Prebuffer'] =      ( [('TimeOffset','Second')],
                                      lambda x: ( str(-float(x)), 'Second' ) ) 

conversionRules['Polarization'] =   ( [('SLPolarization','Radian')],
                                      lambda x: ( float(x) != 0 and str(2.0*math.pi - float(x)) or '0', 'Radian' ) )

conversionRules['SLPolarization'] = ( [('Polarization','Radian')],
                                      lambda x: ( float(x) != 0 and str(2.0*math.pi - float(x)) or '0', 'Radian' ) )

conversionRules['Inclination'] =    ( [('ThetaInclination','Radian')],
                                      lambda x: ( str(float(x) < math.pi and
                                                  math.pi - float(x) or
                                                  3.0*math.pi - float(x)), 'Radian' ) )
                                      
conversionRules['ThetaInclination'] = ( [('Inclination','Radian')],
                                        lambda x: ( str(float(x) < math.pi and
                                                    math.pi - float(x) or
                                                    3.0*math.pi - float(x)), 'Radian' ) )

def eq2ec(a,d):
    alpha = float(a)
    delta = float(d)

    # J2000 epoch
    epsilon = (math.pi/180.0) * 23.439291

    beta = math.asin( math.cos(epsilon) * math.sin(delta) -
                      math.sin(epsilon) * math.cos(delta) * math.sin(alpha) )
                      
    coslambd = math.cos(alpha) * math.cos(delta) / math.cos(beta)
    sinlambd = ( (math.sin(delta) - math.cos(epsilon) * math.sin(beta)) /
                 (math.sin(epsilon) * math.cos(beta)) )
                 
    lambd = math.atan2(sinlambd,coslambd)
    if lambd < 0:
        lambd = lambd + 2.0*math.pi

    return ( (str(beta),'Radian'), (str(lambd),'Radian') )

conversionRules['EclipticLatitude'] = ( [('RightAscension','Radian'),('Declination','Radian')],
                                        lambda x,y: eq2eq(x,y)[0] )

conversionRules['EclipticLongitude'] = ( [('RightAscension','Radian'),('Declination','Radian')],
                                         lambda x,y: eq2eq(x,y)[1] )

def ex2ipr(e,x,s):
    eta = float(e)
    xi = float(x)
    sw = float(s)

    if sw > 0:
        raise NotImplementedError, "convertParameters(): LISA eta/xi configuration with sw > 0 not compatible with PseudoLISA"

    initpos = eta
    initrot = xi + initpos - 1.5*math.pi

    return ( (str(initpos),'Radian'), (str(initrot),'Radian') )


conversionRules['InitialPosition'] = ( [('InitialEta','Radian'),('InitialXi','Radian'),('ArmSwitch','1')],
                                       lambda x,y,z: ex2ipr(x,y,z)[0] )

conversionRules['InitialRotation'] = ( [('InitialEta','Radian'),('InitialXi','Radian'),('ArmSwitch','1')],
                                       lambda x,y,z: ex2ipr(x,y,z)[1] )


def ipr2ex(p,r):
    kappa = float(p)
    lambd = float(r)

    eta = kappa
    xi = lambd - kappa + 1.5*math.pi    

    return ( (str(eta),'Radian'), (str(xi),'Radian'), ('-1.0','1') )
          
conversionRules['InitialEta'] = ( [('InitialPosition','Radian'),('InitialRotation','Radian')],
                                  lambda x,y: ipr2ex(x,y)[0] )

conversionRules['InitialXi'] =  ( [('InitialPosition','Radian'),('InitialRotation','Radian')],
                                  lambda x,y: ipr2ex(x,y)[1] )

conversionRules['ArmSwitch'] =  ( [('InitialPosition','Radian'),('InitialRotation','Radian')],
                                  lambda x,y: ipr2ex(x,y)[2] )                                  


def interp2length(interptype,interpwindow):
    if interptype == 'NearestNeighbor':
        return ('0','1')
    elif interptype == 'LinearExtrapolator':
        return ('-1','1')
    elif interptype == 'Linear':
        return ('1','1')
    elif interptype == 'Lagrange': 
        return (interpwindow,'1')
    else:
        raise NotImplementedError, "convertParameters(): unknown interpolator type %s" % interptype

conversionRules['InterpolatorLength'] = ( [('Interpolator','String'),('InterpolatorWindow','1')],
                                          lambda x,y: interp2length(x,y) )


def length2interp(l):
    interplen = float(l)
    
    if interplen == 0.0:
        return ( ('NearestNeighbor','String'), ('None','1') )
    elif interplen == -1.0:
        return ( ('LinearExtrapolator','String'), ('None','1') )
    elif interplen == 1.0:
        return ( ('Linear','String'), ('None','1') )
    elif interplen > 1.0:
        return ( ('Lagrange','String'), (str(int(interplen)),'1') )
    else:
        raise NotImplementedError, "convertParameters(): unknown interpolator length %s" % interplen

conversionRules['Interpolator'] = ( [('InterpolatorLength','1')],
                                    lambda x: length2interp(x)[0] )

conversionRules['InterpolatorWindow'] = ( [('InterpolatorLength','1')],
                                          lambda x: length2interp(x)[1] )                                    
                                    
                                    
def convertParameters(param,sourceparams):
    if param[0] in sourceparams:
        return sourceparams[param[0]]

    args = []
    
    try:
        for arg in conversionRules[param[0]][0]:
            try:
                args.append(convertUnit(sourceparams[arg[0]][0],sourceparams[arg[0]][1],arg[1])[0])
            except KeyError:
                raise AttributeError, "convertParameters(): need %s to return %s" % (arg[0],param[0])
    except KeyError:
        raise AttributeError
    
    return conversionRules[param[0]][1](*args)
