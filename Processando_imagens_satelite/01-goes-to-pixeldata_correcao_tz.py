#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: andre
"""

# Arguments: [1]folderpath with nc files; [2]:output folder; [3] site code; [4] boxsize; [5] plot phase (Yes or No)
# Usage: 

import os
import sys
import glob
import numpy as np
import scipy.ndimage
import datetime as dt
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pylab as py
import calendar
import warnings

############################### Importando módulos ##########################################
from physconstants import *
from goesconstants import *

############################### Definindo funções #############################################
# Primeira função a ser chamada no loop principal - Define os sites
def getsite (lesite):
    sitelist = np.array(['AF','AT','RB','PR'])
    failoc=len(np.where(sitelist == lesite)[0])
    if (failoc == 0):
        print 'Failed to recognize the site ',lesite
        print 'Possible sites are: ',sitelist
        print 'Exiting now.'
        sys.exit()
    # Define the site location and time zone
    site,mylat,mylon,tz= 'AF',-9.87083333,-56.10388889, -4     #Alta Floresta
    if lesite == 'AT':
        site,mylat,mylon,tz='AT',-2.13333333,-59.0000000, -4   #Torre ATTO
    if lesite == 'RB':
        site,mylat,mylon,tz='RB',-09.956944, -67.868889,-4     #Rio Branco 
    if lesite == 'PR':
        site,mylat,mylon,tz='PR',-1.00000000,-68.08000000, -4  #Site "Pristino" - "PR"  
    return site,mylat,mylon,tz

# Segunda função a ser chamada no loop principal - Constantes de calibração do satélite antes dele ser lançado - vêm do goesconstants.py 
def getprelaunchconstant(fn):
    mVIS,kVIS,C,m2,m3,m4,m6,bb2,bb3,bb4,bb6,factor,F0=mVIS13,kVIS13,C13,m213,m313,m413,m613,bb213,bb313,bb413,bb613,factor13,F0_13
    n2,n3,n4,n6,a2,a3,a4,a6,b2,b3,b4,b6,g2,g3,g4,g6=n213,n313,n413,n613,a213,a313,a413,a613,b213,b313,b413,b613,g213,g313,g413,g613
    SAT=fn.split('/')[len(fn.split('/'))-1].split('.')[0]
    year=fn.split('/')[len(fn.split('/'))-1].split('.')[1]
    julian=int(fn.split('/')[len(fn.split('/'))-1].split('.')[2])
    Bfunction,satlat,satlon=Bfunction13,satlat13,satlon13
    return mVIS,kVIS,C,m2,m3,m4,m6,bb2,bb3,bb4,bb6,factor,F0,n2,n3,n4,n6,a2,a3,a4,a6,b2,b3,b4,b6,g2,g3,g4,g6,Bfunction,satlat,satlon

# Terceira função a ser chamada no loop principal - Constantes de calibração do satélite após entrar em órbita
def getpostlaunchconstant(name):
    fn=name.split('/')[len(name.split('/'))-1]  # example filename fn='goes08.2000.001.171515.BAND_01.nc' 
    mysat    = fn.split('.')[0]                 # pega só o 'goes08' no fn
    myyear   = int(fn.split('.')[1])            # pega só o '2000' no fn
    myjulian = int(fn.split('.')[2])            # pega só o '001' no fn
    g8initdate,g12initdate,g13initdate=dt.datetime(1995,6,9),dt.datetime(2003,4,1),dt.datetime(2010,4,14)
    g8alpha,g8beta,g12alpha,g12beta=1.2366,0.0535,1.0848,0.049
    p0,p1,p2,p3=0.0609347,0.110914,-0.497612,0.115607    
    refdate=dt.datetime(myyear-1,12,31)+dt.timedelta(myjulian)
    myc=1.0      
    if mysat == 'goes13':
        initdate=g13initdate
        mydelta=refdate-initdate
        mydelta=mydelta.total_seconds()/(365.25*24*60*60)
        myc=p0+p1*np.power(mydelta,p2)+np.power(mydelta,p3)
    return myc    

# Quarta função a ser chamada no loop principal - Calcula Refletâncias para os canais 1 (VIS) e 2 (NIR) e também SZA e VZA   
def reflectance ( file, Temp='None', Bfunction='None'):
    band=file.split('.')[4]     # pega só o 'BAND_01' no fn (que é igual ao "file" aqui) 
    procmode='IR'
    if band == 'BAND_01':
        procmode = 'VIS'
    print 'Reading: ', file
    print 'Mode: ', procmode
    fh = Dataset(file, mode='r') #Lendo o netCDF cru, ou seja, lendo o file=goes08.2000.001.171515.BAND_01.nc 
    lon = fh.variables['lon'][:]  #Pegando a Longitude do arquivo .nc
    lat = fh.variables['lat'][:]   #Pegando a Latitude do arquivo .nc
    data = fh.variables['data'][:] #Pegando todos os outros dados do arquivo .nc
    fh.close()     #Fechando o netcdf
    Raw=data/32.   #Convertendo os dados de 16bit para 10bit -> é o "X" no passo 1 de http://www.ospo.noaa.gov/Operations/GOES/calibration/gvar-conversion.html
    year=file.split('.')[1] # pega só o '2000' no fn
    julian=file.split('.')[2] # pega só o '001' (o dia juliano) no fn
    hhmm=file.split('.')[3]   # pega só o '171515' (o horário UTC) no fn
    h=dt.datetime(int(year), 1, 1) + dt.timedelta(int(julian) - 1)
    date=dt.date(h.year,h.month,h.day)
    time=dt.time(int(hhmm[:-4]),int(hhmm[2:-2]),int(hhmm[4:]))   
    SZA, SED = sunzen(date, time, tz, lat, lon) # Chamando a função sunzen
    VZA = viewzen(lat,lon) # Chamando a função viewzen
    mu0 = np.cos(np.deg2rad(SZA)) # mizero = cos(SZA), com SZA sendo convertido de graus para radianos     
    if procmode == 'VIS':
        RefPre=kVIS*(Raw-X0)  #Cálculo de refletância no visivel antes do satélite ser lançado
        RefVIS=RefPre * C     #Cálculo de refletância no visivel com o satélite em órbita, já levando em consideração a constante de degradação C   
        RefVIS=RefVIS*(SED*SED/mu0)  # Refletancia final no VIS - NOAA equation for reflectance, considering negligible path radiance: can be wrong if lots of aerosols above clouds
        cs = RefVIS
    if procmode == 'IR':
        RadIR=(Raw - bb2)/m2 #Cálculo de radiância no NIR (3.75) antes do satélite ser lançado (fórmula em "http://www.ospo.noaa.gov/Operations/GOES/calibration/gvar-conversion.html")
        RadIR=RadIR*factor  #Converting mW/[m2.sr.cm-1] to W/[m2.sr.um] considering GOES-* Channel 2 spectral response function - reposta espectral do canal 2
        Emission = getemission(Temp,Bfunction) #Pegando a emissão de corpo negro do solo, à temperatura Temp e em 3,90 micrômetros - é o B3,7(T) nas fórmulas (3) de Kaufmann-Nakajima_93
        RefIR = (RadIR - Emission)  /     ((t0 * F0 * mu0 / (np.pi*SED*SED)) - Emission) #Cálculo da refletância em 3,9 -> fórmulas (3) de Kaufmann-Nakajima_93
        cs = RefIR
    return cs,SZA,VZA 
    
# Primeira função chamada dentro da função "reflectance" - calculate the sun zenith angle (SZA) and the Sun-Earth distance (SED)   
def sunzen ( date, time, tz, lat, lon ):
    #function to calculate the sun zenith angle and the Sun-Earth distance
    #input: date object such as 2013-02-25 (UTC), time object such as 14:00:00 (UTC)
    #input: time zone such as -4 for Manaus
    #input: lat, lon such as -3.0, -60.0 for Manaus
    #output: sun zenith angle in degrees
    #output: Sun-Earth distance in astronomical units

    ############################### MODIFICAÇÂO IMPLEMENTADA #################################
    tz = lon/15.0                                                                            #
    tz = tz.astype(int)                                                                      #
    ############################### MODIFICAÇÂO IMPLEMENTADA #################################

    days=date-dt.date(1900,1,1)
    days=days.days+2
    time=((time.hour+tz)/24.)+(time.minute/60./24.)+(time.second/60./60./24.)
    
    # Calculate Julian Day and Julian Century
    JD = days + 2415018.5 + time - tz/24.
    JC = (JD - 2451545.0)/36525.0
    
    # Calculate Geometric Mean Anomaly Sun (deg) and Geometric Mean Lon Sun (deg)
    GMAS = 357.52911 + JC*(35999.05029 - (0.0001537)*JC)
    GMLS = np.mod(280.46646 + JC*(36000.76983 + JC*0.0003032),360)
    
    # Calculate Sun Eq of Ctr, Sun True Long (deg),  Sun True Anom (deg) and Sun App Long (deg)
    SEC = np.sin(np.deg2rad(GMAS))*(1.914602-JC*(0.004817+0.000014*JC)) + np.sin(np.deg2rad(2*GMAS))*(0.019993-0.000101*JC) + np.sin(np.deg2rad(3*GMAS))*0.000289
    STL = GMLS + SEC
    STA = GMAS + SEC
    SAL = STL - 0.00569 - 0.00478 * np.sin(np.deg2rad(125.04 - 1934.136*JC))
    
    # Calculate Mean Obliq Ecliptic (deg), Obliq Correction (deg), Eccent Earth Orbit
    MOE = 23 + (26 + ((21.448 - JC*(46.815 + JC*(0.00059 - JC*0.001813))))/60)/60
    OC = MOE + 0.00256*np.cos(np.deg2rad(125.04 - 1934.136*JC))
    EEO = 0.016708634 - JC *(0.000042037 + 0.0000001267*JC)
        
    #Calculate the Equation of Time (min) and the True Solar Time (min)
    vary = np.tan(np.deg2rad(OC/2))*np.tan(np.deg2rad(OC/2))
    EOT = 4*np.rad2deg(vary*np.sin(2*np.deg2rad(GMLS)) - 2*EEO*np.sin(np.deg2rad(GMAS)) + 4*EEO*vary*np.sin(np.deg2rad(GMAS)) * np.cos(2*np.deg2rad(GMLS)) - 0.5*vary*vary*np.sin(4*np.deg2rad(GMLS)) - 1.25*EEO*EEO*np.sin(2*np.deg2rad(GMAS))    )
    TST = np.mod(time*1440 + EOT + 4*lon -60*tz, 1440)
    
    # Calculate the Hour Angle (deg) and the Sun Declination (deg)
    HA = TST/4 - 180
    id=TST<0
    HA[id]=TST[id]/4 + 180
    SD = np.rad2deg(np.arcsin(np.sin(np.deg2rad(OC))*np.sin(np.deg2rad(SAL))))
    
    #Calculate Sun-Earth distance (AUs)
    SED = (1.000001018 * (1 - EEO * EEO))/(1 + EEO*np.cos(np.deg2rad(STA)))    
    
    #Calculate the Sun Zenith Angle (deg) - Ângulo entre o zênite de um observador na Terra e o satélite = SZA (ângulo solar zenital)
    SZA = np.rad2deg(np.arccos(np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(SD)) + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(SD)) *np.cos(np.deg2rad(HA))  ))
    
    return SZA, SED  #### Retorna o ângulo zenital solar "SZA" e a distância Terra-Sol "Sun-Earth Distance, ou SED"

# Segunda função chamada dentro da função "reflectance" - calculate the satellite view zenith angle (VZA) - usa constantes do physconstants.py  
def viewzen ( lat, lon ):
    #function to calculate the satellite view zenith angle 
    #input: lat, lon such as -3.0, -60.0 for Manaus
    #output: sattelite zenith angle in degrees

    # Surface polar angle (latitude) and equatorial angle (longitude), in radians
    PA = np.deg2rad(90.0-lat)
    EA = np.deg2rad(lon)

    # Sattelite polar angle (latitude) and equatorial angle (longitude), in radians
    SPA = np.deg2rad(90.0-satlat)
    SEA = np.deg2rad(satlon)
    
    # X,Y,Z projections of R1 vector, from Earth center to surface pixel
    R1X = ER*np.sin(PA)*np.cos(EA)  # ER = Earth radius in km, defined in physconstants.py
    R1Y = ER*np.sin(PA)*np.sin(EA)
    R1Z = ER*np.cos(PA)

    # X,Y,Z projections of R2 vector, from Earth center to sattelite
    R2X = SR*np.sin(SPA)*np.cos(SEA)  # SR = Sattelite orbit radius in km, defined in physconstants.py
    R2Y = SR*np.sin(SPA)*np.sin(SEA)
    R2Z = SR*np.cos(SPA)

    # X,Y,Z projections of R3 vector, from surface to sattelite
    R3X = R2X - R1X  
    R3Y = R2Y - R1Y
    R3Z = R2Z - R1Z
    
    # X,Y,Z projections of U, unitary vertical surface vector
    UX = np.sin(PA)*np.cos(EA)  
    UY = np.sin(PA)*np.sin(EA)
    UZ = np.cos(PA)
    
    # Calculate the cosine of VZA, then VZA in degrees
    CVZA = (R3X*UX + R3Y*UY + R3Z*UZ)/np.sqrt(R3X*R3X + R3Y*R3Y + R3Z*R3Z)
    VZA  = np.rad2deg(np.arccos(CVZA))
    
    return VZA    #### Retorna o ângulo zenital de visada do satélite, ou 'sattelite view zenith angle = VZA'

# Terceira função chamada dentro da função "reflectance" - Retorna o B3,7(T) usado nas fórmulas (3) de Kaufmann-Nakajima_93   
def getemission(Temp,Bfunction): # de onde vem Temp?
    myB=-1.0*Temp
    myB[myB<0]=mynan #-9999
    l0,l1,l2,l3,l4,h0,h1,h2,h3,h4,h5,h6=Bfunction[0],Bfunction[1],Bfunction[2],Bfunction[3],Bfunction[4],Bfunction[5],Bfunction[6],Bfunction[7],Bfunction[8],Bfunction[9],Bfunction[10],Bfunction[11]
    lotemp=(Temp>154.0)&(Temp<=190.0)
    hitemp=(Temp>190.0)
    myB[lotemp]=l0+l1*(np.power(Temp[lotemp],l2))+l3*Temp[lotemp]+l4*(np.power(Temp[lotemp],2))
    myB[hitemp]=h0+h1*Temp[hitemp]+h2*(np.power(Temp[hitemp],2))+h3*(np.power(Temp[hitemp],3))+h4*(np.power(Temp[hitemp],4))+h5*(np.power(Temp[hitemp],5))+h6*(np.power(Temp[hitemp],6))
    return myB

# Quinta função a ser chamada no loop principal - Temp. brilho canal 4  
def temperature ( file ):
    procmode='IR'
    print 'Reading: ', file
    print 'Mode: ', procmode
    fh = Dataset(file, mode='r')
    data = fh.variables['data'][:]
    fh.close()
    Raw=data/32.   #Convert 16bit to 10bit
    RadIR=(Raw - bb4)/m4
    Teff4=RadIR
    Teff4=Teff4-RadIR+np.float(mynan) #-9999.
    TIR=Teff4
    Teff4=(C2 * n4)/(np.log(1.+(C1 * np.power(n4,3.) /RadIR)))
    TIR=a4+b4*Teff4+g4*( np.power(Teff4,2.) )
    return TIR # Retorna a temperatura de brilho (Canal 4, far infra-red)

# Sexta função a ser chamada no loop principal
def checkvisshape (v,t):
    v,t=np.squeeze(v),np.squeeze(t)
    newV=v
    if (np.shape(v)!=np.shape(t)):
        [av1,av2],[dt1,dt2]=np.shape(v),np.shape(t) # assuming temp matrix has correct dimensions, that v needs to match
        dv1,dv2=4*dt1,4*dt2  # desired dimensions for v matrix
        if (dv1 < av1): # if desired smaller than actual, just trim the matrix
            newV=newV[0:dv1,:]
        if (dv2 < av2):
            newV=newV[:,0:dv2]
        if (dv1 > av1): # if desired bigger than actual, pad with -9999 (mynan)
            current2=np.shape(newV)[1]
            newV=np.vstack((newV,np.zeros((dv1-av1,current2))+mynan))
        if (dv2 > av2):
            current1=np.shape(newV)[0]
            newV=np.hstack((newV,np.zeros((current1,dv2-av2))+mynan))
    return newV

# Sétima função a ser chamada no loop principal - Gera o arquivo "Cloud_Fraction.csv"
def cloudfraction (r063,temp,mylat,mylon,myboxsize,name):
    r063,temp=np.squeeze(r063),np.squeeze(temp)
    lat, lon = getlatlon(name)
    r063=checkvisshape(r063,temp) # after this, r063 will have either same shape as temp or exactly 16x bigger    
    lat=checkvisshape(lat,temp)
    lon=checkvisshape(lon,temp)
    boxfactor=1.0
    [temp,boxfactor]=expandtemp(r063,temp) # this will make sure temperature matrix is resampled to 1km if r063 has that spatial resolution
    dif=np.sqrt(((lat-mylat)*(lat-mylat))+((lon-mylon)*(lon-mylon)))
    loc=np.where(dif==np.min(dif))
    myx,myy=loc[0],loc[1]
    if (np.shape(myx)[0]>1):
        myx=myx[0]
    if (np.shape(myy)[0]>1):
        myy=myy[0]
    llat,llon=lat[myx,myy],lon[myx,myy]
    llat,llon=round(llat,2),round(llon,2)    
    if ((np.abs(llat-mylat)<1.0) and (np.abs(llon-mylon)<1.0) ):
        os.chdir(outdir)
        myname = name.split('/')[len(name.split('/'))-1]          
        boxsize=[int(myboxsize*boxfactor)] #[7,12,20]  Box size around the site, in pixels. Boxsize=7 means [x-7,y-7] to [x+7,y+7], or a square with 15 pixels of side.
        for b in boxsize:
            minx,miny,maxx,maxy=myx-b,myy-b,myx+b,myy+b
            r063b=r063[minx:maxx,miny:maxy]
            tempb=temp[minx:maxx,miny:maxy]
            cloud=(r063b>cloudVISthreshold)&(tempb<cloudTempthreshold)
            noncloud=(r063b<=cloudVISthreshold)|(tempb>=cloudTempthreshold)            
            [d1,d2]=np.shape(r063b)
            totalpixels=np.float(d1)*np.float(d2)
            cloudypixels=np.float(np.sum(cloud))
            noncloudypixels=np.float(np.sum(noncloud))
            CF=mynan #-999
            if totalpixels>0:
                CF = cloudypixels/totalpixels            #print  '{0:.6f}'.format(round(CF,6))
            crefavg,crefstd,ncrefavg,ncrefstd=mynan,mynan,mynan,mynan #-999,-999,-999,-999
            ctmpavg,ctmpstd,nctmpavg,nctmpstd=mynan,mynan,mynan,mynan #-999,-999,-999,-999
            if cloudypixels>1:
                crefavg,crefstd=np.average(r063b[cloud]),np.std(r063b[cloud])
                ctmpavg,ctmpstd=np.average(tempb[cloud]),np.std(tempb[cloud])
            if noncloudypixels>1:
                ncrefavg,ncrefstd=np.average(r063b[noncloud]),np.std(r063b[noncloud])
                nctmpavg,nctmpstd=np.average(tempb[noncloud]),np.std(tempb[noncloud])           
            fh = open(myname[:-11]+'_box'+str(int(b/boxfactor))+'_lat'+str(llat)+'_lon'+str(llon)+'_cloudfraction.csv', 'w')
            ################################# Informações no arquivo "cloudfraction" ##########################################
            outstr=str(np.int(totalpixels))+','+str(np.int(cloudypixels))+','+str(CF)+','+str(crefavg)+','+str(crefstd)+','+str(ctmpavg)+','+str(ctmpstd)+','
            outstr=outstr+str(ncrefavg)+','+str(ncrefstd)+','+str(nctmpavg)+','+str(nctmpstd)+'\n'
            fh.writelines(outstr)
            fh.close()            
    return

# Primeira função chamada dentro da função "cloudfraction"
def expandtemp (v,t):
    factor=1.0
    newtemp=t
    if (np.shape(v)!=np.shape(t)):
        da1,db1=np.shape(v)[0],np.shape(t)[0]
        factor=np.float(da1)/np.float(db1)
        newtemp=scipy.ndimage.zoom(t, factor, order=0)
    return newtemp,factor

# Segunda função chamada dentro da função "cloudfraction"
def getlatlon (file):
    fh = Dataset(file, mode='r')
    lon = fh.variables['lon'][:]
    lat = fh.variables['lat'][:]
    fh.close()
    return lat,lon

# Oitava função a ser chamada no loop principal
def reshapevis (A,B):
    A,B=np.squeeze(A),np.squeeze(B)
    newA=A
    if (np.shape(A)!=np.shape(B)):
        newA=rebinned(A,4)             
    return newA
    
# Primeira função chamada dentro da função "reshapevis"    
def rebinned ( X, factor ):
    X=np.squeeze(X)
    newx=int(np.trunc(np.shape(X)[0]/factor))
    newy=int(np.trunc(np.shape(X)[1]/factor))
    oldx=int(newx*factor)
    oldy=int(newy*factor)
    X=X[0:oldx,0:oldy]
    shape = [newx,newy]
    sh = shape[0],X.shape[0]//shape[0],shape[1],X.shape[1]//shape[1]
    r=X.reshape(sh).mean(-1).mean(1)
    return r    
    
# Nona função a ser chamada no loop principal
def phaseid (r063,r390,temp,name='//////my_ratio_rgb_image.png'):
    r063,r390,temp=np.squeeze(r063),np.squeeze(r390),np.squeeze(temp)
    fn=name.split('/')[len(name.split('/'))-1]
    outfn=fn[:-3]+'-phase.png'
    print 'Phase output:',outfn
    sat=fn.split('.')[0]
    year=fn.split('.')[1]
    julian=fn.split('.')[2]
    hhmm=name.split('.')[3]
    h=dt.datetime(int(year), 1, 1) + dt.timedelta(int(julian) - 1)    
    month=h.month    
    season=getseason(julian,year)
    icelim = icewet
    waterlim = waterwet
    if season == 'DRY':
        icelim = icedry
        waterlim = waterdry
    savereflectanceimg (r063,fn.split('BAND')[0]+'r063')
    savereflectanceimg (r390,fn.split('BAND')[0]+'r390')
    savereflectanceimg (temp,fn.split('BAND')[0]+'temp')
    ratio=r063/r390
    ice   = (ratio>icelim)     & (r063>cloudVISthreshold) & (temp<cloudTempthreshold) & (temp<zeroc)     #ice= (temp<=233) & (r063>=0.4) 
    water = (ratio<waterlim)   & (r063>cloudVISthreshold) & (temp<cloudTempthreshold) & (temp>=zeroc)     #water = (temp >= 273) & (r063 >= 0.4)
    mix   = (ratio >=waterlim) & (ratio<=icelim) & (r063>cloudVISthreshold) & (temp<cloudTempthreshold) & (temp<zeroc) & (temp>tlim)
    R,G,B=r063,r063,r063
    R,G,B=R-r063,R-r063,R-r063
    R[water]=r063[water]
    G[mix]=r063[mix]
    B[ice]=r063[ice]
    plt.close("all")
    f=plt.figure(facecolor='Black')
    ax = f.add_axes([0.1,0.1,0.9,0.9])
    img = np.zeros((r063.shape[0], r063.shape[1], 3), dtype=float)
    img[:,:,0],img[:,:,1],img[:,:,2]=r063,r063,r063
    img[:,:,0] = linear(R)  #, scale_min=0, scale_max=Rmax)     
    img[:,:,1] = linear(G)  #, scale_min=248, scale_max=298)  #np.max(G)/4)
    img[:,:,2] = linear(B)  #, scale_min=248, scale_max=298)   #248-298
    py.clf()
    py.imshow(img, aspect=1.8)    
    if month < 10:
        month='0'+str(month)
    label1=' '+str(sat)+'      '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    clear_frame()
    py.savefig(outfn,dpi=200.,bbox_inches='tight',pad_inches=0.0,transparent=True)
    plt.close("all")
    return

# Primeira função chamada dentro da função "phaseid"
def getseason (julian,year):
    season='WET'
    leap=0
    if (calendar.isleap(int(year))):
        leap=1
    dryini=dryseasonbegin+leap
    dryend=dryseasonend+leap
    if ((int(julian)>=dryini) and (int(julian)<=dryend)):
        season='DRY'
    return season

# Segunda função chamada dentro da função "phaseid"
def savereflectanceimg (data,file):
    procmode, vmin, vmax, myticks, mylabels = ' IR', 0.00, 1.00, [0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], ["0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"]
    bandnum=file.split('.')[len(file.split('.'))-1]
    data=np.squeeze(data)
    img = np.zeros((data.shape[0], data.shape[1], 2), dtype=float)      
    if bandnum == 'r063':
        procmode, vmin, vmax, myticks, mylabels = 'VIS', 0.00 , 1.00, [0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], ["0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"]
    if bandnum == 'temp':
        procmode, vmin, vmax, myticks, mylabels = 'TMP', 200.0, 315.0, [200,210,220, 230,240,250,260,270,280,290,300,310], ["200","210","220","230","240","250","260","270","280","290","300","310"]        
    img = linear(data, scale_min=vmin, scale_max=vmax)
    factor=vmax-vmin
    if factor != 0:
        img=img*factor + vmin
    f=plt.figure(facecolor='Black')
    f.set_size_inches(7,6,forward='True')
    ax = f.add_axes([0.02,0.1,0.9,0.9])
    py.clf()
    py.imshow(np.squeeze(img), aspect=1.8, vmin=vmin, vmax=vmax,cmap='spectral')
    cbar = plt.colorbar(orientation='horizontal', aspect=46.5,fraction=0.02,pad=0.0)
    cbar.set_ticks(myticks)
    cbar.set_ticklabels(mylabels)    
    if bandnum == 'temp':
        cbar.set_label('Temperature (K)', size=10)
    else:
        cbar.set_label('Reflectance', size=10)    
    fn=file+'.png'
    clear_frame()
    print "output: ", fn
    plt.savefig(fn, dpi=200.,bbox_inches='tight',pad_inches=0.0)
    plt.close("all")
    return

# Primeira função chamada dentro da função "savereflectanceimg", que por sua vez está dentro da "phaseid"
def linear(inputArray, scale_min=None, scale_max=None):
	"""Performs linear scaling of the input numpy array.

	@type inputArray: numpy array
	@param inputArray: image data array
	@type scale_min: float
	@param scale_min: minimum data value
	@type scale_max: float
	@param scale_max: maximum data value
	@rtype: numpy array
	@return: image data array
	
	"""			
	imageData=np.array(inputArray, copy=True)	
	if scale_min == None:
		scale_min = imageData.min()
	if scale_max == None:
		scale_max = imageData.max()
	imageData = imageData.clip(min=scale_min, max=scale_max)
	factor=(scale_max - scale_min)
     	if factor == 0:
		factor=1     
	imageData = (imageData -scale_min) / factor
	indices = np.where(imageData < 0)
	imageData[indices] = 0.0
	indices = np.where(imageData > 1)
	imageData[indices] = 1.0	
	return imageData

# Segunda função chamada dentro da função "savereflectanceimg", que por sua vez está dentro da "phaseid"
def clear_frame(ax=None): 
    if ax is None: 
        ax = plt.gca() 
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False) 
    for spine in ax.spines.itervalues(): 
        spine.set_visible(False) 

# Décima e última função a ser chamada no loop principal
def cloudprop (r063,r390,temp,mylat,mylon,myboxsize,SZA,VZA,name='my_ratio_rgb_image.png'):
    r063,r390,temp=np.squeeze(r063),np.squeeze(r390),np.squeeze(temp)
    SZA,VZA=np.squeeze(SZA),np.squeeze(VZA)
    lat, lon = getlatlon(name)
    dif=np.sqrt(((lat-mylat)*(lat-mylat))+((lon-mylon)*(lon-mylon)))
    loc=np.where(dif==np.min(dif))
    myx,myy=loc[0],loc[1]
    if (np.shape(myx)[0]>1):
        myx=myx[0]
    if (np.shape(myy)[0]>1):
        myy=myy[0]
    llat,llon=lat[myx,myy],lon[myx,myy]
    llat,llon=round(llat,2),round(llon,2)    
    if ((np.abs(llat-mylat)<1.0) and (np.abs(llon-mylon)<1.0) ):        
        print 'mylat,mylon',mylat,mylon
        print 'lat,lon',llat,llon
        os.chdir(outdir)
        myname = name.split('/')[len(name.split('/'))-1]    
        boxsize=[myboxsize] #[7,12,20]  Box size around the site, in pixels. Boxsize=7 means [x-7,y-7] to [x+7,y+7], or a square with 15 pixels of side.
        for b in boxsize:
            minx,miny,maxx,maxy=myx-b,myy-b,myx+b,myy+b
            r063b=r063[minx:maxx,miny:maxy]
            r390b=r390[minx:maxx,miny:maxy]
            tempb=temp[minx:maxx,miny:maxy]
            SZAb=SZA[minx:maxx,miny:maxy]
            VZAb=VZA[minx:maxx,miny:maxy]
            cloud=(r063b>cloudVISthreshold)&(tempb<cloudTempthreshold)
            fh = open(myname[:-11]+'_box'+str(b)+'_lat'+str(llat)+'_lon'+str(llon)+'_cloudreflectance063.csv', 'w')
            str1 = ''.join(str(e)+'\n' for e in r063b[cloud])
            fh.writelines(str1)
            fh.close()
            fh = open(myname[:-11]+'_box'+str(b)+'_lat'+str(llat)+'_lon'+str(llon)+'_cloudreflectance390.csv', 'w')
            str1 = ''.join(str(e)+'\n' for e in r390b[cloud])
            fh.writelines(str1)
            fh.close()
            fh = open(myname[:-11]+'_box'+str(b)+'_lat'+str(llat)+'_lon'+str(llon)+'_cloudtemperature.csv', 'w')
            str1 = ''.join(str(e)+'\n' for e in tempb[cloud])
            fh.writelines(str1)
            fh.close()
            fh = open(myname[:-11]+'_box'+str(b)+'_lat'+str(llat)+'_lon'+str(llon)+'_sza.csv', 'w')
            str1 = ''.join(str(e)+'\n' for e in SZAb[cloud])
            fh.writelines(str1)
            fh.close()
            fh = open(myname[:-11]+'_box'+str(b)+'_lat'+str(llat)+'_lon'+str(llon)+'_vza.csv', 'w')
            str1 = ''.join(str(e)+'\n' for e in VZAb[cloud])
            fh.writelines(str1)
            fh.close()            
    else:
        print 'Lat/Lon not found on the data:'
        print 'Required Lat/Lon: ',mylat,mylon
        print 'Closest match Lat/Lon: ',llat,llon
    return

  
###############  MAIN LOOP - LOOP PRINCIPAL ########################
warnings.simplefilter(action = "ignore", category = RuntimeWarning)

# Read input parameters
datadir   = sys.argv[1]
outdir    = sys.argv[2]
lesite    = sys.argv[3]
myboxsize = int(sys.argv[4])
plotphase = sys.argv[5]
plotphase = plotphase.split('=')[1]

# Define directories
curdir=os.getcwd()
os.chdir(outdir)

# Get list of files (Pega todos os caminhos dos arquivos .nc que pertencem a cada canal e guarda estes caminhos em listas, só que não abre os arquivos .nc ainda)
listch1=sorted(glob.glob(datadir+'/*.BAND_01.nc')) #'/home/andre/Desktop/Resultados_2/Dados/12-Dezembro/11UTC/2012/*.BAND_01.nc', por exemplo (canal 1, VIS)
listch2=sorted(glob.glob(datadir+'/*.BAND_02.nc'))
listch4=sorted(glob.glob(datadir+'/*.BAND_04.nc'))

# Get site location
[site,mylat,mylon,tz]=getsite(lesite) # Chamando a 1° função, a que pega a localização do site analisado. Retorna uma lista [site=RB,mylat=-09.956944,mylon=-67.868889,tz=-4], p.e.
print 'Processing site: ',site
print 'Lat/Lon: ',mylat,mylon
print 'Time zone: ',tz

# Main loop over the files
count=0
while count < len(listch1): #Enquanto o count for menor que o número de elementos na lista listch1, faço.... (no caso, para quatro dias em um UTC, terei 8 arquivos .nc para cada canal)
    [mVIS,kVIS,C,m2,m3,m4,m6,bb2,bb3,bb4,bb6,factor,F0,n2,n3,n4,n6,a2,a3,a4,a6,b2,b3,b4,b6,g2,g3,g4,g6,Bfunction,satlat,satlon]=getprelaunchconstant(listch1[count]) #Chamando a 2° função
    C=getpostlaunchconstant(listch1[count]) #Chamando a 3° função
    [r063,SZA,VZA]=reflectance(listch1[count]) #Chamando a 4° função somente para o canal VIS (o 1) e coloca [cs,SZA,VZA] em [r063,SZA,VZA]
    r063,SZA,VZA=np.squeeze(r063),np.squeeze(SZA),np.squeeze(VZA)    
    temp=np.squeeze(temperature(listch4[count])) #Chamando a 5° função
    [r390,SZAt,VZAt]=reflectance(listch2[count],temp,Bfunction) #Chamando a 4° função novamente, só que agora para o canal NIR (o 2) e coloca [cs,SZA,VZA] em [r390,SZA,VZA]
    r390=np.squeeze(r390)
    r063=checkvisshape(r063,temp) # after this, r063 will have either same shape as temp or exactly 16x bigger
    cloudfraction(r063,temp,mylat,mylon,myboxsize,listch1[count])
    r063=reshapevis(r063,r390) # this checks if r063 needs to be rebinned to match r390 shape
    if plotphase=='Yes':
        phaseid(r063,r390,temp,name=listch1[count])
    cloudprop(r063,r390,temp,mylat,mylon,myboxsize,SZA,VZA,name=listch2[count])
    count = count + 1
os.chdir(curdir)
exit()


