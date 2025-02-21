import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
import json
from collections import Counter
import natsort
from scipy.constants import e
GHIN = False

if GHIN:
    sys.path.append(r"C:\Users\vidar\Documents")
    import GitHub.DMU.DMU.utils as dm
    import GitHub.DMU.DMU.plot_utils as dmp
    import GitHub.DMU.DMU.graph_styles as gss
else:
    from DMU import utils as dm
    from DMU import plot_utils as dmp
    from DMU import graph_styles as gss
import os 
##setting a default style:
gssstyle = "PP2_4by3"
bboxstyle = "wide symmetric"
gss.graph_style(gssstyle)
DEFbbox = gss.DEF_BBOX(style=gssstyle, bboxstyle=bboxstyle)
#%%        
#First step is to initialise your DataImportSettings
#Load user settings from a file using CUV
UV = dm.CUV(act = 'init')
if "Processed" not in UV.keys():
    UV["Processed"] = False
plotfile = "MieData.json"
#If you have a same named .txt file with single named variables in a tab delimited format, you can set the following line to true

#If you do not want to use the default .json file, then you can run Init_LDI() on its own in the console, and then you can select a new json file using CUV() -> "ddir" command from the option select..

#Then we want to add folder with your data, and this data can be assorted or not. You do this using DataDir() and then selecting add using 0 in the console

#Once you have added at least one data directory in to your json file, you can now choose a dataset to load.

"""
DEFINING DATA DIRECTORIES AND GETTING THE FILES FROM THE CORRECT ONE
"""
#We get a list of all files in dicts matching the number of extensions we are searching for.
# UV["Processed"] = False # Change once your data is processed properly.

if UV["Processed"] == False:
    DATA = {}
    plotdata = {}
    for key,item in dm.DataDir(act="load").items():
        DIR,DIRPT = [os.path.abspath(item),"abs"]
        
        FLTuple = DList,NList = dm.Get_FileList(DIR,pathtype=DIRPT, ext = (('mat','json')),sorting='numeric')
        plotdata["Filename"]           = []
        plotdata["Dx"]                  = []
        plotdata["Dy"]                  = []
        plotdata["Dz"]                  = []
        plotdata["CS"]                  = {}
        plotdata["CS"]["Qabs"]          = []
        plotdata["CS"]["Qscat"]         = []
        plotdata["CS"]["Qabs_theory"]   = []
        plotdata["CS"]["Qscat_theory"]  = []
        plotdata["lambda"]              = []
        plotdata["FE"]                  = {}
        plotdata["FE"]["lambda_max"]    = []
        plotdata["FE"]["Ere"]            = []
        plotdata["FE"]["Eim"]            = []
        plotdata["FE"]["x"]            = []
        plotdata["FE"]["y"]            = []
        plotdata["FE"]["z"]            = []
        
        for i,matfile in enumerate(FLTuple[0][".mat"]):

            DATA =  dm.MatLoader(FLTuple[0]['.mat'][i],json=True)[0]
            plotdata["Filename"].append(FLTuple[1]['.mat'][i])
            FE = DATA["MieData"]["FieldEnhancement"][0]
            CS = DATA["MieData"]["Csects"][0]
            
            plotdata["Dx"].append(DATA["Dx"])
            plotdata["Dy"].append(DATA["Dy"])
            plotdata["Dz"].append(DATA["Dz"])
            plotdata["CS"]["Qabs"].append(CS["Abssim"])
            plotdata["CS"]["Qabs_theory"].append(CS["Abstheory"])
            plotdata["CS"]["Qscat"].append(CS["Scatsim"])
            plotdata["CS"]["Qscat_theory"].append(CS["Scattheory"])
            
            fpoint                           = int(FE["fpoint"][0]) 
            plotdata["lambda"]               = FE["E"][0]["lambda"][0]
            #plotdata["FE"]["lambda_max"].append(plotdata["lambda"][fpoint])
            #plotdata["FE"]["Ere"].append(np.real(FE["E"][0]["E"]))
            #plotdata["FE"]["Eim"].append(np.imag(FE["E"][0]["E"])) #If you import E field, file will be HUGE - instead, we will retroactively load the E-fields with the largest and smallest scattering cross section.
            #plotdata["FE"]["x"].append(FE["E"][0]["x"])
            #plotdata["FE"]["y"].append(FE["E"][0]["y"])
            #plotdata["FE"]["z"].append(FE["E"][0]["z"])
        
    
                        
        #%%
           
        dm.json_savedata(plotdata,"RAW_"+key+"_"+plotfile,overwrite=True)        
        UV["Processed"] = True
    #%%
FIG = []   
fid = 0 
      
all_data = []
for file in ["RAW_1_MieData.json","RAW_3_MieData.json","RAW_2_MieData.json"]:
    all_data.append(dm.json_loaddata(file))    
PD = []
MAXvals = {"ABS":[],"SCAT":[],"TOT":[],"ALL":[]}

for i,pltdata in enumerate(all_data):  
    PD.append( {"XX":[],"YY":[],"ZZ":[],"SCAT":[],"ABS":[],"NAMES":{}}    )
    
    x = np.unique(pltdata["Dx"]) 
    y = np.unique(pltdata["Dy"]) 
    z = np.unique(pltdata["Dz"])
    rx = len(x); ry = len(y); rz = len(z);

    ind = 0
    SCAT   = np.ones([len(x),len(y),len(z)])
    ABS    = np.ones([len(x),len(y),len(z)])
    COORDS = np.meshgrid(x,y,z)
    #SCAT = np.ones([26,26,5])
    SCAT   = np.array([np.sum(scatval) for scatval in pltdata["CS"]["Qscat"]])
    ABS    = np.array([np.sum(absval) for absval in pltdata["CS"]["Qabs"]])
    NAMES  = np.array([filename for filename in pltdata["Filename"]])  
    
    NAMEMATRIX = np.reshape(NAMES, [rx,ry,rz]).transpose(1, 0, 2)
    SCAT = np.reshape(SCAT, [rx,ry,rz]).transpose(1, 0, 2)
    ABS = np.reshape(ABS, [rx,ry,rz]).transpose(1, 0, 2)
    XX = np.reshape(np.array(pltdata["Dx"]),[rx,ry,rz]).transpose(1, 0, 2)
    YY = np.reshape(np.array(pltdata["Dy"]),[rx,ry,rz]).transpose(1, 0, 2)
    ZZ = np.reshape(np.array(pltdata["Dz"]),[rx,ry,rz]).transpose(1, 0, 2)
    
    PD[i]["SCAT"] = SCAT; PD[i]["ABS"] = ABS; PD[i]["NAMES"] = NAMEMATRIX
    PD[i]["XX"] = XX; PD[i]["YY"] = YY; PD[i]["ZZ"] = ZZ
    
    MAXvals["SCAT"].append([np.min(SCAT),np.max(SCAT)])
    MAXvals["ABS"].append([np.min(ABS),np.max(ABS)])
    MAXvals["TOT"].append([np.min([np.min(SCAT),np.min(ABS)]),np.max([np.max(SCAT),np.max(ABS)])])

cmap = plt.get_cmap("plasma")    
vmin = np.min(MAXvals["TOT"])
vmax=np.max(MAXvals["TOT"])
cnorm = mpl.colors.Normalize(vmin=vmin, 
                             vmax=vmax) 

def find_nearest(array, value):
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return(idx)

#intensity 2d plots, and intensity line plots.
FFIG  = dm.ezplot()
#LFIG  = [dm.ezplot(gspec=[1,2],projection="3d") for i in PD]
LFIG  = [dm.ezplot(projection="3d") for i in PD]
LLFIG = dm.ezplot(projection="3d") 
zind = 1
zinda = [zind,find_nearest(np.unique(PD[1]["ZZ"]),np.unique(PD[0]["ZZ"])[zind]),find_nearest(np.unique(PD[2]["ZZ"]),np.unique(PD[0]["ZZ"])[zind])] 
 
for i,pltdata in enumerate(PD):
    FIG.append(dm.ezplot())
    XXp,YYp,SCATp = [pltdata["XX"][:,:,zind],pltdata["YY"][:,:,zind],pltdata["SCAT"][:,:,zind]]
    zind = zinda[i]
    FIG[fid].ax[0].pcolormesh(XXp,YYp,SCATp, cmap=cmap,norm=cnorm)
    FIG[fid].apply_bbox(bbox=DEFbbox)
    FFIG.ax[0].pcolormesh(XXp,YYp,SCATp, cmap=cmap,norm=cnorm)
    FFIG.ax[0].set_xlabel("Director Dx [m]")
    FFIG.ax[0].set_ylabel("Director Dy [m]")
    FFIG.ax[0].annotate("Director Dz = " + "{:.3e}".format(np.mean(pltdata["ZZ"][:,:,zind])),xy = (0.3,0.9-0.05*fid),xycoords="figure fraction")
    fid+= 1

cmap = plt.get_cmap("cool")    
vmin = np.min([np.min(PDZZ["ZZ"]) for PDZZ in PD])
vmax= np.max([np.max(PDZZ["ZZ"]) for PDZZ in PD])
cnorm = mpl.colors.Normalize(vmin=vmin, 
                             vmax=vmax) 
fid = 0

for i,pltdata in enumerate(PD):
    for zind in range(pltdata["XX"].shape[2]):
        for yi in range(pltdata["YY"][:,:,zind].shape[1]):
            XXp,YYp,SCATp = [pltdata["XX"][:,yi,zind],pltdata["YY"][:,yi,zind],pltdata["SCAT"][:,yi,zind]]
            LLFIG.ax[0].plot(XXp,YYp,SCATp,linewidth=1,c=cmap(cnorm(np.max(pltdata["ZZ"][:,:,zind]))) )
            LLFIG.ax[0].view_init(elev=60, azim=10, roll=0)
            LLFIG.ax[0].set_zlim([0,400])
            LLFIG.ax[0].set_xlabel("Director Dx [m]",fontsize=15)
            LLFIG.ax[0].set_ylabel("Director Dy [m]",fontsize=15)
            LLFIG.ax[0].tick_params(labelsize=10)
            LLFIG.ax[0].yaxis.offsetText.set_fontsize(15)
            LLFIG.ax[0].xaxis.offsetText.set_fontsize(15)
            
            LFIG[fid].ax[0].plot(XXp,YYp,SCATp+1e+9*np.max(pltdata["ZZ"][:,:,zind]),c=cmap(cnorm(np.max(pltdata["ZZ"][:,:,zind]))) )
            LFIG[fid].ax[0].view_init(elev=60, azim=10, roll=0)
            LFIG[fid].ax[0].set_zlim([0,400])
            LFIG[fid].ax[0].set_xlabel("Director Dx [m]",fontsize=15)
            LFIG[fid].ax[0].set_ylabel("Director Dy [m]",fontsize=15)
            LFIG[fid].ax[0].yaxis.offsetText.set_fontsize(15)
            LFIG[fid].ax[0].xaxis.offsetText.set_fontsize(15)
            #More purple, is higher z value

    fid+= 1
#%%
#Determining where the max value 
maxarrays = [np.max(PD[i]["SCAT"]) for i,val in enumerate(PD) ]
pdI = np.where(maxarrays == np.max(maxarrays))[0][0]
xyzi = np.array([val[0] for i,val in enumerate(np.where(PD[pdI]["SCAT"] == np.max(PD[pdI]["SCAT"])))])

XX,YY,ZZ = [PD[pdI]["XX"][*xyzi],PD[pdI]["YY"][*xyzi],PD[pdI]["ZZ"][*xyzi]]

print("Ideal dimensions for the director are as follows:")
print("Dx = " + "{:.3f}".format(XX*1e+9) + " nm")
print("Dy = " + "{:.3f}".format(YY*1e+9) + " nm")
print("Dz = " + "{:.3f}".format(ZZ*1e+9) + " nm")

print("\nWhich represents:")
print("Dx =  lambda/" + "{:.3f}".format(1/(XX*1e+9/910)))
print("Dy = lambda/" + "{:.3f}".format(1/(YY*1e+9/910)))
print("Dz = lambda/" + "{:.3f}".format(1/(ZZ*1e+9/910)))
    
#%%
savefig = False
figformat = [".png",".pdf"]
if savefig == True:
    for fig in FIG:
        for ext in figformat:
            fig.fig.savefig("Figures\\"+fig.file+ext,transparent=True)

