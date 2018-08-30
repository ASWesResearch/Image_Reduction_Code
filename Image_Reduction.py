from astropy.io import fits
import os
from os import system
import numpy as np
import matplotlib.pyplot as plt
#print "TESTING"

def Max_Min_Ave(fname,path):
    """
    fname: str- The name of the FITS file that will have it's Max, Min and Average pixel value displayed
    """
    imagefile=fits.open(path+fname)
    imagedata=imagefile[0].data
    Max=np.max(imagedata)
    Min=np.min(imagedata)
    Ave=np.mean(imagedata)
    print "The Max pixel value is " + str(Max) +" ADU"
    print "The Min pixel value is " + str(Min) +" ADU"
    print "The Average pixel vale is " + str(Ave)+" ADU"
    return Min,Max,Ave

def display(fname,path):
    """
    fname: str- the name of the FITS file you want displayed
    returns: a linear greyscale image of the FITS file
    """
    imagefile=fits.open(path+fname)
    imagedata=imagefile[0].data
    x= np.arange(0,2048)
    y= np.arange(0,2048)
    plt.pcolor(x,y,imagedata,cmap='Greys_r')
    plt.show()

def Average_Image(fname_list,path):
    N=len(fname_list)
    Master=np.empty([N,2048,2048])
    for i in range(0,N):
        fname=fname_list[i]
        imagefile=fits.open(path+fname)
        imagedata=imagefile[0].data
        Master[i]=imagedata
    Max=np.amax(Master,0)
    Sum=np.sum(Master,0)
    Sub=Sum-Max
    Ave=Sub/(N-1)
    return Ave

def Write(A,new_fname,path):
    new_image=A
    new_HDU=fits.PrimaryHDU(new_image)
    new_HDU.writeto(path+new_fname)

def Difference(fname1,fname2,path):
    imagefile1=fits.open(path+fname1)
    imagedata1=imagefile1[0].data
    imagefile2=fits.open(path+fname2)
    imagedata2=imagefile2[0].data
    image_Diff=(imagedata1-imagedata2)
    return image_Diff

def Divide(fname1,fname2,path):
    imagefile1=fits.open(path+fname1)
    imagedata1=imagefile1[0].data
    imagefile2=fits.open(path+fname2)
    imagedata2=imagefile2[0].data
    Ratio=imagedata1/imagedata2
    return Ratio

def Image_Reduction(path,fname_key,Filter,Expo_Time=60,Bin=1,CCD_Temp_Diff_Thresh=5):
    #print "Hello_World"
    Obs_File_LS=os.popen("ls "+path).read()
    #print "Obs_File_LS : ", Obs_File_LS
    Obs_File_LS_L=Obs_File_LS.split("\n")
    print "Obs_File_LS_L : ",Obs_File_LS_L
    Bias_Fname_L=[]
    Dark_Fname_L=[]
    Flat_Fname_L=[]
    Light_Fname_L=[]
    for Obs_Filename_Test in Obs_File_LS_L:
        if((Obs_Filename_Test=='') or ("Avg_Bias" in Obs_Filename_Test) or("Flat_Minus_Avg_Bias" in Obs_Filename_Test) or ("Minus_Dark" in Obs_Filename_Test) or ("Reduced" in Obs_Filename_Test)):
            continue
        print "Obs_Filename_Test : ", Obs_Filename_Test
        #print "Obs_Filename_Test : ", Obs_Filename_Test
        #Obs_Filename_Test_L=Obs_Filename_Test.split("_")
        Obs_Filename_Test_L=Obs_Filename_Test.split(" ")
        #print "Obs_Filename_Test_L : ",Obs_Filename_Test_L
        for Obs_Filename_Segment in Obs_Filename_Test_L:
            if((".fit" in Obs_Filename_Segment) or  (".fits" in Obs_Filename_Segment)):
                Obs_Filename_Test_Corrected=Obs_Filename_Segment
        #print "Obs_Filename_Test_Corrected : ", Obs_Filename_Test_Corrected
        hdul = fits.open(path + Obs_Filename_Test)
        #print hdul
        #hdul.info()
        Date_Obs=hdul[0].header['DATE-OBS']
        #print "Date_Obs : ", Date_Obs
        Exposure_Time=hdul[0].header['EXPTIME']
        #print "Exposure_Time : ", Exposure_Time
        Set_Temp=hdul[0].header['SET-TEMP']
        #print "Set_Temp : ", Set_Temp
        CCD_Temp=hdul[0].header['CCD-TEMP']
        #print "CCD_Temp : ", CCD_Temp
        X_Binning=hdul[0].header['XBINNING']
        #print "X_Binning : ", X_Binning
        Y_Binning=hdul[0].header['YBINNING']
        #print "Y_Binning : ", Y_Binning
        Image_Type=hdul[0].header['IMAGETYP']
        print "Image_Type : ", Image_Type
        CCD_Temp_Diff=CCD_Temp-Set_Temp
        #print "CCD_Temp_Diff : ", CCD_Temp_Diff
        if((Image_Type=="Light Frame") or (Image_Type=="Flat Field")):
            Obs_Filter=hdul[0].header['FILTER']
            print "Obs_Filter : ", Obs_Filter
        if((fname_key in Obs_Filename_Test or ("R" in Obs_Filename_Test)) and (int(X_Binning)==Bin) and (int(X_Binning)==Bin) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh)): #Conditions applicatble to all Image Types
            if(Image_Type=="Bias Frame"): #Conditions for accepting a Bias
                #print "BIAS TEST"
                Bias_Fname_L.append(Obs_Filename_Test) #Conditions for accepting a Dark
            if(Image_Type=="Dark Frame"):
                Dark_Fname_L.append(Obs_Filename_Test)
            if(Image_Type=="Flat Field"):
                Flat_Fname_L.append(Obs_Filename_Test)
            if((Image_Type=="Light Frame") and (Obs_Filter==Filter)):
                Light_Fname_L.append(Obs_Filename_Test)
    print "Bias_Fname_L : ", Bias_Fname_L
    Bias_Avg=Average_Image(Bias_Fname_L,path)
    Bias_Avg_fname=fname_key+"_Avg_Bias.fit"
    #print "Bias_Avg_fname : ", Bias_Avg_fname
    if(Bias_Avg_fname not in Obs_File_LS_L):
        Write(Bias_Avg,Bias_Avg_fname,path)
    #display(Bias_Avg_fname,path)
    Flat_Minus_Avg_Bias_Fname_L=[]
    for Flat_Fname in Flat_Fname_L:
        Flat_Minus_Avg_Bias=Difference(Flat_Fname,Bias_Avg_fname,path)
        Flat_Fname_Fit_L=Flat_Fname.split(".")
        #print "Flat_Fname_Fit_L : ", Flat_Fname_Fit_L
        Flat_Fname_Fit_L_Not_Fits=Flat_Fname_Fit_L[0]
        Flat_Minus_Avg_Bias_Fname=Flat_Fname_Fit_L_Not_Fits+"_Minus_Avg_Bias.fit"
        Flat_Minus_Avg_Bias_Fname_L.append(Flat_Minus_Avg_Bias_Fname)
        if(Flat_Minus_Avg_Bias_Fname not in Obs_File_LS_L):
            Write(Flat_Minus_Avg_Bias,Flat_Minus_Avg_Bias_Fname,path)
    Flat_Minus_Avg_Bias_Avg=Average_Image(Flat_Minus_Avg_Bias_Fname_L,path)
    Flat_Minus_Avg_Bias_Fname_L_Fits=Flat_Minus_Avg_Bias_Fname.split(".")
    #print "Flat_Minus_Avg_Bias_Fname_L_Fits : ", Flat_Minus_Avg_Bias_Fname_L_Fits
    Flat_Minus_Avg_Bias_Avg_Fname=Flat_Minus_Avg_Bias_Fname_L_Fits[0]+"_Avg.fit"
    print "Flat_Minus_Avg_Bias_Avg_Fname : ", Flat_Minus_Avg_Bias_Avg_Fname
    if(Flat_Minus_Avg_Bias_Avg_Fname not in Obs_File_LS_L):
        Write(Flat_Minus_Avg_Bias_Avg,Flat_Minus_Avg_Bias_Avg_Fname,path)
    #display(Flat_Minus_Avg_Bias_Avg_Fname,path)
    Flat_Minus_Avg_Bias_Avg_Min,Flat_Minus_Avg_Bias_Avg_Max,Flat_Minus_Avg_Bias_Avg_Avg=Max_Min_Ave(Flat_Minus_Avg_Bias_Avg_Fname,path)
    print "Flat_Minus_Avg_Bias_Avg_Min : ", Flat_Minus_Avg_Bias_Avg_Min
    print "Flat_Minus_Avg_Bias_Avg_Max : ", Flat_Minus_Avg_Bias_Avg_Max
    print "Flat_Minus_Avg_Bias_Avg_Avg : ", Flat_Minus_Avg_Bias_Avg_Avg
    Flat_Minus_Avg_Bias_Avg_Normalized=Flat_Minus_Avg_Bias_Avg/Flat_Minus_Avg_Bias_Avg_Avg
    Flat_Minus_Avg_Bias_Avg_Fname_L_Fits=Flat_Minus_Avg_Bias_Avg_Fname.split(".")
    Flat_Minus_Avg_Bias_Avg_Fname_Not_Fit=Flat_Minus_Avg_Bias_Avg_Fname_L_Fits[0]
    Flat_Minus_Avg_Bias_Avg_Normalized_Fname=Flat_Minus_Avg_Bias_Avg_Fname_Not_Fit+"_Normalized.fit"
    print "Flat_Minus_Avg_Bias_Avg_Normalized_Fname : ", Flat_Minus_Avg_Bias_Avg_Normalized_Fname
    if(Flat_Minus_Avg_Bias_Avg_Normalized_Fname not in Obs_File_LS_L):
        Write(Flat_Minus_Avg_Bias_Avg_Normalized,Flat_Minus_Avg_Bias_Avg_Normalized_Fname,path)
    #display(Flat_Minus_Avg_Bias_Avg_Normalized_Fname,path)
    for Light_Fname in Light_Fname_L:
        print "Light_Fname : ", Light_Fname
        Light_Fname_Fit_L=Light_Fname.split(".")
        Light_Fname_No_Fit=Light_Fname_Fit_L[0]
        Light_hdul=fits.open(path + Light_Fname)
        Light_Exposure_Time=Light_hdul[0].header['EXPTIME']
        Light_Data=Light_hdul[0].data
        #print "Light_Exposure_Time : ", Light_Exposure_Time
        #print "type(Light_Exposure_Time) : ", type(Light_Exposure_Time)
        Dark_Adj_L=[]
        for Dark_Fname in Dark_Fname_L:
            Dark_hdul=fits.open(path + Dark_Fname)
            Dark_Exposure_Time=Dark_hdul[0].header['EXPTIME']
            #print "Dark_Exposure_Time : ", Dark_Exposure_Time
            #imagefile1=fits.open(path+fname1)
            Dark_Data=Dark_hdul[0].data
            #print "Dark_Data : ", Dark_Data
            Exposure_Time_Ratio=(float(Light_Exposure_Time)/float(Dark_Exposure_Time))
            #print "Exposure_Time_Ratio : ", Exposure_Time_Ratio
            Adjusted_Dark_Data=Exposure_Time_Ratio*Dark_Data
            #print "Adjusted_Dark_Data : ", Adjusted_Dark_Data
            Dark_Adj_L.append(Adjusted_Dark_Data)
        #print "Dark_Adj_L : ", Dark_Adj_L
        Dark_Adj_Avg=np.mean(Dark_Adj_L, axis=0) #Not quite sure if this is correct
        Dark_Adj_Avg_Shape=Dark_Adj_Avg.shape #It has the correct shape though
        #print "Dark_Adj_Avg_Shape : ", Dark_Adj_Avg_Shape
        #print "Dark_Adj_Avg : \n", Dark_Adj_Avg
        #Dark_Adj_Avg_2=np.mean(Dark_Adj_L, axis=1)
        #print "Dark_Adj_Avg : \n", Dark_Adj_Avg_2
        #Dark_Adj_Avg_Shape_2=Dark_Adj_Avg_2.shape
        #print "Dark_Adj_Avg_Shape_2 : ", Dark_Adj_Avg_Shape_2
        Light_Minus_Dark_Data=(Light_Data-Dark_Adj_Avg)
        Light_Minus_Dark_Fname=Light_Fname_No_Fit+"_Minus_Dark.fit"
        print "Light_Minus_Dark_Fname : ", Light_Minus_Dark_Fname
        if(Light_Minus_Dark_Fname not in Obs_File_LS_L):
            Write(Light_Minus_Dark_Data,Light_Minus_Dark_Fname,path)
        Reduced_Light_Data=Divide(Light_Minus_Dark_Fname,Flat_Minus_Avg_Bias_Avg_Normalized_Fname,path)
        Reduced_Light_Fname=Light_Fname_No_Fit+"_Reduced.fit"
        print "Reduced_Light_Fname : ", Reduced_Light_Fname
        if(Reduced_Light_Fname not in Obs_File_LS_L):
            Write(Reduced_Light_Data,Reduced_Light_Fname,path)
        #display(Reduced_Light_Fname,path)






        #EXPTIME
        #hdul.close()

#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_71818/71818/","Image","V")
#72018_Modifed
Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_72018/72018_Modifed_Image_Reduction_Code_Test/","Fireworks","V")
