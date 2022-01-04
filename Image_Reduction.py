from astropy.io import fits
import os
from os import system
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
#import datetime
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

def display(fname,path,X_Size,Y_Size):
    """
    fname: str- the name of the FITS file you want displayed
    returns: a linear greyscale image of the FITS file
    """
    imagefile=fits.open(path+fname)
    imagedata=imagefile[0].data
    x= np.arange(0,X_Size)
    y= np.arange(0,Y_Size)
    plt.pcolor(x,y,imagedata,cmap='Greys_r')
    plt.show()

def Average_Image(fname_list,path,X_Size,Y_Size):
    N=len(fname_list)
    #Master=np.empty([N,X_Size,Y_Size])
    Master=np.empty([N,Y_Size,X_Size])
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

def Image_Reduction(path,fname_key,Filter,Size=False,Date_Observation=False,CCD_Temp_Diff_Thresh=5):
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
        if((Obs_Filename_Test=='') or ("Avg_Bias" in Obs_Filename_Test) or ("Flat_Minus_Avg_Bias" in Obs_Filename_Test) or ("Minus_Dark" in Obs_Filename_Test) or ("Reduced" in Obs_Filename_Test)):
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
        Date_and_Time_Obs=hdul[0].header['DATE-OBS']
        print "Date_and_Time_Obs : ", Date_and_Time_Obs
        #Date_Obs=Date_and_Time_Obs.split("T")[0]
        #print "Date_Obs : ", Date_Obs
        #Day_Obs_Str=Date_Obs.split("-")[2]
        #Time_Obs=Date_Obs=Date_and_Time_Obs.split("T")[1]
        #print "Time_Obs : ", Time_Obs
        #Time_Obs_Hour_String=Time_Obs.split(":")[0]
        #print "Time_Obs_Hour_String : ", Time_Obs_Hour_String
        #Time_Obs_Hour=int(Time_Obs_Hour_String)
        #print "Time_Obs_Hour : ", Time_Obs_Hour

        #if(Date_Observation!=False):
        #if(isinstance(Date_Observation,basestring)):
            #Same_Night_Bool=(((Date_Obs==Date_Observation) and ((Time_Obs_Hour>17) or ()
        #Note: Need to create a boolean variable that is true if the observations ouccur on the same night
        Camera_Type=hdul[0].header['INSTRUME']
        print "Camera_Type : ", Camera_Type
        if(Camera_Type=="Apogee USB/Net"):
            Exposure_Time=hdul[0].header['EXPTIME']
            #print "Exposure_Time : ", Exposure_Time
            Set_Temp=hdul[0].header['SET-TEMP']
            #print "Set_Temp : ", Set_Temp
            CCD_Temp=hdul[0].header['CCD-TEMP']
            #print "CCD_Temp : ", CCD_Temp
            #X_Binning=hdul[0].header['XBINNING']
            #print "X_Binning : ", X_Binning
            #Y_Binning=hdul[0].header['YBINNING']
            #print "Y_Binning : ", Y_Binning
            Image_Type=hdul[0].header['IMAGETYP']
            print "Image_Type : ", Image_Type
            X_Axis_Size=hdul[0].header['NAXIS1'] #Not sure of NAXIS1 is the X-axis
            print "X_Axis_Size : ",  X_Axis_Size
            Y_Axis_Size=hdul[0].header['NAXIS2'] #Not sure of NAXIS2 is the Y-axis
            print "Y_Axis_Size : ",  Y_Axis_Size
            CCD_Temp_Diff=CCD_Temp-Set_Temp
            #print "CCD_Temp_Diff : ", CCD_Temp_Diff
            if((Image_Type=="Light Frame") or (Image_Type=="Flat Field")):
                Obs_Filter=hdul[0].header['FILTER']
                print "Obs_Filter : ", Obs_Filter
            """
            if(Size==False and (X_Axis_Size==Y_Axis_Size)):
                Size=X_Axis_Size
            if(X_Axis_Size!=Y_Axis_Size):
                print "Not a square CCD ! ! !"
                return False
                #return "Not a square CCD ! ! !"
            """
            X_Size=X_Axis_Size
            Y_Size=Y_Axis_Size
            #if((fname_key in Obs_Filename_Test or ("R" in Obs_Filename_Test)) and (int(X_Binning)==Bin) and (int(X_Binning)==Bin) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh)): #Conditions applicatble to all Image Types
            #if((fname_key in Obs_Filename_Test or ("R" in Obs_Filename_Test)) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh) and (Size==X_Axis_Size)): #Conditions applicatble to all Image Types
            #if((fname_key in Obs_Filename_Test) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh) and (Size==X_Axis_Size)): #Conditions applicatble to all Image Types
            #if((fname_key in Obs_Filename_Test) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh)): #Conditions applicatble to all Image Types
            if(CCD_Temp_Diff<CCD_Temp_Diff_Thresh): #Conditions applicatble to all Image Types
                if(Image_Type=="Bias Frame"): #Conditions for accepting a Bias
                    #print "BIAS TEST"
                    Bias_Fname_L.append(Obs_Filename_Test) #Conditions for accepting a Dark
                if(Image_Type=="Dark Frame"):
                    Dark_Fname_L.append(Obs_Filename_Test)
                if(Image_Type=="Flat Field"):
                    Min,Max,Avg=Max_Min_Ave(Obs_Filename_Test,path)
                    if((Avg>10000) and (Avg<30000)):
                        Flat_Fname_L.append(Obs_Filename_Test)
                if(fname_key==False):
                    if((Image_Type=="Light Frame") and (Obs_Filter==Filter)):
                        Light_Fname_L.append(Obs_Filename_Test)
                else:
                    if((Image_Type=="Light Frame") and (Obs_Filter==Filter) and (fname_key in Obs_Filename_Test)):
                        Light_Fname_L.append(Obs_Filename_Test)
        if(Camera_Type=="SBIG CCD"):
            Exposure_Time=hdul[0].header['EXPTIME']
            #print "Exposure_Time : ", Exposure_Time
            #Set_Temp=hdul[0].header['SET-TEMP'] #Does not work for 16" camera (SBIG CCD)
            #print "Set_Temp : ", Set_Temp
            CCD_Temp=hdul[0].header['CCD-TEMP']
            #print "CCD_Temp : ", CCD_Temp
            #X_Binning=hdul[0].header['XBINNING']
            #print "X_Binning : ", X_Binning
            #Y_Binning=hdul[0].header['YBINNING']
            #print "Y_Binning : ", Y_Binning
            #Image_Type=hdul[0].header['IMAGETYP'] #Does not work for 16" camera (SBIG CCD)
            Image_Type=hdul[0].header['FRAME']
            print "Image_Type : ", Image_Type
            #X_Axis_Size=hdul[0].header['NAXIS1'] #Not sure of NAXIS1 is the X-axis #Does not work for 16" camera (SBIG CCD)
            X_Axis_Size=3072
            print "X_Axis_Size : ",  X_Axis_Size
            #Y_Axis_Size=hdul[0].header['NAXIS2'] #Not sure of NAXIS2 is the Y-axis #Does not work for 16" camera (SBIG CCD)
            Y_Axis_Size=2048
            print "Y_Axis_Size : ",  Y_Axis_Size
            #CCD_Temp_Diff=CCD_Temp-Set_Temp #Does not work for 16" camera (SBIG CCD)
            #print "CCD_Temp_Diff : ", CCD_Temp_Diff
            #if((Image_Type=="Light Frame") or (Image_Type=="Flat Field")):
            if((Image_Type=="Light") or (Image_Type=="Flat")):
                Obs_Filter=hdul[0].header['FILTER']
                print "Obs_Filter : ", Obs_Filter
            """
            if(Size==False and (X_Axis_Size==Y_Axis_Size)):
                Size=X_Axis_Size
            if(X_Axis_Size!=Y_Axis_Size):
                print "Not a square CCD ! ! !"
                return False
            """
            X_Size=X_Axis_Size
            Y_Size=Y_Axis_Size
            #if((fname_key in Obs_Filename_Test or ("R" in Obs_Filename_Test)) and (int(X_Binning)==Bin) and (int(X_Binning)==Bin) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh)): #Conditions applicatble to all Image Types
            #if((fname_key in Obs_Filename_Test or ("R" in Obs_Filename_Test)) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh) and (Size==X_Axis_Size)): #Conditions applicatble to all Image Types
            #if((fname_key in Obs_Filename_Test) and (CCD_Temp_Diff<CCD_Temp_Diff_Thresh) and (Size==X_Axis_Size)): #Conditions applicatble to all Image Types
            #if((fname_key in Obs_Filename_Test) and (Size==X_Axis_Size)): #Conditions applicatble to all Image Types
            """
            if((fname_key in Obs_Filename_Test)): #Conditions applicatble to all Image Types
                if(Image_Type=="Bias"): #Conditions for accepting a Bias
                    #print "BIAS TEST"
                    Bias_Fname_L.append(Obs_Filename_Test) #Conditions for accepting a Dark
                if(Image_Type=="Dark"):
                    Dark_Fname_L.append(Obs_Filename_Test)
                if(Image_Type=="Flat"):
                    Flat_Fname_L.append(Obs_Filename_Test)
                if((Image_Type=="Light") and (Obs_Filter==Filter)):
                    Light_Fname_L.append(Obs_Filename_Test)
            """
            if(Image_Type=="Bias"): #Conditions for accepting a Bias
                #print "BIAS TEST"
                Bias_Fname_L.append(Obs_Filename_Test) #Conditions for accepting a Dark
            if(Image_Type=="Dark"):
                Dark_Fname_L.append(Obs_Filename_Test)
            if(Image_Type=="Flat"):
                Min,Max,Avg=Max_Min_Ave(Obs_Filename_Test,path)
                if((Avg>10000) and (Avg<30000)):
                    Flat_Fname_L.append(Obs_Filename_Test)
            if(fname_key==False):
                if((Image_Type=="Light") and (Obs_Filter==Filter)):
                    Light_Fname_L.append(Obs_Filename_Test)
            else:
                if((Image_Type=="Light") and (Obs_Filter==Filter) and (fname_key in Obs_Filename_Test)):
                    Light_Fname_L.append(Obs_Filename_Test)
    print "Bias_Fname_L : ", Bias_Fname_L
    Bias_Avg=Average_Image(Bias_Fname_L,path,X_Size,Y_Size)
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
    Flat_Minus_Avg_Bias_Avg=Average_Image(Flat_Minus_Avg_Bias_Fname_L,path,X_Size,Y_Size)
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
            Exposure_Time_Ratio=(float(Light_Exposure_Time)/float(Dark_Exposure_Time)) #The assumption of linear extropaltion of the Dark frames is not good enough by a long shot. Only Darks the same exposure time as the Lights should be used. The code should be modifed so that the only darks with the same exposure time as the current light are used to correct the current light
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
        #display(Reduced_Light_Fname,path,X_Axis_Size,Y_Axis_Size)






        #EXPTIME
        #hdul.close()

#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_71818/71818/","Image","V")
#72018_Modifed
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_72018/72018_Modifed_Image_Reduction_Code_Test/","Fireworks","V")
#72018_Modifed_Image_Reduction_Code_V2_Test
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_72018/72018_Modifed_Image_Reduction_Code_V2_Test/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_101818/Finding_Fireworks_GDrive/Fireworks_V_101818_Reduced_2/All_Reduced_Test/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_72018/72018_Modifed_Image_Reduction_Code_V3_Test/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_101818/Finding_Fireworks_GDrive/Fireworks_V_101818_Reduced_2/All_Reduced_Test/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_101818/Finding_Fireworks_GDrive/Fireworks_V_101818_Reduced_2/All_Reduced_Test_2/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_102418/Fireworks_V_102418_Reduction/All_Reduced/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_102418/Fireworks_V_102418_Reduction/All_Reduced/","Pleiades","V")
#Andromeda_Light
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_102418/Fireworks_V_102418_Reduction/All_Reduced/","Andromeda_Light","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_102918/FF_300s_Test/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_102918/Fireworks_Test_120s/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_102918/Fireworks_120s/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_V_103018/All_600s_Reduced/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Fireworks_Galaxy/Fireworks_Galaxy_103018/All_Reduced/","Fireworks","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_Reduced/","NGC_3631","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_2681/Galaxy_120718/All_Reduced/","NGC_2681","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_2681/Galaxy_12418/All_Reduced/","NGC_2681","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_020919/All_Reduced/","NGC_3631","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_2681/NGC_2681_121018/All_Reduced/","NGC_2681","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Orion/All_Reduced/","Orion","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Orion/All_Reduced/","Orion","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Orion/All_Reduced/","Orion","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_2681/Galaxy_031219/All_Reduced/","NGC_2681","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5813/Galaxy_031219/All_Reduced/","NGC_5813","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Andromeda","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Andromeda","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Andromeda","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Andromeda","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Orion","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Orion","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Orion","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Orion","H_alpha")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_3631/NGC_3631_013119/All_With_Joyrides_Reduced/","Orion","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/All_Reduced/","NGC_5474","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/Whirlpool_All_Reduced/","Pinwheel","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/Whirlpool_All_Reduced/","Pinwheel","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/Whirlpool_All_Reduced/","Pinwheel","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/M92_All_Reduced/","M92","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/M92_All_Reduced/","M92","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_060319/M92_All_Reduced/","M92","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M57/All_Reduced/","M57","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M57/All_Reduced/","M57","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M57/All_Reduced/","M57","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_5474/Galaxy_072519/All_Reduced/","NGC_5474","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_278/Galaxy_011519/All_Reduced/","NGC_278","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_278/Galaxy_072519/All_Reduced/","NGC_278","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_278/Galaxy_072719/All_Reduced/","NGC_278","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_278/Galaxy_082919/All_Reduced/","NGC_278","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Galaxies/NGC_278/Galaxy_091919/All_Reduced/","NGC_278","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M27/Testing_091619/All_Reduced/","M27","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M27/Testing_091619/All_Reduced/","M27","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M27/Testing_091619/All_Reduced/","M27","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M27/Testing_091619/All_Reduced/","M27","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Stars/Galaxy_092519/All_Reduced/","Stars","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Stars/Galaxy_092519/All_Reduced/","Stars","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Stars/Galaxy_092519/All_Reduced/","Stars","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Stars/Galaxy_092519/All_Reduced/","Stars","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Iris/Galaxy_092519/All_Reduced/","Iris","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Iris/Galaxy_092519/All_Reduced/","Iris","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Iris/Galaxy_092519/All_Reduced/","Iris","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/Iris/Galaxy_092519/All_Reduced/","Iris","R")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M1/Joyride_121919/All_Reduced/","M1","L")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M1/Joyride_121919/All_Reduced/","M1","B")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M1/Joyride_121919/All_Reduced/","M1","V")
#Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/Joyrides/M1/Joyride_121919/All_Reduced/","M1","R")
Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/New_24_Inch_Telescope/Orion_Test/All_Reduced/","M 42","L")
Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/New_24_Inch_Telescope/Orion_Test/All_Reduced/","M 42","B")
Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/New_24_Inch_Telescope/Orion_Test/All_Reduced/","M 42","V")
Image_Reduction("/Network/Servers/vimes.astro.wesleyan.edu/Volumes/vvodata/home/asantini/Desktop/24_Inch_Observations/New_24_Inch_Telescope/Orion_Test/All_Reduced/","M 42","R")
