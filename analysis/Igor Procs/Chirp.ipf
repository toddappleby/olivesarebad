#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function Chirp()

//MLLoadWave /C/Y=4/S=2 "E:Data Analysis_2020:2020_0825:Chirp.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0122:Chirp.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0122:chirpStim.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0122:chirpFFT.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0122:chirpContrast.mat"

wave meanChirp
wave xvalsChirpSave
wave saveResult
wave yMag
wave logY
wave fMag
wave fPower
wave contrastPeaks
wave contrastResponsePeaks
NVAR maxFreq
string freqLabel = num2str(maxFreq)



Display saveResult vs xvalsChirpSave
ModifyGraph width=396,height=216
AppendToGraph/L=spikeRate meanChirp vs xvalsChirpSave
ModifyGraph rgb(meanChirp)=(0,0,0)
ModifyGraph axisEnab(left)={0,0.25}
ModifyGraph axRGB(left)=(65535,65535,65535),tlblRGB(left)=(65535,65535,65535)
ModifyGraph lblPosMode(left)=1,lblMargin(left)=35

Label left "\\f01\\Z09Stimulus\rContrast"

ModifyGraph axisEnab(spikeRate)={0.26,1}
ModifyGraph freePos(spikeRate)=0
ModifyGraph lblPosMode(spikeRate)=1,lblMargin(spikeRate)=18
Label spikeRate "\\f01Spiker Rate (Hz)"

Label bottom "\\f01Time (ms)"

string freqLabelFULL = "Max Freq:" + freqLabel + "Hz" 

SetDrawEnv xcoord= abs,ycoord= abs,textrgb= (16385,16388,65535)
DrawText 140,30,"Frequency Sweep"
SetDrawEnv xcoord= abs,ycoord= abs, fsize=9
DrawText 142,40, freqLabelFULL


SetDrawEnv xcoord= abs,ycoord= abs,textrgb= (36873,14755,58982)
DrawText 271,30,"Contrast Sweep"

//freq and contrast analysis below

AppendToGraph/L=contrast /B=sweep contrastResponsePeaks vs contrastPeaks
ModifyGraph freePos(sweep)=0
ModifyGraph axisEnab(bottom)={0,0.79},axisEnab(sweep)={0.81,1}
ModifyGraph tick(contrast)=1,nticks(contrast)=3,minor(contrast)=1,axisEnab(contrast)={0,0.40},freePos(contrast)={inf,sweep}
ModifyGraph tick(sweep)=1


AppendToGraph/L=frequency /B=sweep2 yMag vs fMag
ModifyGraph axisEnab(sweep2)={0.81,1}
ModifyGraph axisEnab(frequency)={0.6,1}
ModifyGraph nticks(frequency)=2,freePos(frequency)={inf,sweep}
ModifyGraph freePos(sweep2)={0,frequency}
SetAxis frequency 0,10
print maxFreq
print freqLabel
SetAxis sweep2 1,maxFreq
Label frequency "\\f01FFT Magnitude"
ModifyGraph rgb(yMag)=(16385,16388,65535)
ModifyGraph lblPos(frequency)=-5;
ModifyGraph tickEnab(frequency)={5,inf}
ModifyGraph tick(sweep2)=2,lblPos(sweep2)=35
Label sweep2 "\\f01Freq (Hz) (1 to max)";DelayUpdate

Label sweep "\\f01Relative Contrast";DelayUpdate
Label contrast "\\f01Peak Spike Rate (Hz)";DelayUpdate
ModifyGraph lblPos(contrast)=-6
ModifyGraph tick(sweep2)=1,lblPos(sweep2)=31 
ModifyGraph rgb(contrastResponsePeaks)=(36873,14755,58982)

ModifyGraph lblPos(sweep)=30



end
