#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <New Polar Graphs Init>
#include <New Polar Graphs Draw>, version >= 7.05
#include <New Polar Keep On Screen>, version >= 6.13
#include <New Polar Graphs Cursors>
#include <New Polar Graphs>
#include <Graph Utility Procs>, version >=6.2
#include <SaveRestoreWindowCoords>, version >=7

//designed for oriented bar polar graph w/ example PSTH (or spike/raw) at some number of orientations (we'll see what looks good)


function makeOBarPolar()

WMPolarGraphGlobalsInit()

MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0128:oBarsResp.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0128:oBarsPSTH.mat"

//MLLoadWave /C/Y=4/S=2 "E:Data Analysis_2020:2020_0806:oBarsResp.mat"
//MLLoadWave /C/Y=4/S=2 "E:Data Analysis_2020:2020_0806:oBarsPSTH.mat"

wave offRespOut
wave onRespOut
string graphName = "oBarPolar"

//make graph
WMNewPolarGraph("whatisthis",graphName)

//put bar angles here
Make/O/N=12 numOBarAngles
numOBarAngles = {0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345}
WMPolarAppendTrace(graphName,offRespOut,numOBarAngles,360)
WMPolarAppendTrace(graphName,onRespOut,numOBarAngles,360)
//make plot folder directory so can change plot grid/fill
DFREF dfr = root:Packages:WMPolarGraphs:$graphName
SetDataFolder dfr
//turn off minor axes
Variable/G minorGridLineSize=0
//redraw
WMPolarAxesRedrawGraphNow(graphName)

ModifyGraph margin(left)=92,margin(bottom)=92,margin(right)=92,margin(top)=54
ModifyGraph width=256,height=256

//example subfig (labeled..)
SetDrawEnv xcoord= abs,ycoord= abs
DrawLine 75,241,75,295
SetDrawEnv xcoord= abs,ycoord= abs
DrawLine 33,295,75,295

SetDrawEnv xcoord= abs,ycoord= abs,textrot= 90,fsize=9;DelayUpdate
DrawText 77,295,"normalized spike rate";SetDrawEnv fsize= 9

SetDrawEnv xcoord= abs,ycoord= abs,fsize= 9;DelayUpdate
DrawText 47,306,"time"

SetDrawEnv xcoord= abs,ycoord= abs,arrow= 1;DelayUpdate
DrawLine 43,263,48,294;SetDrawEnv arrowfat= 0.30

SetDrawEnv xcoord= abs,ycoord= abs,arrow= 1;DelayUpdate
DrawLine 58,264,63,295;SetDrawEnv arrowfat= 0.30

SetDrawEnv xcoord= abs,ycoord= abs,textrot=-90,fsize=11;DelayUpdate
DrawText 31,262,"onset";

SetDrawEnv xcoord= abs,ycoord= abs,textrot=-90,fsize=11;DelayUpdate
DrawText 48,262,"offset";

SetDrawEnv xcoord= abs,ycoord= abs
DrawText 50,324,"0º"
SetDrawEnv xcoord= abs,ycoord= abs
DrawText 150,324,"45º"
SetDrawEnv xcoord= abs,ycoord= abs
DrawText 250,324,"135º"
SetDrawEnv xcoord= abs,ycoord= abs
DrawText 350,324,"165º"


end

//old shit 
ControlInfo/W=WMPolarGraphPanel mainModifyTracePop
NewDataFolder/O/S root:Packages
NewDataFolder/O/S WMPolarGraphs
NewDataFolder/O/S $graphName


function monkeyBusiness(polarity)
Variable polarity
String graphName= WMPolarTopPolarGraph()

if(polarity) // keeps consistent on/off colors in polar graph for on bars vs off bars 
print "hell yeah!"
Variable isFillToOriginOff = 1,isFillBehindOff=1,fillRedOff=52428,fillGreenOff=1,fillBlueOff=1,fillAlpha=5000
Variable isFillToOriginOn = 1,isFillBehindOn=1,fillRedOn=1,fillGreenOn=16019,fillBlueOn=65535,fillAlphOffa=5000

ModifyGraph rgb(polarY0)=(52428,1,1)
ModifyGraph rgb(polarY1)=(1,16019,65535)

TextBox/C/N=text1 "\\Z24Light Bar"
TextBox/C/N=text1/A=LT/X=-30/Y=-20
TextBox/C/N=text1/F=0

Legend/C/N=text0/J/A=MC "\\s(polarY1) On Resp\r\\s(polarY0) Off Resp"
Legend/C/N=text0/J/X=37/Y=-40/E=2

else

Legend/C/N=text0/J/A=MC "\\s(polarY1) On Resp\r\\s(polarY0) Off Resp"
Legend/C/N=text0/J/X=37/Y=40/E=2

TextBox/C/N=text1 "\\Z24Dark Bar"
TextBox/C/N=text1/A=LT/X=-30/Y=-20
TextBox/C/N=text1/F=0

 //TRYING TO FIX COLOR MISMATCH BETWEEN LEGEND AND PLOT
 isFillToOriginOn = 1 
 isFillBehindOn=1
 fillRedOn=52428
 fillGreenOn=1
 fillBlueOn=1
 fillAlpha=5000
 //isFillToOriginOff = 1 
 //isFillBehindOff=1
 //fillRedOff=52428
 //fillGreenOff=1
 //fillBlueOff=1
 //fillAlpha=5000
 //ontoff
 isFillToOriginOff = 1
 isFillBehindOff=1
 fillRedOff=1
 fillGreenOff=16019
 fillBlueOff=65535
 //fillAlphOna=5000  <<unnecessary diff between on and off
 
 ModifyGraph rgb(polarY1)=(52428,1,1)
 ModifyGraph rgb(polarY0)=(1,16019,65535)

endif


//off resp color (red)
	String polarTraceNameOff= "polarY0 [offRespOut]"
	//ModifyGraph rgb(polarY0)=(52428,1,1) lines now in if statement
	//Variable isFillToOriginOff = 1,isFillBehindOff=1,fillRedOff=52428,fillGreenOff=1,fillBlueOff=1,fillAlpha=5000
	String fillYWaveNameOff ="PolarFillY1",fillXWaveNameOff = "PolarFillX1"
	String df= WMPolarSetPolarTraceSettings(graphName,polarTraceNameOff,isFillToOriginOff,isFillBehindOff,fillRedOff,fillGreenOff,fillBlueOff,fillYWaveNameOff,fillXWaveNameOff,fillAlpha=fillAlpha)
	WMPolarModifyFillToOrigin("oBarPolar","polarY0")

//on resp color (blue)
String polarTraceNameOn= "polarY1 [onRespOut]"

	//ModifyGraph rgb(polarY1)=(1,16019,65535)
	//Variable isFillToOriginOn = 1,isFillBehindOn=1,fillRedOn=1,fillGreenOn=16019,fillBlueOn=65535,fillAlphOffa=5000
	String fillYWaveNameOn ="PolarFillY2",fillXWaveNameOn = "PolarFillX2"
	String dfOn= WMPolarSetPolarTraceSettings(graphName,polarTraceNameOn,isFillToOriginOn,isFillBehindOn,fillRedOn,fillGreenOn,fillBlueOn,fillYWaveNameOn,fillXWaveNameOn,fillAlpha=fillAlpha)
	WMPolarModifyFillToOrigin("oBarPolar","polarY1")



//TextBox/C/N=text1 "\\Z24Dark Bar"
//TextBox/C/N=text1/A=LT/X=-30/Y=-35
//TextBox/C/N=text1/F=0
end


function ohGod()

make/O/N=5000 numPoints = 0
make/O/N=5000 numPoints2 =1
make/O/N=5000 numPoints3 =2
Concatenate/NP/O{numPoints,numPoints2,numpoints3},numPointsFinal

Display /HOST=# ::::psthDataOut0 
RemoveFromGraph psthDataOut0
AppendToGraph/R ::::psthDataOut0
RenameWindow #,psth0
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.001, .78, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut0)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut0)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom " "
SetActiveSubwindow oBarPolar

Display /HOST=# ::::psthDataOut3
RemoveFromGraph psthDataOut3
AppendToGraph/R ::::psthDataOut3
RenameWindow #,psth3
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.24, .78, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut3)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut3)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom " "
SetActiveSubwindow oBarPolar

//thinking maybe 1 too many graphs makes whole fig look cluttered and confusing (when around polar plot at least...)
//Display /HOST=# ::::psthDataOut6
//RemoveFromGraph psthDataOut6
//AppendToGraph/R ::::psthDataOut6
//RenameWindow #,psth6
//SetAxis right *,1
//ModifyGraph width=65, height = 65
//MoveSubWindow fnum =(.39, .001, .1, .1)
//ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
//Label right "Normalized spike rate"
//ModifyGraph zColor(psthDataOut6)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut6)=(0,0,0)
//ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
//Label bottom "\\f0190"
//SetActiveSubwindow oBarPolar


Display /HOST=# ::::psthDataOut9
RemoveFromGraph psthDataOut9
AppendToGraph/R ::::psthDataOut9
RenameWindow #,psth9
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.48, .78, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut9)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut9)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom " "
SetActiveSubwindow oBarPolar

Display /HOST=# ::::psthDataOut11
RemoveFromGraph psthDataOut11
AppendToGraph/R ::::psthDataOut11
RenameWindow #,psth11
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.7, .78, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut11)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut11)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom " "
SetActiveSubwindow oBarPolar


end