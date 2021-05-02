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


function mbCreate()

WMPolarGraphGlobalsInit()

MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0128:mBarsResp.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0128:mBarPSTH.mat"

//MLLoadWave /C/Y=4/S=2 "E:Data Analysis_2020:2020_0806:mBarsResp.mat"
//MLLoadWave /C/Y=4/S=2 "E:Data Analysis_2020:2020_0806:mBarsPSTH.mat"

wave RespOut
string graphName = "mBarPolar"

//make graph
WMNewPolarGraph("whatisthis2",graphName)

//put bar angles here
Make/O/N=12 numOBarAngles
numOBarAngles = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330}
WMPolarAppendTrace(graphName,respOut,numOBarAngles,360)
//make plot folder directory so can change plot grid/fill
DFREF dfr = root:Packages:WMPolarGraphs:$graphName
SetDataFolder dfr
//turn off minor axes
Variable/G minorGridLineSize=0
//redraw
WMPolarAxesRedrawGraphNow(graphName)

ModifyGraph margin(left)=92,margin(bottom)=92,margin(right)=92,margin(top)=97
ModifyGraph width=261,height=256
end

//old shit 
ControlInfo/W=WMPolarGraphPanel mainModifyTracePop
NewDataFolder/O/S root:Packages
NewDataFolder/O/S WMPolarGraphs
NewDataFolder/O/S $graphName


function mbAnnotate(polarity)

Variable polarity
String graphName= WMPolarTopPolarGraph()

if(polarity) // keeps consistent on/off colors in polar graph for on bars vs off bars 
print "hell yeah!"

Variable isFillToOrigin = 1,isFillBehind=1,fillRed=1,fillGreen=16019,fillBlue=65535,fillAlpha=5000

ModifyGraph rgb(polarY0)=(1,16019,65535)

TextBox/C/N=text1 "\\Z18Light Bar"
TextBox/C/N=text1/A=LT/X=-30/Y=-35
TextBox/C/N=text1/F=0

else

TextBox/C/N=text1 "\\Z18Dark Bar"
TextBox/C/N=text1/A=LT/X=-30/Y=-35
TextBox/C/N=text1/F=0

 isFillToOrigin= 1 
 isFillBehind=1
 fillRed=1
 fillGreen=16019
 fillBlue=65535
 fillAlpha=5000
 
 ModifyGraph rgb(polarY0)=(1,16019,65535)


endif


//off resp color (red)
	String polarTraceName= "polarY0 [respOut]"
	//ModifyGraph rgb(polarY0)=(52428,1,1) lines now in if statement
	//Variable isFillToOriginOff = 1,isFillBehindOff=1,fillRedOff=52428,fillGreenOff=1,fillBlueOff=1,fillAlpha=5000
	String fillYWaveName ="PolarFillY0",fillXWaveName = "PolarFillX0"
	String df= WMPolarSetPolarTraceSettings(graphName,polarTraceName,isFillToOrigin,isFillBehind,fillRed,fillGreen,fillBlue,fillYWaveName,fillXWaveName,fillAlpha=fillAlpha)
	WMPolarModifyFillToOrigin("mBarPolar","polarY0")



//TextBox/C/N=text1 "\\Z24Dark Bar"
//TextBox/C/N=text1/A=LT/X=-30/Y=-35
//TextBox/C/N=text1/F=0
end


function mbSubs(regionOut)

WAVE regionOut

if (regionOut[0] ==1 && numpnts(regionOut) == numpnts(::::psthDataOut0))
print "ALL"
make/O/N=(numpnts(regionOut)) numPointsFinal = 1

//text box for what measured
TextBox/C/N=text2 "\\Z18 Whole Trace"
TextBox/C/N=text2/A=LT/X=-32/Y=-28
TextBox/C/N=text2/F=0

elseif (regionOut[0] ==1 && numpnts(regionOut) < numpnts(::::psthDataOut0))
print "onset"

//text box for what measured
TextBox/C/N=text2 "\\Z18 Onset Response"
TextBox/C/N=text2/A=LT/X=-32/Y=-28
TextBox/C/N=text2/F=0

make/O/N=(numpnts(regionOut)) numPoints1 = 1
make/O/N=(numpnts(::::psthDataOut0)-regionOut(numpnts(regionOut)-1)) numPoints2=0

Concatenate/NP/O{numPoints1,numPoints2},numPointsFinal

elseif (regionOut[0] != 1)

//text box for measure
TextBox/C/N=text3 "\\Z18 Offset Response"
TextBox/C/N=text3/A=LT/X=-32/Y=-28
TextBox/C/N=text3/F=0

print "offset or some spicy custom range"

make/O/N=(regionOut[0]) numPoints1 = 0
make/O/N=(regionOut[numpnts(regionOut)]) numPoints2 = 1
make/O/N=(numpnts(::::psthDataOut0)-regionOut[numpnts(regionOut)]) numPoints3=0

Concatenate/NP/O{numPoints1,numPoints2,numPoints3},numPointsFinal
endif



//make/O/N=5000 numPoints = 0
//make/O/N=5000 numPoints3 =2
//Concatenate/NP/O{numPoints,numPoints2,numpoints3},numPointsFinal

Display /HOST=# ::::psthDataOut0 
RemoveFromGraph psthDataOut0
AppendToGraph/R ::::psthDataOut0
RenameWindow #,psth0
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.73, .31, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut0)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut0)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 0"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut2
RemoveFromGraph psthDataOut2
AppendToGraph/R ::::psthDataOut2
RenameWindow #,psth2
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.61, .08, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut2)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut2)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 60"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut3
RemoveFromGraph psthDataOut3
AppendToGraph/R ::::psthDataOut3
RenameWindow #,psth3
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.37, .001, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
Label right "Normalized spike rate"
ModifyGraph zColor(psthDataOut3)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut3)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 90"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut4
RemoveFromGraph psthDataOut4
AppendToGraph/R ::::psthDataOut4
RenameWindow #,psth4
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.12, .08, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut4)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut4)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 120"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut6
RemoveFromGraph psthDataOut6
AppendToGraph/R ::::psthDataOut6
RenameWindow #,psth6
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.001, .315, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut6)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut6)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 180"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut8
RemoveFromGraph psthDataOut8
AppendToGraph/R ::::psthDataOut8
RenameWindow #,psth8
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.12, .72, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut8)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut8)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 240"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut9
RemoveFromGraph psthDataOut9
AppendToGraph/R ::::psthDataOut9
RenameWindow #,psth9
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.39, .78, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut9)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut9)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 270"
SetActiveSubwindow mBarPolar

Display /HOST=# ::::psthDataOut11
RemoveFromGraph psthDataOut11
AppendToGraph/R ::::psthDataOut11
RenameWindow #,psth11
SetAxis right *,1
ModifyGraph width=65, height = 65
MoveSubWindow fnum =(.69, .6, .1, .1)
ModifyGraph axRGB(bottom)=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535),alblRGB(bottom)=(65535,65535,65535)
ModifyGraph zColor(psthDataOut11)={numPointsFinal,1,2,BlueBlackRed,0},zColorMin(psthDataOut11)=(0,0,0)
ModifyGraph axOffset(bottom)=-2,standoff(bottom)=0,alblRGB=(0,0,0)
Label bottom "\\f01Angle: 330"
SetActiveSubwindow mBarPolar

end