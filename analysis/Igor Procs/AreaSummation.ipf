function addAreaSummation()

//MLLoadWave /C/Y=4/S=2 "E:Data Analysis_2020:2020_1118:expandingDarkSpots.mat"
MLLoadWave /C/Y=4/S=2 "Macintosh HD:Users:toddappleby:Documents:Data:Clarinet Exports:2021_0122:expandingDarkSpots.mat"
wave onRespFlip
wave offRespFlip
wave uniqueSpotSizesFlip

Display offRespFlip,onRespFlip vs uniqueSpotSizesFlip
ModifyGraph rgb(onRespFlip)=(0,0,65535)
ModifyGraph rgb(offRespFlip)=(65535,0,0)
Label bottom "Spot size diameter (µm)";DelayUpdate
Label left "Spike count"

Legend/C/N=text0/J "\\s(onRespFlip) On response\r\\s(offRespFlip) Off response"
Legend/C/N=text0/J/F=0
Legend/C/N=text0/J/X=-1.39/Y=50.71

end