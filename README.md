# Cpu-Midas-Journal-800

The Most Beautiful Tool to Understand or Realize the registration on medical images
2D/3D Registration using itk framework

![volume](https://gitlab.com/chiconasa3000/Cpu-midas-journal-800/raw/master/Documentation/Images/approyvol.png "Volumen")

![volume](https://gitlab.com/chiconasa3000/Cpu-midas-journal-800/raw/master/Documentation/Images/approydrr.png "Volumen")

<img style="float: right;" src="https://gitlab.com/chiconasa3000/Cpu-midas-journal-800/raw/master/Documentation/Images/approyvol.png ">

## Test registration process by default using all specific commands

```
./multiTesting -numImag 5 -originVol ../inputData/pelvisSegmIntensityLPI.mha -targetVol ../inputData/pelvisSegmIntensityLPI.mha -compareVols -writeStatistics
```

## Specific Deformation or Relocated in volumetric template (Synthetic Images)
```
./CreateImageSetSimilarity  -v -folderName ../outputData/ImagesDefs -numImages 1 -rx 90.0 -ry 0.0 -rz 0.0 -t 0.0 0.0 0.0 -sg 1.0 -inputVol ../inputData/pelvisSegmIntensityLPI.mha -logFileName logRelocatedVolume.txt
```

## Specific Generation of Virtual Images
```
./genVirtualImage -v -p AP -dc 90 0 0 -foc 0 -1000 0 -scd 0 -res 1 1 -size 337 250 -threshold 0 -o pelvisHealthy_ap -inputVol ../outputData/ImagesDefs/Images/imagenDef_0.mha -logFileName logVirtualImage_ap.txt
```
```
./genVirtualImage -v -p ML -dc 90 90 0 -foc 1000 0 0 -scd 200 -res 1 1 -size 182 250 -threshold 0 -o pelvisHealthy_ml -inputVol ../outputData/ImagesDefs/Images/imagenDef_0.mha -logFileName logVirtualImage_ml.txt
```

## Specific Registration between 2 X-Ray Image and 3D Volume
```
./MultiImageRegistration -movingImage ../inputData/pelvisSegmIntensityLPI.mha -numFixedImages 2 -f0 ../outputData/virtualImages/pelvisHealthy_ap_0.mha 0 -1000 0 -f1 ../outputData/virtualImages/pelvisHealthy_ml_0.mha 1000 0 0 -steptolerance 0.02 -stepsize 4.0 -numLevels 3 -schedule 4 2 1 -outputDirectory ../outputData/resultsReg/ -logfilename logMultiRegistration.txt
```

Enjoy :wink:



