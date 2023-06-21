#!/bin/sh

# Make plots for the Ferrara models.
mkdir -p plots/Ferrara1999Original
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_disk.pdf            --inclinationIndex 0 --spheroidScaleIndex -1 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_spheroidR0p116.pdf  --inclinationIndex 0 --spheroidScaleIndex  0 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_spheroidR0p464.pdf  --inclinationIndex 0 --spheroidScaleIndex  1 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_spheroidR1p856.pdf  --inclinationIndex 0 --spheroidScaleIndex  2 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_disk.pdf           --inclinationIndex 8 --spheroidScaleIndex -1 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_spheroidR0p116.pdf --inclinationIndex 8 --spheroidScaleIndex  0 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_spheroidR0p464.pdf --inclinationIndex 8 --spheroidScaleIndex  1 --showAtlas 1 --showExtrapolation 0
./AVVsOpticalDepth.pl --attenuationFile data/Ferrara1999_MW_hz1.0:Attenuations.hdf5 --plotFile plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_spheroidR1p856.pdf --inclinationIndex 8 --spheroidScaleIndex  2 --showAtlas 1 --showExtrapolation 0
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_disk.pdf            plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc0_disk.pdf           
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_spheroidR0p116.pdf  plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc0_spheroidR0.116.pdf 
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_spheroidR0p464.pdf  plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc0_spheroidR0.464.pdf 
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc0_spheroidR1p856.pdf  plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc0_spheroidR1.856.pdf 
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_disk.pdf           plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc90_disk.pdf           
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_spheroidR0p116.pdf plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc90_spheroidR0.116.pdf 
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_spheroidR0p464.pdf plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc90_spheroidR0.464.pdf 
mv plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1p0_inc90_spheroidR1p856.pdf plots/Ferrara1999Original/AVVsOpticalDepth_MW_hz1.0_inc90_spheroidR1.856.pdf 

exit 0
