

.\bin\ycolorgrade --image tests\greg_zaal_artist_workshop.hdr --exposure 0 --output out\greg_zaal_artist_workshop_01.jpg
.\bin\ycolorgrade --image tests\greg_zaal_artist_workshop.hdr --exposure 1 --filmic --contrast 0.75 --saturation 0.75 --output out\greg_zaal_artist_workshop_02.jpg
.\bin\ycolorgrade --image tests\greg_zaal_artist_workshop.hdr --exposure 0.8 --contrast 0.6 --saturation 0.5 --grain 0.5 --output out\greg_zaal_artist_workshop_03.jpg

.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --exposure -1 --filmic --contrast 0.75 --saturation 0.3 --vignette 0.4 --output out\toa_heftiba_people_01.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --exposure -0.5 --contrast 0.75 --saturation 0 --output out\toa_heftiba_people_02.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --exposure -0.5 --contrast 0.6 --saturation 0.7 --tint-red 0.995 --tint-green 0.946 --tint-blue 0.829 --grain 0.3 --output out\toa_heftiba_people_03.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --mosaic 16 --grid 16 --output out\toa_heftiba_people_04.jpg


.\bin\ycolorgrade --srgb --image tests\toa_heftiba_people.jpg --trueAnaglyph --output out\ExtraCredit_.trueAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\toa_heftiba_people.jpg --grayAnaglyph --output out\ExtraCredit_grayAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\toa_heftiba_people.jpg --colorAnaglyph --output out\ExtraCredit_colorAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\toa_heftiba_people.jpg --halfColorAnaglyph --output out\ExtraCredit_halfColorAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\toa_heftiba_people.jpg --optimizedColorAnaglyph --output out\ExtraCredit_optimizedColorAnaglyph.jpg






.\bin\ycolorgrade --image tests\fotoPerSoften.jpg --srgb  --soften --LivelloSoften 1 --output out\ExtraCredit_soften1.jpg

.\bin\ycolorgrade --image tests\fotoPerSoften.jpg --srgb  --soften --LivelloSoften 2 --output out\ExtraCredit_soften2.jpg
.\bin\ycolorgrade --image tests\fotoPerSoften.jpg --srgb  --soften --LivelloSoften 0 --output out\ExtraCredit_soften0.jpg



.\bin\ycolorgrade --image tests\forest.jpg --srgb  --oilPainting --output out\ExtraCredit_oilPaintingForest20.jpg
.\bin\ycolorgrade --image tests\forest.jpg --srgb  --oilPainting --intensityOil 5 --output out\ExtraCredit_oilPaintingForest5.jpg
.\bin\ycolorgrade --image tests\greg_zaal_artist_workshop.hdr --srgb --oilPainting --output out\ExtraCredit_oilPaintinggreg_zaal_artist_workshop_03.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb --oilPainting --output out\ExtraCredit_oilPaintingPeople.jpg
.\bin\ycolorgrade --image tests\colosseo.jpg --srgb --oilPainting --output out\ExtraCredit_oilPaintingColosseo20.jpg
.\bin\ycolorgrade --image tests\colosseo.jpg --srgb --oilPainting --intensityOil 40 --output out\ExtraCredit_oilPaintingColosseo40.jpg






.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --inverti --output out\ExtraCredit_inversione.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --blur --blurIntensity 3 --output out\ExtraCredit_blur3.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --blur --blurIntensity 10 --output out\ExtraCredit_blur10.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --blur --blurIntensity 20 --output out\ExtraCredit_blur20.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --cut --output out\ExtraCredit_cut.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --edgeDetection --livelloEdgeDetection 1 --output out\ExtraCredit_edgeDetection1.jpg
.\bin\ycolorgrade --image tests\colosseo.jpg --srgb  --edgeDetection --livelloEdgeDetection 1 --output out\ExtraCredit_edgeDetectionColosseo.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --edgeDetection --livelloEdgeDetection 2 --output out\ExtraCredit_edgeDetection2.jpg
.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb  --edgeDetection --livelloEdgeDetection 0 --output out\ExtraCredit_edgeDetection0.jpg





.\bin\ycolorgrade --srgb --image tests\colosseo.jpg --emboss --output out\ExtraCredit_emboss.jpg
.\bin\ycolorgrade --srgb --image tests\greg_zaal_artist_workshop.hdr --enhanceDetails --output out\ExtraCredit_enhanceDetails.jpg
.\bin\ycolorgrade --srgb --image tests\greg_zaal_artist_workshop.hdr --gaussianBlur --output out\ExtraCredit_gaussianBlur.jpg



.\bin\ycolorgrade --srgb --image tests\greg_zaal_artist_workshop.hdr --trueAnaglyph --output out\ExtraCredit_.trueAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\greg_zaal_artist_workshop.hdr --grayAnaglyph --output out\ExtraCredit_grayAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\greg_zaal_artist_workshop.hdr --colorAnaglyph --output out\ExtraCredit_colorAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\greg_zaal_artist_workshop.hdr --halfColorAnaglyph --output out\ExtraCredit_halfColorAnaglyph.jpg
.\bin\ycolorgrade --srgb --image tests\forest.jpg --optimizedColorAnaglyph --output out\ExtraCredit_optimizedColorAnaglyph.jpg



.\bin\ycolorgrade --image tests\fotoPerSoften.jpg --srgb --sketch --output out\ExtraCredit_sketchPersona.jpg
.\bin\ycolorgrade --image tests\colosseo.jpg --srgb --sketch --output out\ExtraCredit_sketchColosseo.jpg


.\bin\ycolorgrade --image tests\toa_heftiba_people.jpg --srgb --CRTTV --output out\ExtraCredit_CRTTVpeople.jpg
.\bin\ycolorgrade --image tests\colosseo.jpg --srgb --CRTTV --output out\ExtraCredit_CRTTVColosseo.jpg






