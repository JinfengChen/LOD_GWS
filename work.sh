echo "cread new position file for LOD score"
python LOD_position.py --input ../input/MPR.cross.uniq.QTL.mr.table

echo "Draw LOD score along with GWAS signals"
python LOD_GWAS.py --input QTL.regions.GWAS.list

echo "specific script for one QTL"
python LOD_GWAS_HeadingDate1.py --input QTL.regions.GWAS.HeadingDate1.list
python LOD_GWAS_HeadingDate2.py --input QTL.regions.GWAS.HeadingDate2.list
python LOD_GWAS_HeadingDate3.py --input QTL.regions.GWAS.HeadingDate3.list

python LOD_GWAS_PlantHeight1.py --input QTL.regions.GWAS.PlantHeight1.list
python LOD_GWAS_PlantHeight2.py --input QTL.regions.GWAS.PlantHeight2.list
python LOD_GWAS_PlantHeight3.py --input QTL.regions.GWAS.PlantHeight3.list

python LOD_GWAS_GrainYield1.py --input QTL.regions.GWAS.GrainYield1.list
python LOD_GWAS_GrainYield2.py --input QTL.regions.GWAS.GrainYield2.list
python LOD_GWAS_GrainYield3.py --input QTL.regions.GWAS.GrainYield3.list
python LOD_GWAS_GrainYield4.py --input QTL.regions.GWAS.GrainYield4.list

python LOD_GWAS_Auxin.py --input QTL.regions.GWAS.Auxin.list
