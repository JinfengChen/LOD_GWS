echo "cread new position file for LOD score"
python LOD_position.py --input ../input/MPR.cross.uniq.QTL.mr.table

echo "Draw LOD score along with GWAS signals"
python LOD_GWAS.py --input QTL.regions.GWAS.list

echo "specific script for one QTL"
python LOD_GWAS_HeadingDate1.py --input QTL.regions.GWAS.HeadingDate1.list
python LOD_GWAS_HeadingDate2.py --input QTL.regions.GWAS.HeadingDate2.list
python LOD_GWAS_HeadingDate3.py --input QTL.regions.GWAS.HeadingDate3.list
