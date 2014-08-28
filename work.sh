echo "cread new position file for LOD score"
python LOD_position.py --input ../input/MPR.cross.uniq.QTL.mr.table

echo "Draw LOD score along with GWAS signals"
python LOD_GWAS.py --input QTL.regions.GWAS.list

