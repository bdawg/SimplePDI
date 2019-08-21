loadFilename = 'allImsCube_MWC480_02-04_20190225_broadband_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAur_01-02_20180108_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_domeflat_750_em1000_20ms_20190320_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAur_01_20190320_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_domeflat_750_em1000_20ms_20190320_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAur_01_20181017_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_HD98800_01_20190320_Open_EmptySlot_0.npz'

loadFilename = 'allImsCube_ABAur_01_20190320_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAur_01_20181017_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAur_02-03_combined_20181022_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAUr_01_20190226_750-50_EmptySlot_0.npz'

loadFilename = 'allImsCube_50pc_HD154445_02_20180626_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_HD154445_01_20190226_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_WDS_02_20190226_750-50_EmptySlot_0.npz'


dataPath = '../SimplePDI_DATA/'


from simplePDI import * 
p = pdiImages(outDir='../SimplePDI_DATA/')
p.loadCube(dataPath, loadFilename)     
p.makeDiffIms(showPlots=False)
# # p.makeDiffIms(showPlots=False, Ibias=200)
p.plotStokesIms(crop=0.5, apertureRad=20)

# Ibias = 200
# p.makeDiffIms(showPlots=False, Ibias=Ibias, deRotate=True, rotatePolz=True)
# p.makeDiffIms(showPlots=False, Ibias=Ibias, deRotate=False, rotatePolz=True)
# p.makeDiffIms(showPlots=False, Ibias=Ibias, deRotate=True, rotatePolz=False)
# p.makeDiffIms(showPlots=False, Ibias=Ibias, deRotate=False, rotatePolz=False)
# p.plotStokesIms(crop=0.6)


plt.ion()