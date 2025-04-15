global CerjanRate = 0.0053
global CerjanGridPoints = 45
CerjanBoundaryCondition(distance2) =  exp(-distance2*CerjanRate*CerjanRate)
#Study on Anisotropy of Shale in Shale Gas Exploration Using Seismic Forward Modeling
#March 2013
#DOI: 10.1002/cjg2.20017
#Weijia SunWeijia SunLi-Yun FuLi-Yun FuGUAN Xi‐ZhuWei WeiWei Wei

# Cerjan et al. suggest N= 20 and a = 0.0015 but Bording pointed out that it's better 
# to work with N = 45 and a = 0.0053