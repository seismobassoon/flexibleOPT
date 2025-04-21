using UnPack

# Here we should be able to even read ak135 etc. for 1D Earth or photos

function readOrMakeDiscretisedModel(vars; option="toyModel")
    materialModels = []
    if option === "toyModel"
        for var in vars
            TXYZDependency = findTXYZDependency(var)
            newstring = "離散化" * split(string(var), "(")[1]
            perturbationType = "1Dto3DcocentricSinePerturbation"
            ParamSizeTXYZ = [ModelSizeTXYZ.ModelSizeT * TXYZDependency[1], ModelSizeTXYZ.ModelSizeX * TXYZDependency[2], ModelSizeTXYZ.ModelSizeY * TXYZDependency[3], ModelSizeTXYZ.ModelSizeZ * TXYZDependency[4]]

            tmpModel = makeHeterogeneousModels(perturbationType, ParamSizeTXYZ, Δnodes)
            materialModels = push!(materialModels, tmpModel)
        end

        return materialModels
    end
end

function readOrMakeDiscretisedModel(vars, filenames; option="2Dphotos")

end

#

function regularGridConstruction(DomainWindow, ModelSizeTXYZ)

    # one method to make grids in model domain
    # NF need to generalise this

    @unpack DomainWindowT, DomainWindowX, DomainWindowY, DomainWindowZ = DomainWindow
    @unpack ModelSizeT, ModelSizeX, ModelSizeY, ModelSizeZ = ModelSizeTXYZ

    # full window in 4D in [s] and [m] or dimensionless TXYZ
    DomainDimension = [DomainWindowT DomainWindowX DomainWindowY DomainWindowZ]

    numberNodes = [ModelSizeT ModelSizeX ModelSizeY ModelSizeZ]

    Δnode = zeros(4)
    for iComp in 1:4
        if numberNodes[iComp] != 0
            Δnode[iComp] = DomainDimension[iComp] / (numberNodes[iComp] - 1)
        end
    end

    Δnodes = (ΔnodeT=Δnode[1], ΔnodeX=Δnode[2], ΔnodeY=Δnode[3], ΔnodeZ=Δnode[4])
    return Δnodes
end


function makeHeterogeneousModels(type, ParamSizeTXYZ, Δnodes; averageValue=1.0, perturbationValue=0.2, wavelength=0.3)
    # wavelength is in m or dimensionless
    #region unpacking, etc.

    @unpack ΔnodeT, ΔnodeX, ΔnodeY, ΔnodeZ = Δnodes

    ParamSizeT = ParamSizeTXYZ[1]
    ParamSizeX = ParamSizeTXYZ[2]
    ParamSizeY = ParamSizeTXYZ[3]
    ParamSizeZ = ParamSizeTXYZ[4]

    # I need to detect the TXYZ dependency in order to construct a model
    spaceDimension = 3
    localDimension = ()

    if ParamSizeT === 0
        ParamSizeT = 1
    end
    if ParamSizeX === 0
        ParamSizeX = 1
        spaceDimension -= 1
    else
        localDimension = (localDimension..., ParamSizeX)
    end

    if ParamSizeY === 0
        ParamSizeY = 1
        spaceDimension -= 1
    else
        localDimension = (localDimension..., ParamSizeY)
    end

    if ParamSizeZ === 0
        ParamSizeZ = 1
        spaceDimension -= 1
    else
        localDimension = (localDimension..., ParamSizeZ)
    end

    centreT = trunc(Int, ParamSizeT / 2)
    centreX = trunc(Int, ParamSizeX / 2)
    centreY = trunc(Int, ParamSizeY / 2)
    centreZ = trunc(Int, ParamSizeZ / 2)

    #endregion

    #region model construction

    modelSpace = Array{Any,spaceDimension}(undef, localDimension)

    if type === "1Dto3DisocentricSinePerturbation"

        tmpArray = Array{Any,3}(undef, ParamSizeX, ParamSizeY, ParamSizeZ)
        for iz in 1:ParamSizeZ
            for iy in 1:ParamSizeY
                for ix in 1:ParamSizeX
                    x = ΔnodeX * (ix - 1 - centreX)
                    y = ΔnodeY * (iy - 1 - centreY)
                    z = ΔnodeZ * (iz - 1 - centreZ)
                    distance = sqrt(x^2 + y^2 + z^2)
                    tmpArray[ix, iy, iz] = averageValue + perturbationValue * sin(distance / wavelength * 2 * π)
                end
            end
        end

        modelSpace = reshape(tmpArray, localDimension)
    end
    #endregion

    return modelSpace

end
