# import ROOT
import numpy as np
# from ROOT import *
import matplotlib.pyplot as plt
import os


import csv

cmsPathRoot = "Inputs/HEPData-CMS.root"


def readCMSRootFile(filepath):
    f_ = TFile.Open(filepath, "OPEN")

    # get xSec values absolute
    foldername = "parton_abs_rho"
    histname = "Hist1D_y1"
    totErr = "Hist1D_y1_e1"
    # statErr = "Hist1D_y1_e2"


def readCMSCSVFiles(filepath):
    # open the xSec table
    name = "parton_abs_rho"

    xSec_abs_binEdges = []
    xSec_abs_val = []
    xSec_abs_err_tot_up = []
    xSec_abs_err_tot_down = []
    # xSec_abs_err_stat_up = []
    # xSec_abs_err_stat_down = []
    xSec_abs_center = []

    with open(filepath+name+".csv", mode ='r')as file:
        csvFile = csv.reader(file)
        for iline,lines in enumerate(csvFile):
            if iline > 8:
                if len(lines)>0:
                    # print(iline, lines)
                    center = float(lines[0])
                    lowEdge = float(lines[1])
                    upEdge = float(lines[2])
                    val = float(lines[3])
                    tot_up = float(lines[4])
                    tot_down = float(lines[5])
                    # stat_up = float(lines[6])
                    # stat_down = float(lines[7])

                    xSec_abs_center.append(center)
                    xSec_abs_val.append(val)
                    xSec_abs_err_tot_up.append(tot_up)
                    xSec_abs_err_tot_down.append(tot_down)
                    # xSec_abs_err_stat_up.append(stat_up)
                    # xSec_abs_err_stat_down.append(stat_down)
                    
                    if not lowEdge in xSec_abs_binEdges:
                        xSec_abs_binEdges.append(lowEdge)
                    if not upEdge in xSec_abs_binEdges:
                        xSec_abs_binEdges.append(upEdge)

    # print (xSec_abs_binEdges)
    # print (xSec_abs_val)
    # print (xSec_abs_err_stat_up)
    # print (xSec_abs_err_stat_down)

    # open the extrapolation uncertainties cov matrix
    name = "parton_abs_rho_covariance_extrapolation"
    covMatrix_extr = [[[] for iVal in xSec_abs_val] for iVal in xSec_abs_val]
    with open(filepath+name+".csv", mode ='r')as file:
        csvFile = csv.reader(file)
        for iline,lines in enumerate(csvFile):
            if iline > 8:
                if len(lines)>0:
                    # print(iline, lines)
                    binEdge_1 = int(float(lines[0]))
                    binEdge_2 = int(float(lines[3]))
                    val = float(lines[-1])
                    # print binEdge_1,binEdge_2, val
                    covMatrix_extr[binEdge_1][binEdge_2] = val

    # open the total uncertainties cov matrix
    name = "parton_abs_rho_covariance_total"
    corrMatrix_total = [[[] for iVal in xSec_abs_val] for iVal in xSec_abs_val]
    with open(filepath+name+".csv", mode ='r')as file:
        csvFile = csv.reader(file)
        for iline,lines in enumerate(csvFile):
            if iline > 8:
                if len(lines)>0:
                    # print(iline, lines)
                    binEdge_1 = int(float(lines[0]))
                    binEdge_2 = int(float(lines[3]))
                    val = float(lines[-1])
                    # print binEdge_1,binEdge_2, val
                    corrMatrix_total[binEdge_1][binEdge_2] = val


    # open the nuisance corr table
    # first get the list of nuisances
                        
    csvContent = []

    name = "fit_np_correlation"
    allNuisancesNames = []
    allNuisancesIdx = []
    with open(filepath+name+".csv", mode ='r') as file:
        csvFile = csv.reader(file)
        for iline,lines in enumerate(csvFile):
            csvContent.append(lines)
            if "Nuisance parameter name" in lines:
                # print(iline, lines)
                allNuisancesNames.append(lines[1])
                allNuisancesIdx.append(iline)

    # print (allNuisancesName)
    # print (allNuisancesIdx)
    corrMatrix = [[[]for iNuisance in allNuisancesNames] for iNuisance in allNuisancesNames]
    for i,idxNuisance_i in enumerate(allNuisancesIdx):
        for j, idxNuisance_j in enumerate(allNuisancesIdx):
            # print (i,j)
            corrLines = csvContent[idxNuisance_i+1+j]
            # corrLines = csvContent[idxNuisance_i:idxNuisance_i+j]
            # print (allNuisancesNames[i], allNuisancesNames[j])
            # print (corrLines)
            corr = float(corrLines[1])
            # print (corr)
            corrMatrix[i][j] = corr
            # corrMatrix

    # print (csvContent)

    # print (corrMatrix)
    # print (len(corrMatrix))
    # print (len(corrMatrix[-1]))



    # now get the rest of all nuisances
    csvContent2 = []
    name = "np_impacts_pulls"
    fit_obs_central = []
    fit_obs_up = []
    fit_obs_down = []
    fit_exp_central = []
    fit_exp_up = []
    fit_exp_down = []
    impact_r1_obs_up = []
    impact_r1_obs_down = []
    impact_r2_obs_up = []
    impact_r2_obs_down = []
    impact_r3_obs_up = []
    impact_r3_obs_down = []
    impact_r4_obs_up = []
    impact_r4_obs_down = []
    # impact_r1_exp_up = []
    # impact_r1_exp_down = []
    # impact_r2_exp_up = []
    # impact_r2_exp_down = []
    # impact_r3_exp_up = []
    # impact_r3_exp_down = []
    # impact_r4_exp_up = []
    # impact_r4_exp_down = []

    ignores = ["rate_ttj0", "rate_ttj1", "rate_ttj2", "rate_ttj3"]
    IdxChanges = []
    with open(filepath+name+".csv", mode ='r') as file:
        csvFile = csv.reader(file)
        for iline,lines in enumerate(csvFile):
            csvContent2.append(lines)
            if "Nuisance parameter name" in lines:
                IdxChanges.append(iline)

    # print (csvContent2)

    # print (csvContent2[IdxChanges[0]+1:IdxChanges[0]+1+len(allNuisancesNames)-4])

    for idxNuisance, nuisanceName in enumerate(allNuisancesNames):
        contentToConsider = csvContent2[IdxChanges[0]+1:IdxChanges[0]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                # print (idxNuisance, nuisanceName,float(entry[1]))
                fit_obs_central.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[1]+1:IdxChanges[1]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                fit_obs_up.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[2]+1:IdxChanges[2]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                fit_obs_down.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[6]+1:IdxChanges[6]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r1_obs_up.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[7]+1:IdxChanges[7]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r1_obs_down.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[8]+1:IdxChanges[8]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r2_obs_up.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[9]+1:IdxChanges[9]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r2_obs_down.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[10]+1:IdxChanges[10]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r3_obs_up.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[11]+1:IdxChanges[11]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r3_obs_down.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[12]+1:IdxChanges[12]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r4_obs_up.append(float(entry[1]))

        contentToConsider = csvContent2[IdxChanges[13]+1:IdxChanges[13]+1+len(allNuisancesNames)-4]
        for entry in contentToConsider:
            if nuisanceName in entry:
                impact_r4_obs_down.append(float(entry[1]))

    # print (fit_obs_central)
    # print (impact_r4_obs_down)

    outFolder = "CMS_ttj_inputfiles"
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    outFileName = "rho_cms_inputs.txt"
    f = open(outFolder+"/"+outFileName, "w")

    f.write("[hessian]\n")
    f.write("\n")
    f.write("[end hessian]\n")
    f.write("\n")


    # now write this all out in convino style.... :/
    # first the estimates (xsec values)
    f.write("[estimates]\n")
    f.write("\n")
    # print 'n_estimates = '+str(len(xSec_abs_val))
    f.write('n_estimates = '+str(len(xSec_abs_val))+"\n")
    f.write("\n")
    for iBin,val in enumerate(xSec_abs_val):
        # print "name_"+str(iBin)+" = cms_abs_xsec_rho_bin"+str(iBin)
        f.write("name_"+str(iBin)+" = cms_abs_xsec_rho_bin"+str(iBin)+"\n")
        # print "value_"+str(iBin)+" = "+str(val)
        f.write( "value_"+str(iBin)+" = "+str(val)+"\n")
        f.write("\n")

    # f.write("\n")
    f.write("[end estimates]")

    f.write("\n")
    f.write("\n")

    f.write("[not fitted]\n")
    f.write("\n")
    # print ""
    # print ""
    # print ""
    # now the funny part, the corr matrix...
    # let's start with the non-fitted extrapolation uncs.
    # print "                         extrapolation_unc"
    f.write("CMS_extrapolation_unc\n")
    for iBin_i,val_i in enumerate(xSec_abs_val):
        for iBin_j,val_j in enumerate(xSec_abs_val):
            if iBin_i == iBin_j:
                covEntry = covMatrix_extr[iBin_i][iBin_j]
                val = np.sqrt(covEntry)
                # print "cms_abs_xsec_rho_bin"+str(iBin_j)+"    "+str(val)
                f.write("cms_abs_xsec_rho_bin"+str(iBin_j)+"    "+str(val)+"\n")
    # print ""
    # print ""
    # print ""
    f.write("\n")
    f.write("[end not fitted]\n")


    f.write("\n")
    f.write("\n")

    f.write("[correlation matrix]\n")
    f.write("\n")
    # and now the real fun, all NP correlations
    for i, idxNuisance_i in enumerate(allNuisancesIdx):
        name_i = allNuisancesNames[i]
        if not "rate_ttj" in name_i:
            constraint_i = (abs(fit_obs_up[i])+abs(fit_obs_down[i]))/2.
            corrs_i = []
            for j, idxNuisance_j in enumerate(allNuisancesIdx):
                if j<=i:
                    corr = corrMatrix[i][j]
                    corrs_i.append(corr)
            strCorr = "   "
            for _ in corrs_i:
                strCorr += str(_)+"   "
            # print name_i+"   ("+str(constraint_i)+")"+strCorr
            f.write(name_i+" ("+str(constraint_i)+")"+strCorr+"\n")
        else:
            idxrate = int(float(name_i[-1]))
            constraint_i = (abs(xSec_abs_err_tot_up[idxrate])+abs(xSec_abs_err_tot_down[idxrate]))/2.
            corrs_i = []
            for j, idxNuisance_j in enumerate(allNuisancesIdx):
                if j<=i:
                    corr = corrMatrix[i][j]
                    corrs_i.append(corr)
            name_i = "cms_abs_xsec_rho_bin"+name_i[-1]
            strCorr = "   "
            for _ in corrs_i:
                strCorr += str(_)+"   "
            # print name_i,"   (",constraint_i,")",corrs_i
            # print name_i+"   ("+str(constraint_i)+")"+strCorr
            f.write(name_i+"   ("+str(constraint_i)+")"+strCorr+"\n")
        
    f.write("\n")
    f.write("[end correlation matrix]\n")

    f.write("\n")
    f.write("[systematics]\n")
    f.write("\n")
    f.write("[end systematics]\n")

    f.close()

# readCMSRootFile(cmsPathRoot)
readCMSCSVFiles("Inputs/")