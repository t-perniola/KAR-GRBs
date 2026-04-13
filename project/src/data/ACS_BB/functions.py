import subprocess, os, glob
import ROOT as M
import numpy as np
import pandas as pd
import sys, os
import matplotlib.pyplot as plt

from gdt.core.data_primitives import TimeBins
from gdt.core.plot import *
from bctools.analysis import BayesianBlocksLightcurve




#######################################################################################
'''
Functions:
- simulate_GRB: simulate a GRB event using Cosima
- extract: extract hits in ACS segments from the simulated file
- generate_lc_with_bkg: generate lightcurve by combining GRB data with background data
- analyze_lc: analyze the lightcurve data with bc-tools
'''
#######################################################################################




def simulate_GRB(source_file):
    """
    Simulate a Gamma-Ray Burst (GRB) event.
    """

    workdir = os.path.dirname(os.path.abspath(source_file))

    try:
        subprocess.run(['cosima', source_file], check=True)

        list_of_files = glob.glob(os.path.join(workdir, '*.sim.gz'))
        output_file = max(list_of_files, key=os.path.getmtime) # Get the most recent .sim.gz file
        
        return output_file
    
    except subprocess.CalledProcessError as e:
        print("Cosima failed!")
        print("Return code:", e.returncode)
        print("Command:", e.cmd)


def extract(GeometryName, CorrectionPath, SimFile, Shared, Mode, CorrThreshold = False, output_name = None, deactivate_threshold = False):
    """
    Extract hits in ACS segments from the file.
    Input:
    - GeometryName: path to geometry
    - CorrectionPath: path to correction file
    - SimFile: path to .sim file
    - Shared: if True, use ASIC-shared configuration
    - Mode: 'original', 'megalib', or 'correction'
    - CorrThreshold: if True, correct threshold based on crystal volume
    - output_name: name of the output file
    - deactivate_threshold: if True, set energy threshold to zero
    """

    # Load MEGAlib into ROOT
    M.gSystem.Load("$(MEGALIB)/lib/libMEGAlib.so")

    # Initialize MEGAlib
    G = M.MGlobal()
    G.Initialize()

    # We are good to go ...

    print("Ligth curve loaded !")

    # Load geometry:
    Geometry = M.MDGeometryQuest()
    if Geometry.ScanSetupFile(M.MString(GeometryName)) == True:
        print("Geometry " + GeometryName + " loaded!")
    else:
        print("Unable to load geometry " + GeometryName + " - Aborting!")
        quit()

    #by default Megalib is noising the energy. Uncomment below if you want to assume perfect energy resolution
    #Geometry.ActivateNoising(False);

    Reader = M.MFileEventsSim(Geometry)

    #Reader.ShowProgress()
    if Reader.Open(M.MString(SimFile)) == False:
        print("Unable to open file " + SimFile + ". Aborting!")
        quit()
    else:
        print("File " + SimFile + " opened successfully!")

    #per hit 
    BGO_Z0_0 = 0.
    BGO_Z0_1 = 0.
    BGO_Z0_2 = 0.
    BGO_Z0_3 = 0.
    BGO_Z0_4 = 0.
    BGO_Z1_4 = 0.
    BGO_Z1_3 = 0.
    BGO_Z1_2 = 0.
    BGO_Z1_1 = 0.
    BGO_Z1_0 = 0.

    BGO_X1_0 = 0.
    BGO_X1_1 = 0.
    BGO_X1_2 = 0.

    BGO_X0_0 = 0.
    BGO_X0_1 = 0.
    BGO_X0_2 = 0.

    BGO_Y1_0 = 0.
    BGO_Y1_1 = 0.
    BGO_Y1_2 = 0.

    BGO_Y0_0 = 0.
    BGO_Y0_1 = 0.
    BGO_Y0_2 = 0.

    time = 0.
    latX = 0.
    lonX = 0.
    latZ = 0.
    lonZ = 0.

    #per event
    Time = []
    Energy_Z1 = []
    Energy_Z0 = []
    Energy_X1 = []
    Energy_X0 = []
    Energy_Y1 = []
    Energy_Y0 = []
    LatX = []
    LonX = []
    LatZ = []
    LonZ = []

    # set the thresholds based on the crystal volumes
    params = np.genfromtxt(os.path.join(CorrectionPath, 'threshold_volume.dat'), unpack=True)
    #reference_threshold = params[3]
    reference_threshold = 80 # keV
    threshold = {}
    bgo_names = ["BGO_X0_0_Crystal", "BGO_X0_1_Crystal", "BGO_X0_2_Crystal", "BGO_X1_0_Crystal", "BGO_X1_1_Crystal", "BGO_X1_2_Crystal", "BGO_Y0_0_Crystal", "BGO_Y0_1_Crystal", "BGO_Y0_2_Crystal", "BGO_Y1_0_Crystal", "BGO_Y1_1_Crystal", "BGO_Y1_2_Crystal", "BGO_Z0_0_Crystal", "BGO_Z0_1_Crystal", "BGO_Z0_2_Crystal", "BGO_Z0_3_Crystal", "BGO_Z0_4_Crystal", "BGO_Z1_0_Crystal", "BGO_Z1_1_Crystal", "BGO_Z1_2_Crystal", "BGO_Z1_3_Crystal", "BGO_Z1_4_Crystal"]
    for bgo_name in bgo_names:
        if Mode == 'correction' and CorrThreshold is True:
            volume = Geometry.GetVolume(M.MString(bgo_name))
            size = volume.GetSize() # cm
            volume_det = size.X()*2 * size.Y()*2 * size.Z()*2 # cm^3
            corr_thre = params[0] * volume_det**2 + params[1] * volume_det + params[2]
            threshold[bgo_name] = reference_threshold * (1 + corr_thre)
        else:
            threshold[bgo_name] = 80 # keV
    # if deactivate_threshold is True, set all thresholds to zero
    if deactivate_threshold: 
        for bgo_name in bgo_names: 
            threshold[bgo_name] = 1e-6 # keV

    FirstHit = True

    while True:
        Event = Reader.GetNextEvent()

        if not Event:
            break
        M.SetOwnership(Event, True)
        
        #print(Event.GetID()) 
        #print(Event.GetTime().GetAsSeconds()) 
        #print(Event.ToString())
        #print(Event.GetNHTs())
        #print(Event.GetHTAt(0))
        #sys.exit() 	  
        if Event.GetNIAs() > 0:
        
            for i in range(Event.GetNHTs()):
                Hit = Event.GetHTAt(i)
                #print(Hit.ToString())
                if Hit.GetDetectorType() == 8:
                    
                    pos = Hit.GetPosition()
                    
                    x=pos.X()
                    y=pos.Y()
                    z=pos.Z()
                    original_energy = Hit.GetOriginalEnergy()
                    #megalib_energy = Hit.GetEnergy()
                    
                    

                    # getting the Voxel ID (x, y, z)
                    global_pos = Hit.GetOriginalPosition()
                    det = Geometry.GetDetector(global_pos)
                    volSeq = Geometry.GetVolumeSequence(global_pos)
                    pos_Voxel3D = volSeq.GetPositionInDetector()
                    Voxel3D = M.MDVoxel3D(det)
                    P = Voxel3D.GetGridPoint(pos_Voxel3D)
                    vx = P.GetXGrid()
                    vy = P.GetYGrid()
                    vz = P.GetZGrid()

                    detector = Geometry.GetDetector(pos).GetName()
                    
                    if volSeq.GetNVolumes() > 0:  # Ensure there is at least one volume
                        volume = volSeq.GetVolumeAt(volSeq.GetNVolumes() - 1)  # Get the last volume
                    size = volume.GetSize() # cm
                    volume_det = size.X()*2 * size.Y()*2 * size.Z()*2 # cm^3
                    #print(f"{detector.GetString()}, volume size: X = {size.X()*2} cm, Y = {size.Y()*2} cm, Z = {size.Z()*2} cm")

                    # applyng correction to energy
                    if Mode == "original":
                        energy = Hit.GetOriginalEnergy()
                    elif Mode == "megalib":
                        energy = Hit.GetEnergy()
                    elif Mode == 'correction':
                        correction_file = os.path.join(CorrectionPath, f'correction_file_{detector.GetString()}.dat')
                        if os.path.exists(correction_file):
                            cmx, cmy, cmz, m, q, a, b, c = np.genfromtxt(correction_file, usecols=(1,2,3,6,7,8,9,10), unpack=True)
                        else:
                            print("Correction file not found: ", correction_file)
                            sys.exit()
                        cond_voxel = ((cmx == vx) & (cmy == vy) & (cmz == vz))
                        m_xy, q_xy, a_xy, b_xy, c_xy = m[cond_voxel][0], q[cond_voxel][0], a[cond_voxel][0], b[cond_voxel][0], c[cond_voxel][0]
                        centroid = m_xy * original_energy + q_xy
                        fwhm = np.sqrt(a_xy**2 + b_xy**2 * original_energy + c_xy**2 * original_energy**2)
                        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
                        energy = np.random.normal(centroid, sigma)
                        #print("Correcting ", original_energy, " keV in ", corrected_energy, " keV")
                    else:
                        energy = Hit.GetEnergy()
                    

                    #bottom
                    if detector.GetString() == "BGO_Z0_0":
                        BGO_Z0_0+=(energy)
                    elif detector.GetString() == "BGO_Z0_1":
                        BGO_Z0_1+=(energy)
                    elif detector.GetString() == "BGO_Z0_2":
                        BGO_Z0_2+=(energy)
                    elif detector.GetString() == "BGO_Z0_3":
                        BGO_Z0_3+=(energy)
                    elif detector.GetString() == "BGO_Z0_4":
                        BGO_Z0_4+=(energy)
                    elif detector.GetString() == "BGO_Z1_4":
                        BGO_Z1_4+=(energy)
                    elif detector.GetString() == "BGO_Z1_3":
                        BGO_Z1_3+=(energy)
                    elif detector.GetString() == "BGO_Z1_2":
                        BGO_Z1_2+=(energy)
                    elif detector.GetString() == "BGO_Z1_1":
                        BGO_Z1_1+=(energy)
                    elif detector.GetString() == "BGO_Z1_0":
                        BGO_Z1_0+=(energy)
                        
                    
                    #Y pannel
                    elif detector.GetString() == "BGO_Y1_0":
                        BGO_Y1_0+=(energy)
                    elif detector.GetString() == "BGO_Y1_1":
                        BGO_Y1_1+=(energy)
                    elif detector.GetString() == "BGO_Y1_2":
                        BGO_Y1_2+=(energy)
                    
                    #Y neg pannel
                    elif detector.GetString() == "BGO_Y0_0":
                        BGO_Y0_0+=(energy)
                    elif detector.GetString() == "BGO_Y0_1":
                        BGO_Y0_1+=(energy)
                    elif detector.GetString() == "BGO_Y0_2":
                        BGO_Y0_2+=(energy)

                    #X pannel
                    elif detector.GetString() == "BGO_X1_0":
                        BGO_X1_0+=(energy)
                    elif detector.GetString() == "BGO_X1_1":
                        BGO_X1_1+=(energy)
                    elif detector.GetString() == "BGO_X1_2":
                        BGO_X1_2+=(energy)
                    
                    #X neg pannel
                    elif detector.GetString() == "BGO_X0_0":
                        BGO_X0_0+=(energy)
                    elif detector.GetString() == "BGO_X0_1":
                        BGO_X0_1+=(energy)
                    elif detector.GetString() == "BGO_X0_2":
                        BGO_X0_2+=(energy)
                    else:
                        print(detector)
                        print(str(x)+" "+str(y)+" "+str(z))
                        print("coordinate not found")
                        sys.exit()
            
                    if FirstHit:
                #time
                        time = ( Event.GetTime().GetAsSeconds())

                        # x axis of space craft pointing at GAL latitude
                        latX=(np.float32(Event.GetGalacticPointingXAxisLatitude()))
                    
                        # x axis of space craft pointing at GAL longitude
                        lonX=(np.float32(Event.GetGalacticPointingXAxisLongitude()))
                    
                        # z axis of space craft pointing at GAL latitude
                        latZ=(np.float32(Event.GetGalacticPointingZAxisLatitude()))
                    
                        # z axis of space craft pointing at GAL longitude
                        lonZ=(np.float32(Event.GetGalacticPointingZAxisLongitude()))

                        FirstHit = False
                
                    #print(time)
                    #print(energyBGO)
                    #print(Event.GetID()) 
                    #sys.exit()
        
                        
                        
        #end of Event	    
        #empty the memory took by the event
        #Check if sum of hits in one crystal >= 80keV
        
        #bot neg
        if BGO_Z0_0 >= threshold["BGO_Z0_0_Crystal"]:
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(BGO_Z0_0))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)        
        if BGO_Z0_1 >= threshold["BGO_Z0_1_Crystal"]:
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(BGO_Z0_1))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        if BGO_Z0_2 >= threshold["BGO_Z0_2_Crystal"]:
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(BGO_Z0_2))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        if BGO_Z0_3 >= threshold["BGO_Z0_3_Crystal"]:
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(BGO_Z0_3))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
            
            
        if BGO_Z0_4 >= threshold["BGO_Z0_4_Crystal"]:
            Energy_Z1.append(np.float32(0))
            
            if Shared == True:
                Energy_Y0.append(np.float32(BGO_Z0_4))  ## ASIC shared
                Energy_Z0.append(np.float32(0))
            else:
                Energy_Z0.append(np.float32(BGO_Z0_4))
                Energy_Y0.append(np.float32(0))
                
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0)) 
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        #bot
        if BGO_Z1_4 >= threshold["BGO_Z1_4_Crystal"]:
            
            if Shared == True:
                Energy_Y1.append(np.float32(BGO_Z1_4)) ## ASIC shared
                Energy_Z1.append(np.float32(0))
            else:
                Energy_Z1.append(np.float32(BGO_Z1_4))
                Energy_Y1.append(np.float32(0))
                
            Energy_Z0.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y0.append(np.float32(0))  
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
            
        if BGO_Z1_3 >= threshold["BGO_Z1_3_Crystal"]:
            Energy_Z0.append(np.float32(0))
            Energy_Z1.append(np.float32(BGO_Z1_3))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        if BGO_Z1_2 >= threshold["BGO_Z1_2_Crystal"]:
            Energy_Z0.append(np.float32(0))
            Energy_Z1.append(np.float32(BGO_Z1_2))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        if BGO_Z1_1 >= threshold["BGO_Z1_1_Crystal"]:
            Energy_Z0.append(np.float32(0))
            Energy_Z1.append(np.float32(BGO_Z1_1))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        if BGO_Z1_0 >= threshold["BGO_Z1_0_Crystal"]:
            Energy_Z0.append(np.float32(0))
            Energy_Z1.append(np.float32(BGO_Z1_0))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_X1_0 >= threshold["BGO_X1_0_Crystal"]:
            Energy_X1.append(np.float32(BGO_X1_0))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_X1_1 >= threshold["BGO_X1_1_Crystal"]:
            Energy_X1.append(np.float32(BGO_X1_1))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_X1_2 >= threshold["BGO_X1_2_Crystal"]:
            Energy_X1.append(np.float32(BGO_X1_2))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_Y1_0 >= threshold["BGO_Y1_0_Crystal"]:
            Energy_Y1.append(np.float32(BGO_Y1_0))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        
        if BGO_Y1_1 >= threshold["BGO_Y1_1_Crystal"]:
            Energy_Y1.append(np.float32(BGO_Y1_1))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        
        if BGO_Y1_2 >= threshold["BGO_Y1_2_Crystal"]:
            Energy_Y1.append(np.float32(BGO_Y1_2))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_X0_0 >= threshold["BGO_X0_0_Crystal"]:
            Energy_X0.append(np.float32(BGO_X0_0))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_X0_1 >= threshold["BGO_X0_1_Crystal"]:
            Energy_X0.append(np.float32(BGO_X0_1))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_X0_2 >= threshold["BGO_X0_2_Crystal"]:
            Energy_X0.append(np.float32(BGO_X0_2))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_Y0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)

        if BGO_Y0_0 >= threshold["BGO_Y0_0_Crystal"]:
            Energy_Y0.append(np.float32(BGO_Y0_0))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        
        if BGO_Y0_1 >= threshold["BGO_Y0_1_Crystal"]:
            Energy_Y0.append(np.float32(BGO_Y0_1))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)
        
        if BGO_Y0_2 >= threshold["BGO_Y0_2_Crystal"]:
            Energy_Y0.append(np.float32(BGO_Y0_2))
            Energy_Z1.append(np.float32(0))
            Energy_Z0.append(np.float32(0))
            Energy_Y1.append(np.float32(0))
            Energy_X1.append(np.float32(0))
            Energy_X0.append(np.float32(0))
            Time.append(time)
            LatX.append(latX)
            LonX.append(lonX)
            LatZ.append(latZ)
            LonZ.append(lonZ)


        BGO_Z0_0 = 0.
        BGO_Z0_1 = 0.
        BGO_Z0_2 = 0.
        BGO_Z0_3 = 0.
        BGO_Z0_4 = 0.
        BGO_Z1_4 = 0.
        BGO_Z1_3 = 0.
        BGO_Z1_2 = 0.
        BGO_Z1_1 = 0.
        BGO_Z1_0 = 0.
        BGO_X1_0 = 0.
        BGO_X1_1 = 0.
        BGO_X1_2 = 0.
        BGO_X0_0 = 0.
        BGO_X0_1 = 0.
        BGO_X0_2 = 0.
        BGO_Y1_0 = 0.
        BGO_Y1_1 = 0.
        BGO_Y1_2 = 0.
        BGO_Y0_0 = 0.
        BGO_Y0_1 = 0.
        BGO_Y0_2 = 0.

        Event = 0
        FirstHit = True

        #end of Event loop

    output_file = save_to_csv(SimFile, Shared, Time, Energy_Z1, Energy_Z0, Energy_X1, Energy_X0, Energy_Y1, Energy_Y0, LatX, LonX, LatZ, LonZ, Mode, output_name)

    return output_file
    

def save_to_csv(SimFile, Shared, Time, Energy_Z1, Energy_Z0, Energy_X1, Energy_X0, Energy_Y1, Energy_Y0, LatX, LonX, LatZ, LonZ, Mode, output_name = None):
    """
    Save the extracted data to a CSV file.
    """

    df = pd.DataFrame({"timestamp[s]": Time,"bgo_z1[keV]":Energy_Z1,"bgo_z0[keV]":Energy_Z0
                    ,"bgo_x1[keV]":Energy_X1,"bgo_x0[keV]":Energy_X0,
                    "bgo_y1[keV]":Energy_Y1,"bgo_y0[keV]":Energy_Y0,
                    "latX":LatX,"lonX":LonX,"latZ":LatZ,"lonZ":LonZ })

    df['datetime']= pd.to_datetime(df["timestamp[s]"],unit="s")
    df.set_index("datetime",inplace=True)

    base_name = os.path.splitext(os.path.basename(SimFile))[0]
    base_name = base_name.split(".")[0]
    suffix = ""

    if Shared:
        suffix += "_shared"

    if Mode == "original":
        suffix += "_original"
    elif Mode == "megalib":
        suffix += "_megalib"
    elif Mode == "correction":
        suffix += "_correction"
    else:
        suffix += "_megalib"

    if output_name == None: output_file = os.path.join(
        os.path.dirname(os.path.abspath(SimFile)),
        base_name + suffix + ".csv"
    )
    else: output_file = os.path.join(
        os.path.dirname(os.path.abspath(SimFile)),
        output_name + ".csv"
    )
    df.to_csv(output_file, index=True) 

    bgo_data_seconds = {}
    for col in df.columns:
        if col.startswith("bgo_"):
            filtered_data = df[df[col] > 0][[col]]
            # Convert datetime index to seconds since epoch
            times_in_seconds = filtered_data.index.map(lambda x: x.timestamp())
            bgo_data_seconds[col] = np.array(times_in_seconds)

    output_file_counts = open(output_file.replace(".csv","_counts.dat"),"w")
    for col, times_in_seconds in bgo_data_seconds.items():
        output_file_counts.write(col+" "+str(len(times_in_seconds))+"\n")
    output_file_counts.close()

    return output_file

def generate_lc_with_bkg(file_input, path_bkg, bin_width, lc_length, pre_trigger_time, t_start=None, t_stop=None):
    """
    Combine GRB data to background data to generate a ligthcurve
    Input:
    - file_input: input file with hits in each ACS segment
    - path_bkg: path to the background files
    - bin_width: bin width for the lightcurve (in seconds)
    - lc_length: length of the lightcurve (in seconds)
    - pre_trigger_time: pre-trigger time to be added to the GRB times (in seconds)
    - t_start: start time for the background selection (in seconds since epoch). If None, a random time is selected.
    - t_stop: stop time for the background selection (in seconds since epoch). If None, a random time is selected.
    """

    # Randomly select background portion
    time_SAA, SAA = [], []
    with open(os.path.join(path_bkg, 'SAAproton_4MeVto2000.dat'), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('DP'):
                columns = line.split()
                time_SAA.append(float(columns[1]))
                SAA.append(float(columns[2]))
    time_SAA, SAA = np.array(time_SAA), np.array(SAA)
    transition_SAA_in = np.where((SAA[:-1] == 0) & (SAA[1:] > 0))[0]
    transition_SAA_out = np.where((SAA[:-1] > 0) & (SAA[1:] == 0))[0]
    transition_times_SAA_in = time_SAA[transition_SAA_in]
    transition_times_SAA_out = time_SAA[transition_SAA_out]
    transition_times_SAA_out = np.append(transition_times_SAA_out, time_SAA[-1])

    t_min = 1835517300.0 # s
    t_max = 1835573699.0 # s
    len_lc = lc_length
    if t_start is None or t_stop is None:
        while True:
            random_time = np.random.uniform(t_min, t_max)
            t_start = random_time
            t_stop = random_time + len_lc
            is_tstart_in_SAA = (transition_times_SAA_in <= t_start) & (t_start <= transition_times_SAA_out)
            is_tstop_in_SAA = (transition_times_SAA_in <= t_stop) & (t_stop <= transition_times_SAA_out)
            if not (np.any(is_tstart_in_SAA) or np.any(is_tstop_in_SAA)):
                break
    else:
        if t_start < t_min or t_stop > t_max:
            raise ValueError("t_start and t_stop must be within the range of SAA data.")
        if t_start >= t_stop:
            raise ValueError("t_start must be less than t_stop.")
        is_tstart_in_SAA = (transition_times_SAA_in <= t_start) & (t_start <= transition_times_SAA_out)
        is_tstop_in_SAA = (transition_times_SAA_in <= t_stop) & (t_stop <= transition_times_SAA_out)
        if np.any(is_tstart_in_SAA) or np.any(is_tstop_in_SAA):
            raise ValueError("t_start and/or t_stop are within SAA periods.")

    components = ['AlbedoNeutrons', 'AlbedoPhotons', 'CosmicPhotons', 'PrimaryAlphas', 'PrimaryElectrons', 'PrimaryPositrons', 'PrimaryProtons', 'SecondaryElectrons', 'SecondaryPositrons', 'SecondaryProtons']
    files_0 = [component + '_BGOhit_Total_0.npz' for component in components]
    files_0.append('SAAprotons_BGOhit_Total_530km_0.npz')
    files_0 = [os.path.join(path_bkg, file) for file in files_0]
    files_1 = [component + '_BGOhit_Total_1.npz' for component in components]
    files_1.append('SAAprotons_BGOhit_Total_530km_1.npz')
    files_1 = [os.path.join(path_bkg, file) for file in files_1]
    z0_times_cut0, z1_times_cut0, x0_times_cut0, x1_times_cut0, y0_times_cut0, y1_times_cut0 = [], [], [], [], [], []
    z0_times_cut1, z1_times_cut1, x0_times_cut1, x1_times_cut1, y0_times_cut1, y1_times_cut1 = [], [], [], [], [], []
    for file0, file1 in zip(files_0, files_1):
        # -- Channel 0 -- #
        # Load data
        data = np.load(file0)
        z0_events = data['z0_events']
        z1_events = data['z1_events']
        x0_events = data['x0_events']
        x1_events = data['x1_events']
        y0_events = data['y0_events']
        y1_events = data['y1_events']

        # Take times from events
        z0_times = z0_events[:, 0]
        z1_times = z1_events[:, 0]
        x0_times = x0_events[:, 0]
        x1_times = x1_events[:, 0]
        y0_times = y0_events[:, 0]
        y1_times = y1_events[:, 0]

        # Cut events
        z0_events_cut0 = z0_events[(t_start <= z0_times) & (z0_times < t_stop)]
        z1_events_cut0 = z1_events[(t_start <= z1_times) & (z1_times < t_stop)]
        x0_events_cut0 = x0_events[(t_start <= x0_times) & (x0_times < t_stop)]
        x1_events_cut0 = x1_events[(t_start <= x1_times) & (x1_times < t_stop)]
        y0_events_cut0 = y0_events[(t_start <= y0_times) & (y0_times < t_stop)]
        y1_events_cut0 = y1_events[(t_start <= y1_times) & (y1_times < t_stop)]

        # Cut times
        z0_times_cut0.extend(z0_times[(t_start <= z0_times) & (z0_times < t_stop)])
        z1_times_cut0.extend(z1_times[(t_start <= z1_times) & (z1_times < t_stop)])
        x0_times_cut0.extend(x0_times[(t_start <= x0_times) & (x0_times < t_stop)])
        x1_times_cut0.extend(x1_times[(t_start <= x1_times) & (x1_times < t_stop)])
        y0_times_cut0.extend(y0_times[(t_start <= y0_times) & (y0_times < t_stop)])
        y1_times_cut0.extend(y1_times[(t_start <= y1_times) & (y1_times < t_stop)])

        # -- Channel 1 -- #
        # Load data
        data = np.load(file1)
        z0_events = data['z0_events']
        z1_events = data['z1_events']
        x0_events = data['x0_events']
        x1_events = data['x1_events']
        y0_events = data['y0_events']
        y1_events = data['y1_events']

        # Take times from events
        z0_times = z0_events[:, 0]
        z1_times = z1_events[:, 0]
        x0_times = x0_events[:, 0]
        x1_times = x1_events[:, 0]
        y0_times = y0_events[:, 0]
        y1_times = y1_events[:, 0]

        # Cut events
        z0_events_cut1 = z0_events[(t_start <= z0_times) & (z0_times < t_stop)]
        z1_events_cut1 = z1_events[(t_start <= z1_times) & (z1_times < t_stop)]
        x0_events_cut1 = x0_events[(t_start <= x0_times) & (x0_times < t_stop)]
        x1_events_cut1 = x1_events[(t_start <= x1_times) & (x1_times < t_stop)]
        y0_events_cut1 = y0_events[(t_start <= y0_times) & (y0_times < t_stop)]
        y1_events_cut1 = y1_events[(t_start <= y1_times) & (y1_times < t_stop)]

        # Cut times
        z0_times_cut1.extend(z0_times[(t_start <= z0_times) & (z0_times < t_stop)])
        z1_times_cut1.extend(z1_times[(t_start <= z1_times) & (z1_times < t_stop)])
        x0_times_cut1.extend(x0_times[(t_start <= x0_times) & (x0_times < t_stop)])
        x1_times_cut1.extend(x1_times[(t_start <= x1_times) & (x1_times < t_stop)])
        y0_times_cut1.extend(y0_times[(t_start <= y0_times) & (y0_times < t_stop)])
        y1_times_cut1.extend(y1_times[(t_start <= y1_times) & (y1_times < t_stop)])

    # Shift GRB times
    df = pd.read_csv(file_input)
    time, z0, z1, x0, x1, y0, y1 = df.iloc[:, 1:8].values.T
    time = time + t_start + pre_trigger_time

    # Remove GRB events if t > tstop
    time = time[time < t_stop]
    z0 = z0[time < t_stop]
    z1 = z1[time < t_stop]
    x0 = x0[time < t_stop]
    x1 = x1[time < t_stop]
    y0 = y0[time < t_stop]
    y1 = y1[time < t_stop]

    # Filter times for each segment
    times_z0_0 = time[(z0 > 80) & (z0 < 2000)]
    times_z1_0 = time[(z1 > 80) & (z1 < 2000)]
    times_x0_0 = time[(x0 > 80) & (x0 < 2000)]
    times_x1_0 = time[(x1 > 80) & (x1 < 2000)]
    times_y0_0 = time[(y0 > 80) & (y0 < 2000)]
    times_y1_0 = time[(y1 > 80) & (y1 < 2000)]
    times_z0_1 = time[z0 > 2000]
    times_z1_1 = time[z1 > 2000]
    times_x0_1 = time[x0 > 2000]
    times_x1_1 = time[x1 > 2000]
    times_y0_1 = time[y0 > 2000]
    times_y1_1 = time[y1 > 2000]

    # Filter energies for each segment
    energies_z0_0 = z0[(z0 > 80) & (z0 < 2000)]
    energies_z1_0 = z1[(z1 > 80) & (z1 < 2000)]
    energies_x0_0 = x0[(x0 > 80) & (x0 < 2000)]
    energies_x1_0 = x1[(x1 > 80) & (x1 < 2000)]
    energies_y0_0 = y0[(y0 > 80) & (y0 < 2000)]
    energies_y1_0 = y1[(y1 > 80) & (y1 < 2000)]
    energies_z0_1 = z0[z0 > 2000]
    energies_z1_1 = z1[z1 > 2000]
    energies_x0_1 = x0[x0 > 2000]
    energies_x1_1 = x1[x1 > 2000]
    energies_y0_1 = y0[y0 > 2000]
    energies_y1_1 = y1[y1 > 2000]

    # Generate event lists for GRB
    z0_events_signal0 = np.array([(t, e) for t, e in zip(times_z0_0, energies_z0_0)])
    z1_events_signal0 = np.array([(t, e) for t, e in zip(times_z1_0, energies_z1_0)])
    x0_events_signal0 = np.array([(t, e) for t, e in zip(times_x0_0, energies_x0_0)])
    x1_events_signal0 = np.array([(t, e) for t, e in zip(times_x1_0, energies_x1_0)])
    y0_events_signal0 = np.array([(t, e) for t, e in zip(times_y0_0, energies_y0_0)])
    y1_events_signal0 = np.array([(t, e) for t, e in zip(times_y1_0, energies_y1_0)])
    z0_events_signal1 = np.array([(t, e) for t, e in zip(times_z0_1, energies_z0_1)])
    z1_events_signal1 = np.array([(t, e) for t, e in zip(times_z1_1, energies_z1_1)])
    x0_events_signal1 = np.array([(t, e) for t, e in zip(times_x0_1, energies_x0_1)])
    x1_events_signal1 = np.array([(t, e) for t, e in zip(times_x1_1, energies_x1_1)])
    y0_events_signal1 = np.array([(t, e) for t, e in zip(times_y0_1, energies_y0_1)])
    y1_events_signal1 = np.array([(t, e) for t, e in zip(times_y1_1, energies_y1_1)])

    # Generate event lists for signal + background
    def safe_concat(a, b):
        return a if b.size == 0 else np.concatenate((a, b), axis=0)
    z0_events_0 = safe_concat(z0_events_cut0, z0_events_signal0)
    z1_events_0 = safe_concat(z1_events_cut0, z1_events_signal0)
    x0_events_0 = safe_concat(x0_events_cut0, x0_events_signal0)
    x1_events_0 = safe_concat(x1_events_cut0, x1_events_signal0)
    y0_events_0 = safe_concat(y0_events_cut0, y0_events_signal0)
    y1_events_0 = safe_concat(y1_events_cut0, y1_events_signal0)
    z0_events_1 = safe_concat(z0_events_cut1, z0_events_signal1)
    z1_events_1 = safe_concat(z1_events_cut1, z1_events_signal1)
    x0_events_1 = safe_concat(x0_events_cut1, x0_events_signal1)
    x1_events_1 = safe_concat(x1_events_cut1, x1_events_signal1)
    y0_events_1 = safe_concat(y0_events_cut1, y0_events_signal1)
    y1_events_1 = safe_concat(y1_events_cut1, y1_events_signal1)

    times_evt0, z0_evt0, z1_evt0, x0_evt0, x1_evt0, y0_evt0, y1_evt0 = [], [], [], [], [], [], []
    for z00 in z0_events_0:
        times_evt0.append(z00[0])
        z0_evt0.append(z00[1])
        z1_evt0.append(0)
        x0_evt0.append(0)
        x1_evt0.append(0)
        y0_evt0.append(0)
        y1_evt0.append(0)
    for z10 in z1_events_0:
        times_evt0.append(z10[0])
        z0_evt0.append(0)
        z1_evt0.append(z10[1])
        x0_evt0.append(0)
        x1_evt0.append(0)
        y0_evt0.append(0)
        y1_evt0.append(0)
    for x00 in x0_events_0:
        times_evt0.append(x00[0])
        z0_evt0.append(0)
        z1_evt0.append(0)
        x0_evt0.append(x00[1])
        x1_evt0.append(0)
        y0_evt0.append(0)
        y1_evt0.append(0)
    for x10 in x1_events_0:
        times_evt0.append(x10[0])
        z0_evt0.append(0)
        z1_evt0.append(0)
        x0_evt0.append(0)
        x1_evt0.append(x10[1])
        y0_evt0.append(0)
        y1_evt0.append(0)
    for y00 in y0_events_0:
        times_evt0.append(y00[0])
        z0_evt0.append(0)
        z1_evt0.append(0)
        x0_evt0.append(0)
        x1_evt0.append(0)
        y0_evt0.append(y00[1])
        y1_evt0.append(0)
    for y10 in y1_events_0:
        times_evt0.append(y10[0])
        z0_evt0.append(0)
        z1_evt0.append(0)
        x0_evt0.append(0)
        x1_evt0.append(0)
        y0_evt0.append(0)
        y1_evt0.append(y10[1])
    
    times_evt1, z0_evt1, z1_evt1, x0_evt1, x1_evt1, y0_evt1, y1_evt1 = [], [], [], [], [], [], []
    for z01 in z0_events_1:
        times_evt1.append(z01[0])
        z0_evt1.append(z01[1])
        z1_evt1.append(0)
        x0_evt1.append(0)
        x1_evt1.append(0)
        y0_evt1.append(0)
        y1_evt1.append(0)
    for z11 in z1_events_1:
        times_evt1.append(z11[0])
        z0_evt1.append(0)
        z1_evt1.append(z11[1])
        x0_evt1.append(0)
        x1_evt1.append(0)
        y0_evt1.append(0)
        y1_evt1.append(0)
    for x01 in x0_events_1:
        times_evt1.append(x01[0])
        z0_evt1.append(0)
        z1_evt1.append(0)
        x0_evt1.append(x01[1])
        x1_evt1.append(0)
        y0_evt1.append(0)
        y1_evt1.append(0)
    for x11 in x1_events_1:
        times_evt1.append(x11[0])
        z0_evt1.append(0)
        z1_evt1.append(0)
        x0_evt1.append(0)
        x1_evt1.append(x11[1])
        y0_evt1.append(0)
        y1_evt1.append(0)
    for y01 in y0_events_1:
        times_evt1.append(y01[0])
        z0_evt1.append(0)
        z1_evt1.append(0)
        x0_evt1.append(0)
        x1_evt1.append(0)
        y0_evt1.append(y01[1])
        y1_evt1.append(0)
    for y11 in y1_events_1:
        times_evt1.append(y11[0])
        z0_evt1.append(0)
        z1_evt1.append(0)
        x0_evt1.append(0)
        x1_evt1.append(0)
        y0_evt1.append(0)
        y1_evt1.append(y11[1])
    
    # Convert lists to numpy arrays
    times_evt0, z0_evt0, z1_evt0, x0_evt0, x1_evt0, y0_evt0, y1_evt0 = np.array(times_evt0), np.array(z0_evt0), np.array(z1_evt0), np.array(x0_evt0), np.array(x1_evt0), np.array(y0_evt0), np.array(y1_evt0)
    times_evt1, z0_evt1, z1_evt1, x0_evt1, x1_evt1, y0_evt1, y1_evt1 = np.array(times_evt1), np.array(z0_evt1), np.array(z1_evt1), np.array(x0_evt1), np.array(x1_evt1), np.array(y0_evt1), np.array(y1_evt1)

    # Sort all events by time
    sorted_indices = np.argsort(times_evt0)
    times_evt0 = times_evt0[sorted_indices]
    z0_evt0    = z0_evt0[sorted_indices]
    z1_evt0    = z1_evt0[sorted_indices]
    x0_evt0    = x0_evt0[sorted_indices]
    x1_evt0    = x1_evt0[sorted_indices]
    y0_evt0    = y0_evt0[sorted_indices]
    y1_evt0    = y1_evt0[sorted_indices]

    sorted_indices = np.argsort(times_evt1)
    times_evt1 = times_evt1[sorted_indices]
    z0_evt1    = z0_evt1[sorted_indices]
    z1_evt1    = z1_evt1[sorted_indices]
    x0_evt1    = x0_evt1[sorted_indices]
    x1_evt1    = x1_evt1[sorted_indices]
    y0_evt1    = y0_evt1[sorted_indices]
    y1_evt1    = y1_evt1[sorted_indices]

    # Bin GRB to lightcurve
    bin_edges = np.arange(t_start, t_stop + bin_width, bin_width)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

    counts_GRB_z0_0, _ = np.histogram(times_z0_0, bins=bin_edges) # signal
    counts_GRB_z1_0, _ = np.histogram(times_z1_0, bins=bin_edges) # signal
    counts_GRB_x0_0, _ = np.histogram(times_x0_0, bins=bin_edges) # signal
    counts_GRB_x1_0, _ = np.histogram(times_x1_0, bins=bin_edges) # signal
    counts_GRB_y0_0, _ = np.histogram(times_y0_0, bins=bin_edges) # signal
    counts_GRB_y1_0, _ = np.histogram(times_y1_0, bins=bin_edges) # signal
    counts_bkg_z0_0, _ = np.histogram(z0_times_cut0, bins=bin_edges) # background
    counts_bkg_z1_0, _ = np.histogram(z1_times_cut0, bins=bin_edges) # background
    counts_bkg_x0_0, _ = np.histogram(x0_times_cut0, bins=bin_edges) # background
    counts_bkg_x1_0, _ = np.histogram(x1_times_cut0, bins=bin_edges) # background
    counts_bkg_y0_0, _ = np.histogram(y0_times_cut0, bins=bin_edges) # background
    counts_bkg_y1_0, _ = np.histogram(y1_times_cut0, bins=bin_edges) # background

    counts_GRB_z0_1, _ = np.histogram(times_z0_1, bins=bin_edges) # signal
    counts_GRB_z1_1, _ = np.histogram(times_z1_1, bins=bin_edges) # signal
    counts_GRB_x0_1, _ = np.histogram(times_x0_1, bins=bin_edges) # signal
    counts_GRB_x1_1, _ = np.histogram(times_x1_1, bins=bin_edges) # signal
    counts_GRB_y0_1, _ = np.histogram(times_y0_1, bins=bin_edges) # signal
    counts_GRB_y1_1, _ = np.histogram(times_y1_1, bins=bin_edges) # signal
    counts_bkg_z0_1, _ = np.histogram(z0_times_cut1, bins=bin_edges) # background
    counts_bkg_z1_1, _ = np.histogram(z1_times_cut1, bins=bin_edges) # background
    counts_bkg_x0_1, _ = np.histogram(x0_times_cut1, bins=bin_edges) # background
    counts_bkg_x1_1, _ = np.histogram(x1_times_cut1, bins=bin_edges) # background
    counts_bkg_y0_1, _ = np.histogram(y0_times_cut1, bins=bin_edges) # background
    counts_bkg_y1_1, _ = np.histogram(y1_times_cut1, bins=bin_edges) # background


    # Sum signal and bkg
    signal_z0_0 = counts_GRB_z0_0 + counts_bkg_z0_0
    signal_z1_0 = counts_GRB_z1_0 + counts_bkg_z1_0
    signal_x0_0 = counts_GRB_x0_0 + counts_bkg_x0_0
    signal_x1_0 = counts_GRB_x1_0 + counts_bkg_x1_0
    signal_y0_0 = counts_GRB_y0_0 + counts_bkg_y0_0
    signal_y1_0 = counts_GRB_y1_0 + counts_bkg_y1_0

    signal_z0_1 = counts_GRB_z0_1 + counts_bkg_z0_1
    signal_z1_1 = counts_GRB_z1_1 + counts_bkg_z1_1
    signal_x0_1 = counts_GRB_x0_1 + counts_bkg_x0_1
    signal_x1_1 = counts_GRB_x1_1 + counts_bkg_x1_1
    signal_y0_1 = counts_GRB_y0_1 + counts_bkg_y0_1
    signal_y1_1 = counts_GRB_y1_1 + counts_bkg_y1_1


    # -- Output -- #

    base_name = os.path.basename(file_input).replace('.csv', '')
    output_file_events_c0 = os.path.join(os.path.dirname(file_input), base_name + '_c0.evt')
    output_file_events_c1 = os.path.join(os.path.dirname(file_input), base_name + '_c1.evt')
    output_file_lc_c0 = os.path.join(os.path.dirname(file_input), base_name + '_c0.lc')
    output_file_lc_c1 = os.path.join(os.path.dirname(file_input), base_name + '_c1.lc')

    # Event list
    with open(output_file_events_c0, 'w') as f:
        f.write('# Time[s]    z0[keV]    z1[keV]    x0[keV]    x1[keV]    y0[keV]    y1[keV]\n')
        for i in range(len(times_evt0)):
            f.write(f"{times_evt0[i]:<20}  "
                    f"{z0_evt0[i]:<20}  "
                    f"{z1_evt0[i]:<20}  "
                    f"{x0_evt0[i]:<20}  "
                    f"{x1_evt0[i]:<20}  "
                    f"{y0_evt0[i]:<20}  "
                    f"{y1_evt0[i]:<20}\n")
    with open(output_file_events_c1, 'w') as f:
        f.write('# Time[s]    z0[keV]    z1[keV]    x0[keV]    x1[keV]    y0[keV]    y1[keV]\n')
        for i in range(len(times_evt1)):
            f.write(f"{times_evt1[i]:<20}  "
                    f"{z0_evt1[i]:<20}  "
                    f"{z1_evt1[i]:<20}  "
                    f"{x0_evt1[i]:<20}  "
                    f"{x1_evt1[i]:<20}  "
                    f"{y0_evt1[i]:<20}  "
                    f"{y1_evt1[i]:<20}\n")
    
    # Lightcurve
    with open(output_file_lc_c0, 'w') as f:
        f.write('# Time[s]    z0    z1    x0    x1    y0    y1\n')
        for i in range(len(bin_centers)):
            f.write(f"{bin_centers[i]:<20}  "
                    f"{signal_z0_0[i]:<20}  "
                    f"{signal_z1_0[i]:<20}  "
                    f"{signal_x0_0[i]:<20}  "
                    f"{signal_x1_0[i]:<20}  "
                    f"{signal_y0_0[i]:<20}  "
                    f"{signal_y1_0[i]:<20}\n")
    with open(output_file_lc_c1, 'w') as f:
        f.write('# Time[s]    z0    z1    x0    x1    y0    y1\n')
        for i in range(len(bin_centers)):
            f.write(f"{bin_centers[i]:<20}  "
                    f"{signal_z0_1[i]:<20}  "
                    f"{signal_z1_1[i]:<20}  "
                    f"{signal_x0_1[i]:<20}  "
                    f"{signal_x1_1[i]:<20}  "
                    f"{signal_y0_1[i]:<20}  "
                    f"{signal_y1_1[i]:<20}\n")


    # Create a 3x1 grid of subplots
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15, 8))

    axs = axs.flatten()

    # Plot each signal in a different subplot
    axs[0].errorbar(bin_centers[signal_z0_0 > 0], signal_z0_0[signal_z0_0 > 0], 
                    yerr = np.sqrt(signal_z0_0[signal_z0_0 > 0]), fmt='.', label='bottom1')
    axs[0].set_title('bottom1')
    axs[0].set_xlabel('Time [s]')
    axs[0].set_ylabel('Counts')

    axs[1].errorbar(bin_centers[signal_z1_0 > 0], signal_z1_0[signal_z1_0 > 0], 
                    yerr = np.sqrt(signal_z1_0[signal_z1_0 > 0]), fmt='.', label='bottom2')
    axs[1].set_title('bottom2')
    axs[1].set_xlabel('Time [s]')
    axs[1].set_ylabel('Counts')

    axs[2].errorbar(bin_centers[signal_x0_0 > 0], signal_x0_0[signal_x0_0 > 0], 
                    yerr = np.sqrt(signal_x0_0[signal_x0_0 > 0]), fmt='.', label='x0')
    axs[2].set_title('x0')
    axs[2].set_xlabel('Time [s]')
    axs[2].set_ylabel('Counts')

    axs[3].errorbar(bin_centers[signal_x1_0 > 0], signal_x1_0[signal_x1_0 > 0], 
                    yerr = np.sqrt(signal_x1_0[signal_x1_0 > 0]), fmt='.', label='x1')
    axs[3].set_title('x1')
    axs[3].set_xlabel('Time [s]')
    axs[3].set_ylabel('Counts')

    axs[4].errorbar(bin_centers[signal_y0_0 > 0], signal_y0_0[signal_y0_0 > 0], 
                    yerr = np.sqrt(signal_y0_0[signal_y0_0 > 0]), fmt='.', label='y0')
    axs[4].set_title('y0')
    axs[4].set_xlabel('Time [s]')
    axs[4].set_ylabel('Counts')

    axs[5].errorbar(bin_centers[signal_y1_0 > 0], signal_y1_0[signal_y1_0 > 0], 
                    yerr = np.sqrt(signal_y1_0[signal_y1_0 > 0]), fmt='.', label='y1')
    axs[5].set_title('y1')
    axs[5].set_xlabel('Time [s]')
    axs[5].set_ylabel('Counts')

    fig.tight_layout()

    fig.savefig(os.path.join(os.path.dirname(file_input), base_name+'_c0.pdf'))
    #plt.show()

    # Create a 3x1 grid of subplots
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15, 8))

    axs = axs.flatten()

    # Plot each signal in a different subplot
    axs[0].errorbar(bin_centers[signal_z0_1 > 0], signal_z0_1[signal_z0_1 > 0], 
                    yerr = np.sqrt(signal_z0_1[signal_z0_1 > 0]), fmt='.', label='bottom1')
    axs[0].set_title('bottom1')
    axs[0].set_xlabel('Time [s]')
    axs[0].set_ylabel('Counts')

    axs[1].errorbar(bin_centers[signal_z1_1 > 0], signal_z1_1[signal_z1_1 > 0], 
                    yerr = np.sqrt(signal_z1_1[signal_z1_1 > 0]), fmt='.', label='bottom2')
    axs[1].set_title('bottom2')
    axs[1].set_xlabel('Time [s]')
    axs[1].set_ylabel('Counts')

    axs[2].errorbar(bin_centers[signal_x0_1 > 0], signal_x0_1[signal_x0_1 > 0], 
                    yerr = np.sqrt(signal_x0_1[signal_x0_1 > 0]), fmt='.', label='x0')
    axs[2].set_title('x0')
    axs[2].set_xlabel('Time [s]')
    axs[2].set_ylabel('Counts')

    axs[3].errorbar(bin_centers[signal_x1_1 > 0], signal_x1_1[signal_x1_1 > 0], 
                    yerr = np.sqrt(signal_x1_1[signal_x1_1 > 0]), fmt='.', label='x1')
    axs[3].set_title('x1')
    axs[3].set_xlabel('Time [s]')
    axs[3].set_ylabel('Counts')

    axs[4].errorbar(bin_centers[signal_y0_1 > 0], signal_y0_1[signal_y0_1 > 0], 
                    yerr = np.sqrt(signal_y0_1[signal_y0_1 > 0]), fmt='.', label='y0')
    axs[4].set_title('y0')
    axs[4].set_xlabel('Time [s]')
    axs[4].set_ylabel('Counts')

    axs[5].errorbar(bin_centers[signal_y1_1 > 0], signal_y1_1[signal_y1_1 > 0], 
                    yerr = np.sqrt(signal_y1_1[signal_y1_1 > 0]), fmt='.', label='y1')
    axs[5].set_title('y1')
    axs[5].set_xlabel('Time [s]')
    axs[5].set_ylabel('Counts')

    fig.tight_layout()
    fig.savefig(os.path.join(os.path.dirname(file_input), base_name+'_c1.pdf'))

    plt.close('all')

    #plt.show()

    return t_start, t_stop

def analyze_lc(lightcurve, p0=0.05, isRate=False, panels=['z0', 'z1', 'x0', 'x1', 'y0', 'y1']):
    """
    Analyze the light curve data.
    Input:
        - lightcurve: data file with times and counts
        - p0: false alarm probability for the Bayesian Block algorithm
        - isRate: if True, the lightcurve files contains rates instead of counts
        - panels: list of detectors for which we have lightcurves
    """

    data = np.genfromtxt(lightcurve, unpack=True)
    time = data[0]
    bin_width = time[1] - time[0]
    signal = {}
    for i, panel in enumerate(panels):
        signal[panel] = data[i+1]
        if isRate: signal[panel] = signal[panel] * bin_width

    # Construct light curve object
    bin_width = time[1] - time[0]
    bin_edges = np.zeros(len(time) + 1)
    bin_edges[1:-1] = (time[:-1] + time[1:]) / 2
    bin_edges[0] = time[0] - (time[1] - time[0]) / 2
    bin_edges[-1] = time[-1] + (time[-1] - time[-2]) / 2
    lo_edges, hi_edges = bin_edges[:-1], bin_edges[1:]
    exposure = np.full(len(time), bin_width)


    lc = {} # lc per panel
    for panel in panels:
        signal_panel = signal[panel]
        lc[panel] = TimeBins(signal_panel, lo_edges, hi_edges, exposure)
    if len(panels) > 1: lc_psum = TimeBins.sum([lcs for lcs in lc.values()])
    else: lc_psum = lc[panels[0]]


    # Apply Bayesian Blocks algorithm
    lc_sel = lc_psum
    try: 
        bb_lc = BayesianBlocksLightcurve(lc_sel)
        bb_lc.compute_bayesian_blocks(p0=p0)
        signal_range = bb_lc.signal_range
        t90 = bb_lc.duration(quantile = .9)
        t90_error = bb_lc.duration_error(.9, nsamples = 100)
    except: 
        print('WARNING: ')
        return lc_sel, None, -9999, -9999, -9999, -9999, -9999, -9999, -9999, None, None
    fig,ax = plt.subplots()
    fig2,ax2 = plt.subplots()
    ax = bb_lc.plot(ax=ax)
    ax2 = bb_lc.plot(ax=ax2)

    # Li&Ma calculation of the significance
    signal_lc = lc_sel.slice(bb_lc.signal_range.tstart, bb_lc.signal_range.tstop)
    bkg_lc = lc_sel.slice(bb_lc.signal_range.tstop, lc_sel.centroids[-1])
    bkg_lc = lc_sel.slice(lc_sel.centroids[0], bb_lc.signal_range.tstart)
    t_on = np.sum(signal_lc.exposure)
    t_off = np.sum(bkg_lc.exposure)
    alpha = t_on / t_off
    N_on = np.sum(signal_lc.rates * signal_lc.exposure)
    N_off = np.sum(bkg_lc.rates * bkg_lc.exposure)
    S = np.sqrt(2) * ( N_on * np.log( ((1+alpha)/alpha) * (N_on/(N_on+N_off)) ) + N_off * np.log( (1+alpha) * (N_off/(N_on+N_off)) ) ) ** 0.5

    # Li&Ma calculation of the peak significance
    significance = []
    for rate, exp in zip(signal_lc.rates, signal_lc.exposure):
        N_on = rate * exp
        t_on = exp
        alpha = t_on / t_off
        S_bin = np.sqrt(2) * ( N_on * np.log( ((1+alpha)/alpha) * (N_on/(N_on+N_off)) ) + N_off * np.log( (1+alpha) * (N_off/(N_on+N_off)) ) ) ** 0.5
        significance.append(S_bin)
    significance = np.array(significance)
    S_peak = np.max(significance)

    result = (lc_sel, bb_lc, signal_range.tstart, signal_range.tstop, t90, t90_error[0], t90_error[1], S, S_peak, ax, ax2)
    plt.close(fig)
    plt.close(fig2)

    # Save results
    GRB_dir = os.path.abspath(os.path.dirname(lightcurve))
    GRBname = os.path.basename(lightcurve).split('_')[0]
    if '_original' in lightcurve: mode = "original"
    elif '_megalib' in lightcurve: mode = "megalib"
    elif '_correction' in lightcurve: mode = "correction"
    else: mode = None
    save_results(GRB_dir, GRBname, signal_range.tstart, signal_range.tstop, t90, t90_error[0], t90_error[1], S, S_peak, mode)
    save_figure(ax, ax2, GRB_dir, GRBname, signal_range.tstart, signal_range.tstop, mode)

    return result

def save_results(GRB_dir, GRBname, signal_range1, signal_range2, t90, t90_error1, t90_error2, S, S_peak, Mode):
    """
    Save the results of the light curve analysis.
    """

    if Mode == "original": mode = "_original"
    elif Mode == "megalib": mode = "_megalib"
    elif Mode == "correction": mode = "_correction"
    else: mode = ""

    file_results_GRB = os.path.join(GRB_dir, f"bb_results"+mode+".dat")

    # Get GBM T90, flux and direction
    log_file = [f for f in os.listdir(GRB_dir) if f.endswith('.log')][0]
    source_file = [f for f in os.listdir(GRB_dir) if f.endswith('.source')][0]
    with open(log_file, "r") as f:
        for line in f:
            if "INFO: T90:" in line:
                items = line.split()
                t90_gbm = float(items[items.index("T90:") + 1])
    with open(source_file, 'r') as f:
        for line in f:
            if '.Flux' in line:
                line = line.strip()
                flux = float(line.split()[-1])
            if '.Beam' in line:
                line = line.strip()
                theta, phi = float(line.split()[-2]), float(line.split()[-1])
    
    with open(file_results_GRB, 'w') as f:
        f.write("# GRB_name flux[ph/s/cm2] theta[deg] phi[deg] t90_gbm[s] t_start[s] t_stop[s] t90[s] t90_error_low[s] t90_error_up[s] S S_peak\n")
        f.write(f"{GRBname} {flux} {theta} {phi} {t90_gbm} {signal_range1} {signal_range2} {t90} {t90_error1} {t90_error2} {S} {S_peak}\n")

def save_figure(ax, axzoom, GRB_dir, GRBname, signal_range1, signal_range2, Mode):
    """
    Save the figure of the light curve analysis.
    """

    if Mode == "original": 
        mode = "_original"
        title = "(Original)"
    elif Mode == "megalib":
        mode = "_megalib"
        title = "(MEGAlib correction)"
    elif Mode == "correction":
        mode = "_correction"
        title = "(Correction matrix)"
    else:
        mode = ""
        title = ""

    ax.set_title(f"{GRBname} {title}")
    axzoom.set_title(f"{GRBname} {title}")
    axzoom.set_xlim(signal_range1-5, signal_range2+5)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Rate [Hz]")
    axzoom.set_xlabel("Time [s]")
    axzoom.set_ylabel("Rate [Hz]")
    ax.legend()
    axzoom.legend()
    ax.grid(alpha=0.5)
    axzoom.grid(alpha=0.5)
    fig = ax.figure
    figzoom = axzoom.figure
    fig.tight_layout()
    figzoom.tight_layout()

    fig.savefig(os.path.join(GRB_dir, "bb_results"+mode+".pdf"))
    figzoom.savefig(os.path.join(GRB_dir, "bb_results_zoom"+mode+".pdf"))

    plt.close(fig)
    plt.close(figzoom)

