from odbAccess import *
import numpy as np
import os
from . import Constants as C
import re
from itertools import zip_longest
import csv
import time


def PostProcess(jobName, fileName, materialParameters):
    odbName = jobName + ".odb"
    odb = openOdb(path=odbName, readOnly=True)

    #### ---------------- ####
    #    Get displacements
    #### ---------------- ####
    # contactRegionName = "s_Surf-1"
    # surface_names = odb.rootAssembly.surfaces.keys()
    # surface_region_name = next((s for s in surface_names if "S_" in s), None)

    outputFileName = fileName + "_Results.csv"
    outputFolder = "SimDataOutputs"

    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

    outputFilePath = os.path.join(outputFolder, outputFileName)

    allNodesInContactRegion = odb.rootAssembly.surfaces[
        C.slave_surface_name.upper()
    ].nodes[0]
    unique_nodes = {node.label: node for node in allNodesInContactRegion}.values()

    undeformedCoords = [
        (
            node.label,
            node.coordinates[0],
            node.coordinates[1] - C.ys2,
            node.coordinates[2] - C.dpo_z,
        )
        for node in unique_nodes
    ]

    sortedUndeformedCoords = sorted(
        undeformedCoords, key=lambda coord: (coord[3], -coord[1], coord[2])
    )

    displacementField = odb.steps.values()[-1].frames[-1].fieldOutputs["U"]

    dispSubset = displacementField.getSubset(
        region=odb.rootAssembly.nodeSets[C.contact_region_nodes_name.upper()]
    )
    displacements = {
        value.nodeLabel: np.array(value.data) for value in dispSubset.values
    }

    deformedCoords = []
    for nodeLabel, x, y, z in sortedUndeformedCoords:
        disp = displacements.get(nodeLabel, np.array([0.0, 0.0, 0.0]))
        deformedCoords.append((nodeLabel, x + disp[0], y + disp[1], z + disp[2]))

    #### ---------------- ####
    # Get forces and energies
    #### ---------------- ####

    steps = odb.steps.keys()
    history_regions = odb.steps[steps[0]].historyRegions.keys()
    history_indenter_region_name = next(
        (s for s in history_regions if "Node" in s), None
    )
    history_substrate_region_name = next(
        (s for s in history_regions if "Element" in s), None
    )

    def get_history_data(step_name, region_name):
        history_region = odb.steps[step_name].historyRegions[region_name]
        data = []
        time = []
        for key in history_region.historyOutputs.keys():
            output_data = np.array(history_region.historyOutputs[key].data).T
            if key == "RF1":
                time.append(output_data[0, :])
            data.append(output_data[1, :])
        return np.array(time), np.array(data)

    time_array, [rf1, rf2, rf3] = get_history_data(
        steps[0], history_indenter_region_name
    )

    _, [IE, KE] = get_history_data(steps[0], history_substrate_region_name)

    # get wallclock time
    wallclock_time = None
    sta_file = jobName + ".sta"
    with open(sta_file, "r") as f:
        for line in f:
            if "WALLCLOCK TIME" in line:
                # Extract numeric value (last number in the line)
                match = re.search(r"([\d\.]+)\s*$", line)
                if match:
                    wallclock_time = float(match.group(1))
                    break
    if wallclock_time is None:
        raise ValueError(f"WALLCLOCK TIME not found in {sta_file}")

    with open(outputFilePath, "w") as f:
        writer = csv.writer(f)
        f.write("# Simulated using Abaqus 2025 licenced by Aarhus Uniersity\n")
        f.write("# Made by Peter Thorhauge Moellmann\n")
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        f.write(f"# Simulation date and time: {time_str}\n")
        f.write("# ----------------------------\n")
        f.write(
            f"# Indenter type: Rockwell with tip radius {C.tip_radius}mm and cone angle {C.cone_angle} degrees\n"
        )
        mat_str = ", ".join([f"{k}={v}" for k, v in materialParameters.items()])
        f.write(f"# Material parameters: {mat_str}\n")
        f.write(f"# WallclockTime={wallclock_time:.2f} s\n")

        writer.writerow(
            [
                "Time",
                "RF1",
                "RF2",
                "RF3",
                "IE",
                "KE",
                "NodeLabel",
                "x_undeformed",
                "y_undeformed",
                "z_undeformed",
                "x_deformed",
                "y_deformed",
                "z_deformed",
            ]
        )

        nodeLabels = []
        xUndeformed = []
        yUndeformed = []
        zUndeformed = []
        xDeformed = []
        yDeformed = []
        zDeformed = []
        for (nodeLabel, x_und, y_und, z_und), (_, x_def, y_def, z_def) in zip(
            sortedUndeformedCoords, deformedCoords
        ):
            nodeLabels.append(nodeLabel)
            xUndeformed.append(x_und)
            yUndeformed.append(y_und)
            zUndeformed.append(z_und)
            xDeformed.append(x_def)
            yDeformed.append(y_def)
            zDeformed.append(z_def)

        rows = zip_longest(
            time_array.reshape(
                -1,
            ),
            rf1,
            rf2,
            rf3,
            IE,
            KE,
            nodeLabels,
            xUndeformed,
            yUndeformed,
            zUndeformed,
            xDeformed,
            yDeformed,
            zDeformed,
            fillvalue="",
        )

        writer.writerows(rows)

    print(f"CSV results written: {outputFilePath}")

    odb.close()
