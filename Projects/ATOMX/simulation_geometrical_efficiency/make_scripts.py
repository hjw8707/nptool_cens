import sys

sys.path.append("../common")
from atomx_simulation_script_writter import atomx_simulation_script_writter

if __name__ == "__main__":
    #print(atomx_simulation_script_writter.__doc__)

    a = atomx_simulation_script_writter("ATOMX_geometrical_efficiency")
    a.set_path("/mnt/nptool/atomx/data", "input", "/mnt/nptool/atomx/figures")
    a.make_project_configuration()
    a.make_geant4_vis_file()

    nbeams = 100000
    nz = 20
    z1 = -140
    z2 = +140
    dz = (z2-z1)/nz
    sz = dz/6.
    all_batch = open("all_batch.sh","w")
    for iz in range(nz):
        print()
        z0 = z1 + dz*iz
        name = f"ATOMX_SiArray_{iz}_z{z0}"
        detector = a.make_detector_file(make_gas=False)
        reaction = a.make_eventgen_file_isotropic(particle_name="proton", energy1=15, energy2=15, angle1=0, angle2=180, z0=z0, sz=sz, eventgen_name=f"ATOMX_{iz}_z{z0}")
        batch_file, viewer_file = a.make_script(detector.name, reaction.name, run_name=name, macro_for_ana="root -q -b draw_hit.C", nbeams=nbeams)
        print(f"sh {batch_file.name}",file=all_batch)
    print()
    print("root -q -b draw_summary.C", file=all_batch)
    print(all_batch.name)
