## Note that this script has been only tested on SARS-CoV-2 spike proteins. ##

from pymol import cmd
from pymol import stored
import pandas as pd
import numpy as np
from datetime import datetime
import time
from argparse import ArgumentParser

def main(pdbs_file):
    wild, d614g, alpha, beta, gamma, delta, kappa, omicron, unknown_variant, potential_ace2_bound, not_trimeric, missing_chs = [],[],[],[],[],[],[],[],[],[],[],[]
    with open(pdbs_file, 'r') as f:
        for l in f:
            pdbs = l.split(',')
            
    # start analysis & save logs to .txt file.
    with open(f"{datetime}_PDB_preprocess_log.txt", "a") as f:
        for code in pdbs:
            print(f'\nLoading {code} in PyMOL...')
            print(f'\nLoading {code} in PyMOL...', file=f)
            cmd.load(f'{code}.cif', code)
            print(f'\nCheck if the structure contains chain A/B/C...')
            print(f'\nCheck if the structure contains chain A/B/C...', file=f)
            chs = cmd.get_chains(code)
            print(f'\n{code} has {",".join(chs)} chains.')
            print(f'\n{code} has {",".join(chs)} chains.', file=f)
            if 'A' in chs and 'B' in chs and 'C' in chs:
                print(f'\nCheck if the structure is trimeric & contains S1/S2 domains...')
                print(f'\nCheck if the structure is trimeric & contains S1/S2 domains...', file=f)
                resi614_xyz = cmd.get_coords(f'{code} and resi 614 and name CA', 1) # for S1
                resi797_xyz = cmd.get_coords(f'{code} and resi 797 and name CA', 1) # for S2 (in FP)
                print(f'\nresi614_xyz:{resi614_xyz}')
                print(f'\nresi797_xyz:{resi614_xyz}')
                print(f'\nresi614_xyz:{resi614_xyz}', file=f)
                print(f'\nresi797_xyz:{resi614_xyz}', file=f)
                if resi614_xyz is not None and resi797_xyz is not None and len(resi614_xyz) == 3 and len(resi797_xyz) == 3:
                    print('\nTrimeric. Continue to determine the variant types...')
                    print('\nTrimeric. Continue to determine the variant types...', file=f)
                    stored.resn = []
                    cmd.iterate(f'{code} and chain A and resi 614 and name CA', 'stored.resn.append(resn)')
        #             print(f'\nstored.resn614_chA:{stored.resn}\n)
                    if len(stored.resn) != 0 and stored.resn[0] == 'ASP': # wild-type
                        wild.append(code)
                    elif len(stored.resn) == 0:
                        print('Residue 614 is not present in chain A. Move on to the next structure.')
                        print('Residue 614 is not present in chain A. Move on to the next structure.', file=f)
                        continue
                    else: # could be d614g, alpha-kappa variants.
                        # get all the needed residue names for later analysis.
                        stored.resn611, stored.resn612, stored.resn614\
                        , stored.resn498, stored.resn499, stored.resn501\
                        , stored.resn19, stored.resn20, stored.resn26, stored.resn138, stored.resn190, stored.resn80, stored.resn215, stored.resn156\
                        , stored.resn154, stored.resn567, stored.resn678, stored.resn67, stored.resn339 = [],[],[],[],[]\
                        , [],[],[],[],[]\
                        , [],[],[],[],[]\
                        , [],[],[],[]

                        cmd.iterate(f'{code} and chain A and resi 611 and name CA', 'stored.resn611.append(resn)')
                        print(f'\nResidue name of resi611:{stored.resn611}')
                        print(f'\nResidue name of resi611:{stored.resn611}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 612 and name CA', 'stored.resn612.append(resn)')
                        print(f'\nResidue name of resi612:{stored.resn612}')
                        print(f'\nResidue name of resi612:{stored.resn612}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 614 and name CA', 'stored.resn614.append(resn)')
                        print(f'\nResidue name of resi614:{stored.resn614}')
                        print(f'\nResidue name of resi614:{stored.resn614}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 498 and name CA', 'stored.resn498.append(resn)')
                        print(f'\nResidue name of resi498:{stored.resn498}')
                        print(f'\nResidue name of resi498:{stored.resn498}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 499 and name CA', 'stored.resn499.append(resn)')
                        print(f'\nResidue name of resi499:{stored.resn499}')
                        print(f'\nResidue name of resi499:{stored.resn499}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 501 and name CA', 'stored.resn501.append(resn)')
                        print(f'\nResidue name of resi501:{stored.resn501}')
                        print(f'\nResidue name of resi501:{stored.resn501}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 19 and name CA', 'stored.resn19.append(resn)')
                        print(f'\nResidue name of resi19:{stored.resn19}')
                        print(f'\nResidue name of resi19:{stored.resn19}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 20 and name CA', 'stored.resn20.append(resn)')
                        print(f'\nResidue name of resi20:{stored.resn20}')
                        print(f'\nResidue name of resi20:{stored.resn20}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 26 and name CA', 'stored.resn26.append(resn)')
                        print(f'\nResidue name of resi26:{stored.resn26}')
                        print(f'\nResidue name of resi26:{stored.resn26}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 138 and name CA', 'stored.resn138.append(resn)')
                        print(f'\nResidue name of resi138:{stored.resn138}')
                        print(f'\nResidue name of resi138:{stored.resn138}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 190 and name CA', 'stored.resn190.append(resn)')
                        print(f'\nResidue name of resi190:{stored.resn190}')
                        print(f'\nResidue name of resi190:{stored.resn190}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 80 and name CA', 'stored.resn80.append(resn)')
                        print(f'\nResidue name of resi80:{stored.resn80}')
                        print(f'\nResidue name of resi80:{stored.resn80}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 215 and name CA', 'stored.resn215.append(resn)')
                        print(f'\nResidue name of resi215:{stored.resn215}')
                        print(f'\nResidue name of resi215:{stored.resn215}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 154 and name CA', 'stored.resn154.append(resn)')
                        print(f'\nResidue name of resi154:{stored.resn154}')
                        print(f'\nResidue name of resi154:{stored.resn154}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 156 and name CA', 'stored.resn156.append(resn)')
                        print(f'\nResidue name of resi156:{stored.resn156}')
                        print(f'\nResidue name of resi156:{stored.resn156}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 567 and name CA', 'stored.resn567.append(resn)')
                        print(f'\nResidue name of resi567:{stored.resn567}')
                        print(f'\nResidue name of resi567:{stored.resn567}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 678 and name CA', 'stored.resn678.append(resn)')
                        print(f'\nResidue name of resi678:{stored.resn678}')
                        print(f'\nResidue name of resi678:{stored.resn678}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 67 and name CA', 'stored.resn67.append(resn)')
                        print(f'\nResidue name of resi67:{stored.resn67}')
                        print(f'\nResidue name of resi67:{stored.resn67}', file=f)
                        cmd.iterate(f'{code} and chain A and resi 339 and name CA', 'stored.resn339.append(resn)')
                        print(f'\nResidue name of resi339:{stored.resn339}')
                        print(f'\nResidue name of resi339:{stored.resn339}', file=f)

                        if len(stored.resn611)!=0 and len(stored.resn498)!=0 and len(stored.resn567)!=0 and len(stored.resn678)!=0 and stored.resn611[0] == 'GLY' and stored.resn498[0] == 'TYR' and stored.resn567[0] == 'ASP' and stored.resn678[0] == 'HIS':
                            alpha.append(code)
                        elif len(stored.resn611)!=0 and len(stored.resn498)!=0 and len(stored.resn80)!=0 and len(stored.resn215)!=0 and stored.resn611[0] == 'GLY' and stored.resn498[0] == 'TYR' and stored.resn80[0] == 'ALA' and stored.resn215[0] == 'GLY':
                            beta.append(code)
                        elif len(stored.resn614)!=0 and len(stored.resn501)!=0 and len(stored.resn20)!=0 and len(stored.resn26)!=0 and len(stored.resn138)!=0 and len(stored.resn190)!=0 and stored.resn614[0] == 'GLY' and stored.resn501[0] == 'TYR' and stored.resn20[0] == 'ASN' and stored.resn26[0] == 'SER' \
                            and stored.resn138[0] == 'TYR' and stored.resn190[0] == 'SER':
                            gamma.append(code)
                        elif len(stored.resn612)!=0 and len(stored.resn499)!=0 and len(stored.resn19)!=0 and len(stored.resn156)!=0 and  stored.resn612[0] == 'GLY' and stored.resn499[0] == 'ASN' and stored.resn19[0] == 'ARG' and stored.resn156[0] == 'GLY':
                            delta.append(code)
                        elif len(stored.resn614)!=0 and len(stored.resn501)!=0 and len(stored.resn154)!=0 and stored.resn614[0] == 'GLY' and stored.resn501[0] == 'ASN' and stored.resn154[0] == 'LYS':
                            kappa.append(code)
                        elif len(stored.resn614)!=0 and len(stored.resn501)!=0 and len(stored.resn67)!=0 and len(stored.resn339)!=0 and stored.resn614[0] == 'GLY' and stored.resn501[0] == 'TYR' and stored.resn67[0] == 'VAL' and stored.resn339[0] == 'ASP':
                            omicron.append(code)
                        elif len(stored.resn614)!=0 and len(stored.resn501)!=0 and len(stored.resn154)!=0 and stored.resn614[0] == 'GLY' and stored.resn501[0] == 'ASN' and stored.resn154[0] == 'GLU':
                            d614g.append(code)
                        else:
                            print('\nCannot determine the variant type.')
                            print('\nCannot determine the variant type.', file=f)
                            unknown_variant.append(code)
                elif resi614_xyz is not None and resi797_xyz is not None and len(resi614_xyz) > 3 and len(resi797_xyz) == 3: # since ACE2 also has resi 614
                    print('\nPotential ACE2-bound structure.')
                    print('\nPotential ACE2-bound structure.', file=f)
                    potential_ace2_bound.append(code)
                else:
                    print('\nNot trimeric. Move on to the next structure.')
                    print('\nNot trimeric. Move on to the next structure.', file=f)
                    not_trimeric.append(code)
            else:
                print('\nA/B/C chains are not present simultaneously. Move on to the next structure.')
                print('\nA/B/C chains are not present simultaneously. Move on to the next structure.', file=f)
                missing_chs.append(code)

        print(f"\n<Summary>\nWild({len(wild)}):\n{','.join(wild)}\nD614G({len(d614g)}):\n{','.join(d614g)}\nAlpha({len(alpha)}):\n{','.join(alpha)}\nBeta({len(beta)}):\n{','.join(beta)}\nGamma({len(gamma)}):\n{','.join(gamma)}\nDelta({len(delta)}):\n{','.join(delta)}\nKappa({len(kappa)}):\n{','.join(kappa)}\nOmicron({len(omicron)}):\n{','.join(omicron)}\nUnknown Variant({len(unknown_variant)}):\n{','.join(unknown_variant)}\nPotential ACE2-bound({len(potential_ace2_bound)}):\n{','.join(potential_ace2_bound)}\nNot trimeric({len(not_trimeric)}):\n{','.join(not_trimeric)}\nMissing A/B/C chains({len(missing_chs)}):\n{','.join(missing_chs)}")                   
        print(f"\n<Summary>\nWild({len(wild)}):\n{','.join(wild)}\nD614G({len(d614g)}):\n{','.join(d614g)}\nAlpha({len(alpha)}):\n{','.join(alpha)}\nBeta({len(beta)}):\n{','.join(beta)}\nGamma({len(gamma)}):\n{','.join(gamma)}\nDelta({len(delta)}):\n{','.join(delta)}\nKappa({len(kappa)}):\n{','.join(kappa)}\nOmicron({len(omicron)}):\n{','.join(omicron)}\nUnknown Variant({len(unknown_variant)}):\n{','.join(unknown_variant)}\nPotential ACE2-bound({len(potential_ace2_bound)}):\n{','.join(potential_ace2_bound)}\nNot trimeric({len(not_trimeric)}):\n{','.join(not_trimeric)}\nMissing A/B/C chains({len(missing_chs)}):\n{','.join(missing_chs)}", file=f)
        
        print('\n-----Elapsed time: %.2f seconds.-----\n'%(time.time()-start))
        print('\n-----Elapsed time: %.2f seconds.-----\n'%(time.time()-start), file=f)

if __name__=="__main__":
    parser=ArgumentParser()
    parser.add_argument('--PDBids', dest='PDBids', help='The directory to a text file containing all the PDB IDs for the analysis. Comma-separated format. No spaces.')
    
    # extract the args values.
    args = parser.parse_args()
    pdbs_file = args.PDBids 
    
    # execute the main func.
    start = time.time()
    datetime = datetime.now().strftime('%Y%m%d-%H%M%S')
    main(pdbs_file)