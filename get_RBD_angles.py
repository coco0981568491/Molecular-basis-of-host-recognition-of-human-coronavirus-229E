## Note that this script has been only tested on SARS-CoV-2 spike proteins. ##

from pymol import cmd
import numpy as np
import pandas as pd
from numpy import linalg as la
import math
from datetime import datetime
import time
# import smtplib
# from email.mime.multipart import MIMEMultipart
# from email.mime.text import MIMEText
# from email.mime.base import MIMEBase
# from email import encoders
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

# Record start time & date.
start = time.time()
datetime = datetime.now().strftime('%Y%m%d-%H%M%S')

# Delete existing objects from PDB_preprocessor.py analysis.
cmd.delete("all")
print("Delete all existing objects within PyMOL from PDB_preprocessor.py analysis...\n")

# info. for final table.
# 10 files below rbd_start_resi_absent:['7edj', '7v7z', '7v81', '7v82', '7v83', '7v88', '7v89', '7v8a', '7v85', '7v86']
# 6wps has C chain for Ab. 6xcn B chain for Ab. 7jv6. 7jvc
# these files were already removed from pdb_class_dict.

# Please use "FindAb" func. in RBD_Angle_PostProcessing notebook to complete the following groups.
# pdb_class_dict = {'wild_apo':['6vxx', '6vyb', '6x29', '6x2a', '6x2b', '6x2c', '6x6p', '6x79', '6xf5', '6xf6', '6xkl', '6xlu', '6xm0', '6xm3', '6xm4', '6xm5', '6xr8', '6z97', '6zb4', '6zb5', '6zge', '6zgg', '6zgi', '6zow', '6zox', '6zoy', '6zoz', '6zp0', '6zp1', '6zp2', '6zp7', '7a93', '7cab', '7cn9', '7ddd', '7ddn', '7df3', '7dk3', '7dwy', '7dwz', '7dx0', '7dx3', '7e7b', '7e7d', '7fcd', '7fce', '7jji', '7jwy', '7kdg', '7kdh', '7kj5', '7l7k', '7m0j', '7mjg', '7mtc', '7mte', '7nt9', '7nta', '7od3', '7odl', '7p7b', '7tla', '7tlb', '7tlc', '7tld']\
#                     , 'wild_ace2':[]\
#                     , 'wild_ab':['6wpt', '6xcm', '6z43', '6zdh', '6zhd', '6zxn', '7a25', '7a29', '7akd', '7b18', '7c2l', '7cac', '7cai', '7cak', '7cwl', '7cwn', '7cyp', '7d0b', '7d0c', '7d0d', '7dk4', '7dk5', '7dk6', '7dk7', '7e3k', '7e3l', '7e5r', '7e5s', '7e8c', '7ej4', '7ej5', '7fae', '7faf', '7jv4', '7jw0', '7jzl', '7jzn', '7k43', '7k4n', '7k8s', '7k8t', '7k8u', '7k8v', '7k8w', '7k8x', '7k8z', '7k90', '7k9h', '7k9j', '7kkk', '7kkl', '7kl9', '7kmk', '7kml', '7kqb', '7kqe', '7ks9', '7ksg', '7kxj', '7kxk', '7l02', '7l06', '7l09', '7l2d', '7l2e', '7l2f', '7l3n', '7l56', '7l57', '7l58', '7laa', '7lab', '7lcn', '7ld1', '7ljr', '7lqv', '7lqw', '7lrt', '7ls9', '7lxy', '7lxz', '7ly2', '7m6e', '7m6f', '7m6g', '7m6h', '7m6i', '7mjh', '7mjj', '7mjk', '7mkl', '7mm0', '7my3', '7n0g', '7n0h', '7n5h', '7n8h', '7n9e', '7n9t', '7nd3', '7nd4', '7nd5', '7nd6', '7nd7', '7nd8', '7nd9', '7nda', '7ndb', '7ndc', '7ndd', '7ntc', '7oan', '7p40', '7p78', '7p79', '7r8m', '7r8n', '7r8o', '7ra8', '7rkv', '7rw2', '7s0c', '7s0d', '7sc1', '7sn3', '7tat', '7v2a', '7vnc', '7vnd', '7vne']\
#                     , 'd614g_apo':['6zwv', '7bnm', '7bnn', '7bno', '7dx1', '7dx2', '7dzw', '7eaz', '7eb0', '7eb3', '7eb4', '7eb5', '7kdi', '7kdj', '7kdk', '7kdl', '7ke4', '7ke6', '7ke7', '7ke8', '7ke9', '7kea', '7keb', '7kec', '7krq', '7krr', '7krs', '7lwi', '7lwj', '7lwk', '7lwl', '7lwm', '7lwn', '7lwo', '7lwp', '7sbk', '7sbl', '7sbo', '7sxr', '7sxs']\
#                     , 'd614g_ace2':[]\
#                     , 'd614g_ab':['7dzx', '7dzy', '7eh5', '7so9']\
#                     , 'alpha_apo':['7edf', '7edg', '7edi']\
#                     , 'alpha_ace2':['7edj']\
#                     , 'alpha_ab':[]\
#                     , 'beta_apo':['7v76', '7v77', '7v8c']\
#                     , 'beta_ace2':['7v7z']\
#                     , 'beta_ab':['7q6e', '7q9f', '7q9g', '7q9i', '7q9j', '7q9k', '7q9m', '7q9p']\
#                     , 'gamma_apo':['7m8k', '7sbs', '7sbt', '7v78', '7v79', '7v7a']\
#                     , 'gamma_ace2':['7v81', '7v82', '7v83']\
#                     , 'gamma_ab':[]\
#                     , 'delta_apo':['7v7n', '7v7o', '7v7p', '7v7q', '7v7r', '7v7s', '7v7t', '7v7u', '7v7v', '7w92', '7w94']\
#                     , 'delta_ace2':['7v88', '7v89', '7v8a']\
#                     , 'delta_ab':['7w9e']\
#                     , 'kappa_apo':['7sbp', '7sbq', '7sbr', '7v7d', '7v7e', '7v7f', '7v7g', '7vxe', '7vxi']\
#                     , 'kappa_ace2':['7v85', '7v86']\
#                     , 'kappa_ab':['7sob', '7soe']\
#                     , 'omicron_apo':['7t9j', '7tb4', '7wk2', '7wk3']\
#                     , 'omicron_ace2':[]\
#                     , 'omicron_ab':['7tm0', '7wk9', '7wka']}

pdb_class_dict = {'wild_apo':['6vsb', '6vxx', '6vyb', '6x29', '6x2a', '6x2b', '6x2c', '6x6p', '6x79', '6xf5', '6xf6', '6xkl', '6xlu', '6xm0', '6xm3', '6xm4', '6xm5', '6xr8', '6z97', '6zb4', '6zb5', '6zge', '6zgg', '6zgh', '6zgi', '6zow', '6zox', '6zoy', '6zoz', '6zp0', '6zp1', '6zp2', '6zp7', '7a93', '7cab', '7cn9', '7ddd', '7ddn', '7df3', '7dk3', '7dwy', '7dwz', '7dx0', '7dx3', '7e7b', '7e7d', '7fb0', '7fb1', '7fb3', '7fb4', '7fcd', '7fce', '7jji', '7jwy', '7kdg', '7kdh', '7kj5', '7l7k', '7m0j', '7mjg', '7mtc', '7mtd', '7mte', '7nt9', '7nta', '7od3', '7odl', '7p7b', '7q1z', '7qur', '7qus', '7ru1', '7ru2', '7s6i', '7tla', '7tlb', '7tlc', '7tld', '7v20', '7vhj', '7vhk', '7vrv', '7vrw', '7wgv', '7wgx', '7wgy', '7wgz', '7xtz', '7xu0', '7xu1', '7xu2', '7xu3', '7xu6', '7yqu', '7z3z', '8h3d', '8h3e','7n9b']\
                    , 'wild_ace2':['7df4','7dx5','7dx6','7dx7','7dx8','7dx9','7kms','7kmz','7knb','7kne','7knh','7kni']\
                    , 'wild_ab':['6wps', '6wpt', '6xcm', '6xcn', '6xey', '6z43', '6zdh', '6zhd', '6zxn', '7a25', '7a29', '7a94', '7a95', '7a96', '7a97', '7a98', '7akd', '7b18', '7byr', '7c2l', '7cac', '7cai', '7cak', '7chh', '7cwl', '7cwm', '7cwn', '7cwt', '7cwu', '7cyp', '7czp', '7czq', '7czr', '7czs', '7czt', '7czu', '7czv', '7czw', '7czx', '7czy', '7czz', '7d00', '7d03', '7d0b', '7d0c', '7d0d', '7dk4', '7dk5', '7dk6', '7dk7', '7e3k', '7e3l', '7e5r', '7e5s', '7e8c', '7e9n', '7e9o', '7e9q', '7ej4', '7ej5', '7enf', '7epx', '7fae', '7faf', '7jv4', '7jv6', '7jvc', '7jw0', '7jwb', '7jzl', '7jzn', '7k43', '7k4n', '7k8s', '7k8t', '7k8u', '7k8v', '7k8w', '7k8x', '7k8z', '7k90', '7k9h', '7k9j', '7kj2', '7kj3', '7kj4', '7kkk', '7kkl', '7kl9', '7kmk', '7kml', '7kqe', '7ks9', '7ksg', '7kxj', '7kxk', '7l02', '7l06', '7l09', '7l2d', '7l2e', '7l2f', '7l3n', '7l56', '7l57', '7l58', '7laa', '7lab', '7lcn', '7ld1', '7ljr', '7lqv', '7lqw', '7lrt', '7ls9', '7lss', '7lxy', '7lxz', '7ly2', '7m6e', '7m6f', '7m6g', '7m6h', '7m6i', '7mjh', '7mjj', '7mjk', '7mkl', '7mm0', '7mw2', '7mw3', '7mw4', '7mw5', '7mw6', '7my3', '7n0g', '7n0h', '7n5h', '7n8h', '7n9e', '7n9t', '7nd3', '7nd4', '7nd5', '7nd6', '7nd7', '7nd8', '7nd9', '7nda', '7ndb', '7ndc', '7ndd', '7ntc', '7ny5', '7oan', '7p40', '7p78', '7p79', '7q0a', '7r40', '7r4i', '7r4q', '7r4r', '7r8m', '7r8n', '7r8o', '7ra8', '7rbv', '7rkv', '7rq6', '7ru3', '7ru5', '7rw2', '7s0c', '7s0d', '7s6j', '7s6k', '7s6l', '7sc1', '7sn3', '7swx', '7tat', '7tb8', '7tpr', '7tyz', '7u0p', '7u0q', '7u0x', '7uap', '7uar', '7uhc', '7uz4', '7uz5', '7uz6', '7uz7', '7uz9', '7uza', '7v23', '7v2a', '7vnc', '7vnd', '7vne', '7whb', '7whd', '7ws0', '7ws1', '7x7n', '7xmx', '7xmz', '7y7j', '7ydy', '7ye5', '7ye9', '7z6v', '7z7x', '7z85', '7z86', '7z9q', '7z9r', '7zce', '8cxn', '8cxq', '8cy6', '8cy7', '8cy9', '8cya', '8cyb', '8cyc', '8cyd', '8d0z', '8dxs', '8hc4']\
                    , 'd614g_apo':['6xs6', '6zwv', '7bnm', '7bnn', '7bno', '7dx1', '7dx2', '7dzw', '7eaz', '7eb0', '7eb3', '7eb4', '7eb5', '7kdi', '7kdj', '7kdk', '7kdl', '7ke4', '7ke6', '7ke7', '7ke8', '7ke9', '7kea', '7keb', '7kec', '7krq', '7krr', '7krs', '7lwi', '7lwk', '7lwl', '7lwm', '7lwn', '7lwo', '7lwp', '7r18', '7r19', '7r1b', '7sbk', '7sbl', '7sbo', '7sxr', '7sxs', '7t67', '7tgx', '7tgy', '7w92', '7w94', '7xu4', '7xu5', '8dlt', '7lwj']\
                    , 'd614g_ace2':['7sxx']\
                    , 'd614g_ab':['7dzx', '7dzy', '7eh5', '7so9', '7t3m', '7vq0', '7w9e', '7wo5', '7woa', '7wob', '7ws3', '7wuh', '7zss', '8dlw', '8dlx', '8dlz', '8hhx']\
                    , 'alpha_apo':['7edf', '7edg', '7edh', '7edi','7fet','7lws','7lwt','7lwu','7lwv','7r13','7r14','7r15','8dli']\
                    , 'alpha_ace2':['7edj','8dlj','7fem','7r1a']\
                    , 'alpha_ab':[]\
                    , 'beta_apo':['7v76', '7v77', '7v8c','7lyk','7lyl','7lym','7lyn','7lyo','7lyp','7lyq','7r16','7r17','8dll','7vx1','7wev']\
                    , 'beta_ace2':['7v7z','7vxd','7vxf','7vxk','7vxm','8dlm']\
                    , 'beta_ab':['7q6e', '7q9f', '7q9g', '7q9i', '7q9j', '7q9k', '7q9m', '7q9p', '7zr7', '7zr9', '7zrc','7fjn','7fjo','7wcz','7wd0','7wd7','7wd9','7wdf']\
                    , 'gamma_apo':['7m8k', '7sbs', '7sbt', '7v78', '7v79', '7v7a', '8dlo']\
                    , 'gamma_ace2':['7v81','7v82','7v83','8dlp']\
                    , 'gamma_ab':[]\
                    , 'delta_apo':['7v7n', '7v7o', '7v7p', '7v7q', '7v7r', '7v7s', '7v7t', '7v7u', '7v7v', '7wg8', '7wg9','7tey','7tou','7tov','7tox','7toy','7toz','7toz','7tp0','7tp1','7tp2','7tp7','7tp8','7tp9','7tpa','7tpc','7tpe','7tpf','7tph','7vhh']\
                    , 'delta_ace2':['7tex','7v88','7v89','7v8a','7w99','7w9b','7w9c']\
                    , 'delta_ab':['8hhy','7wcp','7wwl','7wwm']\
                    , 'kappa_apo':['7sbp', '7sbq', '7sbr', '7tf3', '7v7d', '7v7e', '7v7f', '7v7g', '7vxe', '7vxi']\
                    , 'kappa_ace2':['7tf0','7v85','7v86','7vx9','7vxa','7vxb','7vxc']\
                    , 'kappa_ab':['7sob', '7soe']\
                    , 'omicron_apo':['7t9j', '7tb4', '7tei', '7tf8', '7tge', '7thk', '7tl1', '7tl9', '7tnw', '7to4', '7wg6', '7wg7', '7wk2', '7wk3', '7wz1', '7xiy', '7y9s','7qo7','8gs6','8dm1','7tgw','7ub0','7ub5','7ub6','7wp9','7xiw','7xix','7xnq','7xnr','7xns','7yqt','7yqv','7yqw']\
                    , 'omicron_ace2':['7wk4','7wk5','7wvp','7wvq','7xch','7xid','7y9z','7yr2','7yr3','8dm5','7wpa','7xo7','7xo8']\
                    , 'omicron_ab':['7t9k', '7tca', '7tcc', '7tm0', '7upl', '7uz8', '7wgb', '7whi', '7whj', '7whk', '7wjy', '7wk9', '7wka', '7woq', '7wor', '7wos', '7wou', '7wov', '7wpd', '7wpf', '7wrh', '7ws4', '7ws5', '7ws8', '7ws9', '7wwi', '7wwj', '7xco', '7xst', '8dm9', '8err', '8hc2', '8hc3', '8hc9', '8hca', '8hcb', '8hhz','7qti','8dzi','8dzh','8dm3','7we8','7we9','7wea','7wec','7wpe','7wti','7wtk','7xdb','7xic','7xod','7yqx','7yqy','7yqz','7yr1']}

classes = list(pdb_class_dict.keys())
chids = ['A', 'B', 'C']
abnormal_chid = []
dfs = []
sheet_names = []

# Start analysis & save logs to .txt file.
with open(f"{datetime}_get_RBD_angles_log.txt", "a") as f:
    print(f"\nclasses:\n{classes}\n")
    print(f"\nclasses:\n{classes}\n", file=f)
    for c in classes:
        if len(pdb_class_dict[c]) != 0:
            print(f'\nclass:{c}')
            print(f'\nclass pdbs:\n{pdb_class_dict[c]}\n')
            print(f'\nclass:{c}', file=f)
            print(f'\nclass pdbs:\n{pdb_class_dict[c]}\n', file=f)
            rbd_angles_chA, rbd_angles_chB, rbd_angles_chC = [], [], []
            pdbs = []
            sheet_names.append(c) # only write existing data for clarity.
            for code in pdb_class_dict[c]:
                print(f'\n[Protein: {code} ({c})...]\n')
                print(f'\n[Loading {code} in PyMOL...]\n')
                print(f'\n[Protein: {code} ({c})...]\n', file=f)
                print(f'\n[Loading {code} in PyMOL...]\n', file=f)
                objname = code
                cmd.load(f'{code}.cif', objname)
                print(f'\n[Start processing {objname}...]\n')
                print(f'\n[Start processing {objname}...]\n', file=f)
                cmd.remove(f'{objname} and het')

                # note that the resi330 may shift due to deletions in different variants.
                if 'alpha' in c or 'beta' in c:
                    rbd_start_resi = 327
                elif 'delta' in c:
                    rbd_start_resi = 328
                elif 'wild' in c or 'd614g' in c or 'gamma' in c or 'kappa' in c or 'omicron' in c:
                    rbd_start_resi = 330

                # make sure the A/B/C chains are only for spike, not for Ab.
                resi330_xyz_chA = cmd.get_coords(f'{objname} and resi {rbd_start_resi} and name CA and chain A', 1)
                resi330_xyz_chB = cmd.get_coords(f'{objname} and resi {rbd_start_resi} and name CA and chain B', 1)
                resi330_xyz_chC = cmd.get_coords(f'{objname} and resi {rbd_start_resi} and name CA and chain C', 1)
                resi520_xyz_chA = cmd.get_coords(f'{objname} and resi {rbd_start_resi+190} and name CA and chain A', 1)
                resi520_xyz_chB = cmd.get_coords(f'{objname} and resi {rbd_start_resi+190} and name CA and chain B', 1)
                resi520_xyz_chC = cmd.get_coords(f'{objname} and resi {rbd_start_resi+190} and name CA and chain C', 1)
                print(f'\nresi330_xyz_chA:{resi330_xyz_chA}')
                print(f'\nresi330_xyz_chB:{resi330_xyz_chB}')
                print(f'\nresi330_xyz_chC:{resi330_xyz_chC}')
                print(f'\nresi520_xyz_chA:{resi520_xyz_chA}')
                print(f'\nresi520_xyz_chB:{resi520_xyz_chB}')
                print(f'\nresi520_xyz_chC:{resi520_xyz_chC}')
                print(f'\nresi330_xyz_chA:{resi330_xyz_chA}', file=f)
                print(f'\nresi330_xyz_chB:{resi330_xyz_chB}', file=f)
                print(f'\nresi330_xyz_chC:{resi330_xyz_chC}', file=f)
                print(f'\nresi520_xyz_chA:{resi520_xyz_chA}', file=f)
                print(f'\nresi520_xyz_chB:{resi520_xyz_chB}', file=f)
                print(f'\nresi520_xyz_chC:{resi520_xyz_chC}', file=f)
                if resi330_xyz_chA is not None and resi330_xyz_chB is not None and resi330_xyz_chC is not None and resi520_xyz_chA is not None and resi520_xyz_chB is not None and resi520_xyz_chC is not None:
                    pdbs.append(code)
                    # get 2 vectors (zx & zy) to calculate the normal vector of the plane.
                    resi330_31 = resi330_xyz_chA - resi330_xyz_chC
                    resi330_32 = resi330_xyz_chB - resi330_xyz_chC
                    normal_v = np.cross(resi330_31, resi330_32) 
                    com_rbd_a = cmd.centerofmass(f'{objname} and resi {rbd_start_resi}-{rbd_start_resi+190} and chain A')
                    com_rbd_b = cmd.centerofmass(f'{objname} and resi {rbd_start_resi}-{rbd_start_resi+190} and chain B')
                    com_rbd_c = cmd.centerofmass(f'{objname} and resi {rbd_start_resi}-{rbd_start_resi+190} and chain C')
                    cmd.pseudoatom(f'{objname}_chA_rbd_com', pos=f'[{com_rbd_a[0]},{com_rbd_a[1]},{com_rbd_a[2]}]')
                    cmd.pseudoatom(f'{objname}_chB_rbd_com', pos=f'[{com_rbd_b[0]},{com_rbd_b[1]},{com_rbd_b[2]}]')
                    cmd.pseudoatom(f'{objname}_chC_rbd_com', pos=f'[{com_rbd_c[0]},{com_rbd_c[1]},{com_rbd_c[2]}]')
                    print(f'{objname}_com_rbd_a:{com_rbd_a}')
                    print(f'{objname}_com_rbd_b:{com_rbd_b}')
                    print(f'{objname}_com_rbd_c:{com_rbd_c}')
                    print(f'{objname}_com_rbd_a:{com_rbd_a}', file=f)
                    print(f'{objname}_com_rbd_b:{com_rbd_b}', file=f)
                    print(f'{objname}_com_rbd_c:{com_rbd_c}', file=f)

                    # calculate all the rbd vectors to get the angles btw vector and the normal vector using resi 330 as the origin.
                    ch_cnt = 1
                    for ch in chids:
                        ori = cmd.get_coords(f'{objname} and resi {rbd_start_resi} and name CA and chain {ch}', 1)[0]
                        cmd.pseudoatom(f'{objname}_rbd_ori_{ch_cnt}', pos=f'[{ori[0]},{ori[1]},{ori[2]}]')
                        rbd_com = cmd.centerofmass(f'{objname} and resi {rbd_start_resi}-{rbd_start_resi+190} and chain {ch}')
                        rbd_v = rbd_com - ori
                        angle = np.around(math.degrees(np.arccos(np.inner(rbd_v, normal_v)/(la.norm(rbd_v)*la.norm(normal_v)))), decimals=1) 
                        if angle > 90:
                            angle = 180 - angle

                        if ch == 'A':
                            rbd_angles_chA.append(angle)
                        elif ch == 'B':
                            rbd_angles_chB.append(angle)
                        elif ch == 'C':
                            rbd_angles_chC.append(angle)
                        ch_cnt += 1
                else:
                    print('Not all A/B/C chains are for spike. Move on to the next structure.')
                    print('Not all A/B/C chains are for spike. Move on to the next structure.', file=f)
                    pdbs.append(code)
                    rbd_angles_chA.append(np.nan)
                    rbd_angles_chB.append(np.nan)
                    rbd_angles_chC.append(np.nan)
                    abnormal_chid.append(code)

            # making the summary table for each qualified data.
            df = pd.DataFrame({'PDB codes':pdbs, '[chainA]RBD angle':rbd_angles_chA, '[chainB]RBD angle':rbd_angles_chB, '[chainC]RBD angle':rbd_angles_chC})
            dfs.append(df)

        else:
            print(f'\nclass:{c}')
            print(f'\nclass pdbs:\n{pdb_class_dict[c]}\n')
            print(f'\nEmpty list. Move on to the next class.\n')
            print(f'\nclass:{c}', file=f)
            print(f'\nclass pdbs:\n{pdb_class_dict[c]}\n', file=f)
            print(f'\nEmpty list. Move on to the next class.\n', file=f)

    # final summary excel file.
    c_cnt = 0
    with pd.ExcelWriter(f'{datetime}_RBD_Angles.xlsx') as writer:
        for df in dfs:
            df.to_excel(writer, sheet_name = f'{sheet_names[c_cnt]}', index = False)
            c_cnt += 1

    print(f'{len(abnormal_chid)} files abnormal_chid:\n{abnormal_chid}\n')
    print(f'{len(abnormal_chid)} files abnormal_chid:\n{abnormal_chid}\n', file=f)

    # Elapsed time.
    print('\n-----Elapsed time: %.2f seconds.-----\n'%(time.time()-start))
    print('\n-----Elapsed time: %.2f seconds.-----\n'%(time.time()-start), file=f)
    
# Send email when finish processing.
# fromaddr = 'coco0981568491@gmail.com'
# toaddr = ['coco0981568491@gmail.com']

# fromaddr = 'r10b46011@ntu.edu.tw'
# toaddr = ['r10b46011@ntu.edu.tw']

# fromaddr = 'as0200589@gate.sinica.edu.tw'
# toaddr = 'as0200589@gate.sinica.edu.tw'

# msg = MIMEMultipart()

# msg['From'] = fromaddr
# # msg['To'] = ", ".join(toaddr)
# msg['To'] = toaddr
# msg['Subject'] = "[get_rbd_angles processing results]"

# body = "Please refer to the attached Excel file for the processing results.\nThanks!"

# msg.attach(MIMEText(body, 'plain'))

# filename = f'{datetime}_RBD_Angles.xlsx'
# attachment = open(f"/Users/asibc512/RBD_Angle/20230326_RBDangle_analysis/{filename}","rb")

# part = MIMEBase('application', 'octet-stream')
# part.set_payload((attachment).read())
# encoders.encode_base64(part)
# part.add_header('Content-Disposition', "attachment; filename= %s" % filename)
# msg.attach(part)

# # server = smtplib.SMTP('smtp.gmail.com', 587)
# # server = smtplib.SMTP_SSL("smtps.ntu.edu.tw", 465)
# server = smtplib.SMTP('smtp.iis.sinica.edu.tw', 587)

# server.starttls()
# server.login(fromaddr,'Cocoh94m4vu3') #Type Password
# text = msg.as_string()
# server.sendmail(fromaddr, toaddr, text)
# server.quit()

# Upload results to google drive when finish processing.
gauth = GoogleAuth()           
# gauth.LocalWebserverAuth()       
drive = GoogleDrive(gauth)
upload_file = f'/Users/asibc512/RBD_Angle/20230326_RBDangle_analysis/{datetime}_RBD_Angles.xlsx'
gfile = drive.CreateFile({'parents': [{'id': '1mLiIxNOMLaQKmpz05Kua0X0C5Wc7FhNQ'}]})
# Read file and set it as the content of this instance.
gfile.SetContentFile(upload_file)
gfile.Upload() # Upload the file.