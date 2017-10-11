#!/usr/bin/env python


import glob
import os
import stat
import sys
import argparse
from distutils import spawn
try:
    from pyrelion import MetaData, Item
except:
    print('You must have the "pyrelion" python library in your $PYTHONPATH')
    sys.exit()
import math


class ArgumentParser():
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                              description="Converts a relion star file with ctf estimation parameters to a format readable by imod/etomo."
                                                            "The micrographs in the star file must be in tilt order (-ve to +ve) and separated by tiltseries")
        general_args = self.parser.add_argument_group('General arguments')


        add_g = general_args.add_argument  # shortcut

        add_g('-i', '--imod_edf_files', default=None,
              help='A wild card expression (IN QUOTES!) to multiple imod project files (with the ".edf" extension)')
        add_g('--ctf_star',  help='An ordered relion star file with ctf parameters estimated. Micrographs must be ordered (-ve to +ve) and separated by tilt series.')
        add_g('-v', '--version', default='3', help='Which version of the imod defocus file would you like to write? (2 or 3?) Version 3 files include astigmatism values and will only work with later versions of imod (imod_4.9+). You can also use option "star" to write separate relion star files to each folder or "novactf" to write ctffind4-like defocus files to each folder. Astigmatism angles will not be rotated for these 2 options.')


        if len(sys.argv) == 1:  # if no args print usage.
            self.usage()
            sys.exit()

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print "Error: " + '\n'.join(msgs)
        print " "
        sys.exit(2)

    def validate(self, args):
        # Find the files
        if glob.glob(args.imod_edf_files) == []:
            self.error('No imod edf files found. "%s"' % (args.imod_edf_files))
            sys.exit(2)
        if not os.path.isfile(args.ctf_star):
            self.error('No star file found. "%s"' % (args.ctf_star))
            sys.exit(2)
        if args.version not in accepted_versions:
            self.error('version %s not supported. Must be from this list; %s.' % (args.version, ', '.join(accepted_versions)))
            sys.exit(2)





#Nitpicky details
accepted_versions = ('2', '3', 'novactf', 'star')
star_append = '_defocus.star'
novactf_append = '_novactf_defocus.txt'
ctffind4_header_line = '# Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; #4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; #7 - spacing (in Angstroms) up to which CTF rings were fit successfully' #may not be needed



def split_by_folder(md):
    current_dir = os.path.dirname(md._data[0].rlnMicrographName) if os.path.dirname(md._data[0].rlnMicrographName) != "" else '.'
    temp_md = MetaData()
    md_dict = {}
    for p in md:
        folder = os.path.dirname(p.rlnMicrographName) if os.path.dirname(p.rlnMicrographName) != "" else '.'
        if folder == current_dir:
            temp_md.addItem(p)
        else:
            temp_md.addLabels(md.getLabels())
            md_dict[current_dir] = temp_md
            current_dir = folder
            temp_md = MetaData()
            temp_md.addItem(p)
    temp_md.addLabels(md.getLabels())
    md_dict[current_dir] = temp_md
    return md_dict



def read_tilt_angles(tlt_path):
    f = open(tlt_path,'r')
    lines = f.readlines()
    tilt_list = []
    for i in range(0,len(lines)):
        tilt_list.append(float(lines[i][0:-1]))
    f.close()
    return tilt_list

def read_xf(xf_path):
    f = open(xf_path,'r')
    lines = f.readlines()
    transform_list = []
    for i in range(0,len(lines)):
        transform = lines[i][2:-1]
        t_list = transform.split()
        transform_list.append([float(i) for i in t_list])
    f.close()
    return transform_list


def xfToRot(a11, a12, a21, a22, dx, dy):
    #This function uses translated code from IMOD's amat_to_rotmagstr.c function to find the .ali to .st angle from lines in the .xf file
    #detect axis rotations
    dtheta = math.degrees(math.atan2(a22, a12) - math.atan2(a21, a11))
    if dtheta > 180:
        dtheta = dtheta - 360
    if dtheta <= -180:
        dtheta = dtheta + 360
    if dtheta < 0:
        a12 = -a12
        a22 = -a22
    #find Rot (theta)
    if a21 != a12 or a22 != -a11:
        theta = math.degrees(math.atan2((a21 - a12), (a22 + a11)))
    else:
        theta = 0
    return theta


def get_imod_file_root(tomo_name, imod_folder):
    if tomo_name == None:
        tomo_name = os.path.basename(imod_folder)
    slash_noslash = '' if imod_folder == "" else '/'
    imod_file_root = '%s%s%s' % (imod_folder, slash_noslash, tomo_name)
    return imod_file_root


def write_imod_defocus_file(md, imod_folder, tomo_name, version):
    imod_file_root = get_imod_file_root(tomo_name, imod_folder)
    defocus_file_name = '%s.defocus' % (imod_file_root)
    tlt_file_name = '%s.tlt' % (imod_file_root)
    xf_file_name = '%s.xf' % (imod_file_root)
    if not os.path.isfile(tlt_file_name):
        print('Cannot find "%s". Make sure the tilt series has been aligned before executing this script!' % tlt_file_name)
        return
    tilt_list = read_tilt_angles(tlt_file_name)

    if version == 3:
        if not os.path.isfile(xf_file_name):
            print('Cannot find "%s". Make sure the tilt series has been aligned before executing this script!' % xf_file_name)
            return
        transform_list = read_xf(xf_file_name)
    os.rename(defocus_file_name, defocus_file_name + '~') if os.path.isfile(defocus_file_name) else None
    f = open(defocus_file_name, 'w')
    if len(md) != len(tilt_list):
        print('Number of tilts in .tlt does not match star file!')

    if version == 3:
        defocus_file_flag = 1
        phase_shifts_included = False
        if 'rlnPhaseShift' in md.getLabels():
            phase_shifts_included = True
            defocus_file_flag = defocus_file_flag+4
        f.write('%d 0  0. 0. 0 3\n' % defocus_file_flag)
        for i, (p, tlt, xform) in enumerate(zip(md, tilt_list, transform_list)):
            defocusU = p.rlnDefocusU / 10
            defocusV = p.rlnDefocusV / 10
            st2ali_angle = xfToRot(xform[0], xform[1], xform[2], xform[3], xform[4], xform[5])
            defocusAng = p.rlnDefocusAngle+st2ali_angle
            if defocusU > defocusV: #lowest defocus value first
                tempV = defocusU
                defocusU = defocusV
                defocusV = tempV
                defocusAng = defocusAng - 90
            line = '%d\t%d\t%f\t%f\t%f\t%f\t%f' % (i + 1, i + 1, tlt, tlt, defocusU, defocusV, defocusAng)
            if phase_shifts_included:
                zero_to_360_phase_shift = float(p.rlnPhaseShift) % 360 #imod does not deal with negative phase shifts well i think so convert them.
                line = '%s\t%f' % (line, zero_to_360_phase_shift)
            f.write(line + '\n')
    elif version ==2:
        for i, (p, tlt) in enumerate(zip(md, tilt_list)):
            avg_defocus = (0.5 * (p.rlnDefocusU + p.rlnDefocusV)) / 10  # in nm
            line = '%d\t%d\t%f\t%f\t%f' % (i + 1, i + 1, tlt, tlt, avg_defocus)
            line = line + '\t2' if i == 0 else line
            f.write(line + '\n')
    f.close()
    print('An imod defocus file (version %d) has been created at %s' % (version, defocus_file_name))


def write_separate_stars(md, imod_folder, tomo_name, append):
    imod_file_root = get_imod_file_root(tomo_name, imod_folder)
    star_path = '%s%s' % (imod_file_root, star_append)
    md.write(star_path)
    print('Star file written to %s' % star_path)

def write_nova_ctf_defocus(md, imod_folder, tomo_name, append):
    imod_file_root = get_imod_file_root(tomo_name, imod_folder)
    defocus_path = '%s%s' % (imod_file_root, novactf_append)
    f = open(defocus_path, 'w')
    f.write(ctffind4_header_line+'\n')
    for i,p in enumerate(md):
        mic_num = i+1
        defocusU = p.rlnDefocusU
        defocusV = p.rlnDefocusV
        defocusAngle = p.rlnDefocusAngle
        if 'rlnPhaseShift' in md.getLabels():
            phase_shift = math.radians(float(p.rlnPhaseShift))
        else:
            phase_shift = 0.0
        cc = float(p.rlnCtfFigureOfMerit) if 'rlnCtfFigureOfMerit' in md.getLabels() else 0
        ctf_maxres = float(p.rlnCtfMaxResolution) if 'rlnCtfMaxResolution' in md.getLabels() else 0
        line = '%d\t%f\t%f\t%f\t%f\t%f\t%f' % (mic_num, defocusU, defocusV, defocusAngle, phase_shift, cc, ctf_maxres)
        f.write(line+'\n')
    f.close()
    print('NovaCTF defocus file (in the style of ctffind4) written to %s' % defocus_path)




def main(imod_edf_files, ctf_star, version):
    edf_list = sorted(glob.glob(imod_edf_files))
    edf_folders = [os.path.dirname(i) for i in edf_list]
    edf_folders = ['.'] if edf_folders == [''] else edf_folders
    edf_names = [os.path.splitext(os.path.basename(i))[0] for i in edf_list]
    ctf_md = MetaData(ctf_star)
    md_dict = split_by_folder(ctf_md)
    for edf_folder, edf_name in zip(edf_folders,edf_names):
        try:
            md = md_dict[edf_folder]
        except:
            print('imod folder "%s" not found in the star file.' % edf_folder)
            continue
        if version == '2' or version == '3':
            write_imod_defocus_file(md, edf_folder, edf_name, int(version))
        elif version == 'novactf':
            write_nova_ctf_defocus(md, edf_folder, edf_name, novactf_append)
        elif version == 'star':
            write_separate_stars(md, edf_folder, edf_name, star_append)

if __name__ == "__main__":
    argparser = ArgumentParser()
    args = argparser.parser.parse_args()
    argparser.validate(args)

    main(args.imod_edf_files, args.ctf_star, args.version)