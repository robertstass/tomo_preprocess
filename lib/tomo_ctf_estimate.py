#!/usr/bin/env python

import glob
import os
import stat
import sys
import argparse
from distutils import spawn
from tomo_motioncor2 import *
try:
    from pyrelion import MetaData, Item
except:
    print('You must have the "pyrelion/localrec" python library in your $PYTHONPATH')
    sys.exit()




class ArgumentParser():
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                              description="Run's relion's ctffind4 function on tilt series data. "
                                                          "Output's an ordered star file that can be later converted to a "
                                                            "format readable by imod/etomo using tomo_write_imod_defocus.py.")
        general_args = self.parser.add_argument_group('General arguments')
        tilt_args = self.parser.add_argument_group('Tilt info arguments')
        ctf_args = self.parser.add_argument_group('Ctf estimation arguments')
        advanced_args = self.parser.add_argument_group('Advanced arguments (optional)')

        add_g = general_args.add_argument  # shortcut
        add_t = tilt_args.add_argument
        add_c = ctf_args.add_argument
        add_a = advanced_args.add_argument

        add_g('--input_folders', default=None,
              help='A wild card expression (IN QUOTES!) to some folders for batch operations of this script on multiple tilt series')
        add_g('-i', '--input_files', default='?????_??.??.??.mrc',
            help='A wild card expression (IN QUOTES!) to the movie stacks. Note that these should sort alphabetically in the order that they were taken.')
        add_g('--tomo_name', default=default_tomo_name,
            help='An arbitrary name for the output stack, (if using --input_folders the name of the folder is given instead)')


        add_t('--tilt_scheme',
             help='A tilt scheme used from the following list; %s. (more detailed info in the README.txt)' % (', '.join(accepted_tilt_schemes)))
        add_t('--min_angle', type=float,
              help='The minimum (most negative) angle used in the tilt series')
        add_t('--angle_step', type=float, help='The angular step between tilt images')
        add_t('--starting_angle', type=float, default=default_starting_tilt_angle, help='(optional) The starting tilt angle. Only used for bidirectional and dose symmetric tilt schemes.')


        add_c('--rln_version', type=int, default=default_rln_version, help='Relion version 1 or 2 (1 < 1.6 and 2 = 2+)')
        add_c('--ctf_software', type=str, default=default_ctf_software, help='Which ctf estimation software to use. You must input one of the following options %s. NB gctf only works with relion 2+.' % ', '.join(ctf_software_options))
        add_c('--ctf_star', default=default_ctf_star, help='Output star file for the ctf estimation parameters. (NB The extension is ignored to create a directory for relion 2+).')

        add_c('--pixel_size', type=float, default=None, help='Pixel size in angstrom (A)')
        add_c('--CS', type=int, default=default_CS, help="Relion's CS parameter: Spherical aberration (mm)")
        add_c('--HT', type=int, default=default_HT, help="Relion's HT parameter: Voltage (kV)")
        add_c('--AmpCnst', type=float, default=default_AmpCnst, help="Relion's AmpCnst parameter: Fraction of amplitude contrast.")
        add_c('--Box', type=int, default=default_Box, help="Relion's Box parameter: Box size for ctf estimation")
        add_c('--ResMin', type=int, default=default_ResMin, help="Relion's ResMin parameter: Minimum resolution for ctf estimation (A)")
        add_c('--ResMax', type=int, default=default_ResMax, help="Relion's  parameter: Maximum resolution for ctf estimation (A)")
        add_c('--dFMin', type=int, default=default_dFMin, help="Relion's dFMin parameter: Minimum defocus value for ctf estimation (A)")
        add_c('--dFMax', type=int, default=default_dFMax, help="Relion's dFMax parameter: Maximum defocus value for ctf estimation (A)")
        add_c('--FStep', type=int, default=default_FStep, help="Relion's FStep parameter: Defocus search step size")
        add_c('--dAst', type=int, default=default_dAst, help="Relion's dAst parameter: Amount of astigmatism (A)")
        add_c('--ctfWin', type=int, default=default_ctfWin, help="Relion's ctfWin parameter: a squared window of this size at the center of the micrograph will be used to estimate the CTF. set -1 to use whole micrograph")
        add_c('--cores', type=int, default=default_cores, help="Number of cpu cores to use (ctffind only)")
        add_c('--gpu', type=int, default=default_gpu, help="Which gpu's to use. (for gctf only)")
        add_c('--ctf_exe', type=str, default=default_ctf_exe, help="Path to the ctffind/gctf executable.")


        add_a('--only_do_unfinished', action='store_true',
            help="Relion's only_do_unfinished parameter: Only run ctf estimation on those micrographs that do not have logfile with final values.")
        add_a('--custom_tilt_order', default=None, help="A comma delimited list of integers denoting the order that tilts were taken. This can be used with the tilt info options.")
        add_a('--use_tilt_order_files', action='store_true',
              help="Supply tilt information as a 'tilt.order' file in the same directory as the input files. Format: One line per tilt with integers denoting the order they were taken. A second column can optionally be included with dose values (after movie recorded). If the file does not exist the values specified in the 'Tilt info arguments' or by '--custom_tilt_order' will be used.")
        add_a('--write_tilt_order_files', action='store_true', help='Write out tilt.order files based on the given tilt arguments. Existing files will not be overwritten.')
        add_a('--mic_star', default=None,help='Input a star file with micrographs in tilt order and split by tomogram (if so, no tilt scheme options are required.')
        add_a('--only_make_sorted_mic_star', action='store_true', help='Only create the micrograph star file so that ctf estimation can be done with the relion gui.')
        add_a('--only_print_command', action='store_true', help="Only print the ctffind/gctf command to the terminal (don't execute it)")

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

        if args.mic_star == None:
            #Find the files
            if args.input_folders == None and glob.glob(args.input_files) == []:
                self.error('No files found. %s' % (args.input_files))
                sys.exit(2)
            if args.input_folders != None and glob.glob(args.input_folders) == []:
                self.error('No folders found. %s' % (args.input_folders))
                sys.exit(2)

            tilt_info_validation(self, args)

        ctf_estimation_argument_validation(self, args)

        if sys.version_info < (2, 7):
            self.error("Python version 2.7 or later is required.")



def ctf_estimation_argument_validation(self, args):
    if not args.only_make_sorted_mic_star:
        if args.rln_version not in rln_versions:
            self.error('Relion version must be 1 or 2 (1 < 1.6 and 2 = 2+). Not %d.' % (args.rln_version))
        if args.ctf_software not in ctf_software_options:
            self.error('Ctf estimation software %s not supported. Must be from this list; %s.' % (
                args.ctf_software, ', '.join(ctf_software_options)))
        if args.rln_version == 1 and args.ctf_software in rln2_only_ctf_software:
            self.error('Ctf estimation software %s is not supported by relion version 1' % args.ctf_software)

        if args.ctf_software == 'gctf' and args.cores != 1:
            print('Note that the --cores option is ignored for gctf.')
        elif args.ctf_software != 'gctf' and args.gpu != 0:
            print('The --gpu option is ignored for ctffind.')


        if args.pixel_size == None:
            self.error('Must input a pixel size (A)')

        if not args.only_print_command:
            # Find the software dependencies
            if not (spawn.find_executable("relion_run_ctffind")):
                self.error("Relion (relion_run_ctffind) not found.",
                           "Make sure relion programs are in $PATH.")
            if not os.path.exists(args.ctf_exe):
                self.error('Cannot find the ctf executable at path %s' % args.ctf_exe)


#Nitpicky details
default_micrograph_star_file = 'micrographs.star'
default_ctf_star = 'ctf.star'
ctf_software_options = ('ctffind4', 'ctffind3', 'gctf')
rln2_only_ctf_software = (ctf_software_options[2])
rln_versions = (1,2)
default_ctf_software = ctf_software_options[0]
default_rln_version = rln_versions[1]
ctf_estimation_command_log_file = 'tomo_ctf_estimation_command_log.txt'

#Ctffind defaults
default_CS = 2
default_HT = 300
default_AmpCnst = 0.1
default_XMAG = 10000
default_Box = 256
default_ResMin = 50
default_ResMax = 7
default_dFMin = 20000
default_dFMax = 50000
default_FStep = 500
default_dAst = 100
default_ctfWin = -1
default_ctf_exe = "/raid/44/juha/apps/ctffind4/ctffind"
default_cores = 1
default_gpu = 0



XMAG = default_XMAG #default to 10000 so you can just use the pixel size as DStep. This is what relion2.0 does by default.




def write_sorted_mic_star(folder, input_files,  use_tilt_order_files, custom_tilt_order, tilt_scheme, min_angle, angle_step, starting_tilt_angle, tomo_name, write_tilt_order_files):
    file_list = glob.glob(folder + '/' + input_files)
    if sort_by_name:
        file_list = sorted(file_list)
    else:
        file_list.sort(key=os.path.getmtime)
    total_tilts = len(file_list)

    order_list, doses_list, tilt_order_path = parse_tilt_order_and_dose(folder, False, False, use_tilt_order_files,
                                                                        custom_tilt_order, total_tilts, tilt_scheme,
                                                                        min_angle, angle_step,
                                                                        starting_tilt_angle, tomo_name,
                                                                        None, None, write_tilt_order_files)
    md = MetaData()
    for order in order_list:
        file = file_list[order-1]
        p = Item()
        p.rlnMicrographName = file
        md.addItem(p)
    return md


def run_ctffind(micrograph_star_file, ctf_star, pixel_size, ctfWin, CS, HT, AmpCnst, Box, ResMin, ResMax, dFMin, dFMax, FStep, dAst, ctf_exe, cores, only_print_command, rln_version, ctf_software, gpu, only_do_unfinished):
    mpirun_ctffind = 'mpirun -n %d `which relion_run_ctffind_mpi`' % cores
    nompi_ctffind = '`which relion_run_ctffind`'
    ctffind_string = '--i %s --XMAG %d --DStep %f --ctfWin %d --CS %f --HT %d --AmpCnst %f --Box %d --ResMin %d --ResMax %d --dFMin %d --dFMax %d --FStep %d --dAst %d' % (
                        micrograph_star_file, XMAG, pixel_size, ctfWin, CS, HT, AmpCnst, Box, ResMin, ResMax, dFMin, dFMax,FStep, dAst)
    ctffind_string = ctffind_string + ' --only_do_unfinished' if only_do_unfinished else ctffind_string
    if rln_version == 1:
        if ctf_software == 'ctffind3':
            ctffind_string = '%s %s --o %s --ctffind_exe "%s"' % (mpirun_ctffind, ctffind_string, ctf_star, ctf_exe)
        if ctf_software == 'ctffind4':
            ctf_exe = ctf_exe+' --omp-num-threads 1 --old-school-input'
            ctffind_string = '%s %s --o %s --ctffind_exe "%s"' % (mpirun_ctffind, ctffind_string, ctf_star, ctf_exe)
    if rln_version == 2:
        ctf_dir = os.path.splitext(ctf_star)[0]
        if ctf_software == 'ctffind3':
            ctffind_string = '%s %s --o %s --ctffind_exe "%s"' % (mpirun_ctffind, ctffind_string, ctf_dir, ctf_exe)
        if ctf_software == 'ctffind4':
            ctffind_string = '%s %s --o %s --ctffind_exe "%s" --is_ctffind4' % (mpirun_ctffind, ctffind_string, ctf_dir, ctf_exe)
        if ctf_software == 'gctf':
            ctffind_string = '%s %s --o %s --use_gctf --gctf_exe "%s" --angpix %f --gpu %d' % (nompi_ctffind, ctffind_string, ctf_dir, ctf_exe, pixel_size, gpu)
    #ctffind_string = c)tffind_string+' &'
    f = open(ctf_estimation_command_log_file, 'w')
    f.write(ctffind_string+'\n')
    f.close()
    if only_print_command == True:
        print('Not executing: %s' % ctffind_string)
    else:
        print('Executing: %s' % ctffind_string)
        os.system(ctffind_string)




def main(input_files, input_folders, mic_star, tomo_name, tilt_scheme, min_angle, angle_step, pixel_size, only_make_sorted_mic_star, only_do_unfinished,
         use_tilt_order_files, custom_tilt_order, write_tilt_order_files, starting_tilt_angle, only_print_command, rln_version, ctf_software, CS, HT, AmpCnst, Box, ResMin, ResMax, dFMin, dFMax, FStep, dAst, ctfWin, cores, gpu, ctf_exe, ctf_star):

    if mic_star == None:
        if input_folders != None:
            folder_list = sorted(glob.glob(input_folders))
            folder_list = [dir for dir in folder_list if os.path.isdir(dir)]
        else:
            folder = '/'.join(input_files.split('/')[0:-1])
            folder_list = ['.'] if folder == '' else [folder]
            input_files = input_files.split('/')[-1]

        mic_md = MetaData()

        for folder in folder_list:
            if input_folders != None:
                tomo_name = folder.split('/')[-1]

            #Main script for individual tilt series
            tomo_mic_md = write_sorted_mic_star(folder, input_files, use_tilt_order_files, custom_tilt_order, tilt_scheme, min_angle, angle_step, starting_tilt_angle, tomo_name, write_tilt_order_files)
            mic_md.addData(tomo_mic_md)

        mic_md.addLabels('rlnMicrographName')
        mic_md.write(default_micrograph_star_file)

        mic_star = default_micrograph_star_file

    if only_make_sorted_mic_star == True:
        print('Only making sorted micrograph star... quitting.')
        return

    run_ctffind(mic_star, ctf_star, pixel_size, ctfWin, CS, HT, AmpCnst, Box, ResMin, ResMax, dFMin, dFMax,
                FStep, dAst, ctf_exe, cores, only_print_command, rln_version, ctf_software, gpu, only_do_unfinished)

    print("")
    print("You can inspect the ctf estimation quality by running 'relion_display' on the resulting star file with the '--display rlnCtfImage option'")
    print("Once you have aligned the tilt series in imod, use the 'tomo_write_imod_defocus' script to write the estimated ctf parameters into imod compatible .defocus files.")


if __name__ == "__main__":
    argparser = ArgumentParser()
    args = argparser.parser.parse_args()
    argparser.validate(args)

    # main(motioncor2,frames,input_files,input_folders,tomo_name,tilt_scheme,dose_per_movie,min_angle,angle_step,pixel_size,binning,throw,trunc,gpu, only_make_batch_file, only_do_unfinished, skip_building_full_stack)
    main(args.input_files, args.input_folders, args.mic_star, args.tomo_name, args.tilt_scheme, args.min_angle, args.angle_step, args.pixel_size, args.only_make_sorted_mic_star, args.only_do_unfinished, args.use_tilt_order_files, args.custom_tilt_order, args.write_tilt_order_files, args.starting_angle, args.only_print_command, args.rln_version, args.ctf_software,
         args.CS, args.HT, args.AmpCnst, args.Box, args.ResMin, args.ResMax, args.dFMin, args.dFMax, args.FStep, args.dAst, args.ctfWin, args.cores, args.gpu, args.ctf_exe, args.ctf_star)


