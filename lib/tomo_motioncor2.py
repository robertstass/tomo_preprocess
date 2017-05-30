#!/usr/bin/env python

import glob
import os
import stat
import sys
import argparse
from distutils import spawn


class ArgumentParser():
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                              description="Applies motion correction to a tiltseries using motioncor2 "
                                                          "outputing a .st image stack in the correct order. Requires motioncor2 and bsoft to be sourced. "
                                                          "There are also options to apply dose weighting. "
                                                          "A typical example of how to use this script; "
                                                        "tomo_motioncor2 --input_folders \"??\" -i \"?????_??.??.??.mrc\" --tilt_scheme  continuous_positive --min_angle -30 --angle_step 3 -apix 1.35 --dose_per_movie 4.0 --do_custom_doseweighting ."
                                                        "(another useful argument initially is --only_make_batch_file that will create a batch file that can be checked before actually running) "
                                                        "(More detailed info in the README.txt)")
        general_args = self.parser.add_argument_group('General arguments')
        tilt_args = self.parser.add_argument_group('Tilt info arguments')
        motioncor_args = self.parser.add_argument_group('Motioncor2 arguments (optional)')
        dw_args = self.parser.add_argument_group('Dose weighting arguments (optional)')
        advanced_args = self.parser.add_argument_group('Advanced arguments (optional)')

        add_g = general_args.add_argument  # shortcut
        add_t = tilt_args.add_argument
        add_m = motioncor_args.add_argument
        add_d = dw_args.add_argument
        add_a = advanced_args.add_argument

        add_g('--input_folders', default=None,
              help='A wild card expression (IN QUOTES!) to some folders for batch operations of this script on multiple tilt series')
        add_g('-i', '--input_files', default='?????_??.??.??.mrc',
            help='A wild card expression (IN QUOTES!) to the movie stacks. Note that these should sort alphabetically in the order that they were taken.')
        add_g('--tomo_name', default='tomo',
            help='An arbitrary name for the output stack, (if using --input_folders the name of the folder is given instead)')

        add_t('--tilt_scheme',
             help='A tilt scheme used from the following list; %s. (more detailed info in the README.txt)' % (', '.join(accepted_tilt_schemes)))
        add_t('--min_angle', type=float,
              help='The minimum (most negative) angle used in the tilt series')
        add_t('--angle_step', type=float, help='The angular step between tilt images')
        add_t('--starting_angle', type=float, default=default_starting_tilt_angle, help='(optional) The starting tilt angle. Only used for bidirectional and dose symmetric tilt schemes.')

        add_m('--motioncor2', default=default_motioncor2_exe, help='Location of the motioncor2 executable')
        add_m('--throw', default=0, type=int, help='Discard this number of frames from the start of each movie')
        add_m('--trunc', default=0, type=int, help='Discard this number of frames from the end of each movie')
        add_m('--binning', default=1, type=int, help='Bin the images')
        add_m('--gpu', default=0, type=int, help='The id of the gpu to use')
        add_m('--patch', default=default_patch, type=int, help="motioncor2 Patch option. x will be given as 'x x'")
        add_m('--iterations', default=default_iterations, type=int, help='motioncor2 Iter option.')


        add_d('-dose', '--dose_per_movie', type=float, default=None,
              help='Dose applied per movie. (this is divided by the number of frames to input into motioncor2')
        add_d('--pre_dose', default=0, type=float, help='Pre exposure before tomogram taken.')
        add_d('--frames', type=int, help='The number of frames per movie')
        add_d('-apix', '--pixel_size', type=float, help='The pixel size of the images in angstrom')
        add_d('--do_motioncor2_doseweighting', action='store_true',
              help="Use motioncor2's dose weighting (with -InitDose). Note this doesn't currently work properly!")
        add_d('--do_custom_doseweighting', action='store_true',
              help="Use this option to use the tomo_dose_filter.py script for dose weighting each tilt series")

        add_a('--only_do_unfinished', action='store_true',
            help='Only run motion correction on those that do not have a motion corrected image already. (by deleting bad examples, you can easily rerun the script with different settings for the deleted images only)')
        add_a('--custom_tilt_order', default=None, help="A comma delimited list of integers denoting the order that tilts were taken. This can be used with the tilt info options.")
        add_a('--use_tilt_order_files', action='store_true',
              help="Supply tilt information as a 'tilt.order' file in the same directory as the input files. Format: One line per tilt with integers denoting the order they were taken. A second column can optionally be included with dose values (after movie recorded). If the file does not exist the values specified in the 'Tilt info arguments' or by '--custom_tilt_order' will be used.")
        add_a('--write_tilt_order_files', action='store_true', help='Write out tilt.order files based on the given tilt arguments. Existing files will not be overwritten.')
        add_a('--only_make_batch_file', action='store_true', help='Only create the list of commands to execute but do not execute them. These can be modified and run separately.')


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
        #Find the files
        if args.input_folders == None and glob.glob(args.input_files) == []:
            self.error('Error: No files found.' % (args.input_files))
            sys.exit(2)
        if args.input_folders != None and glob.glob(args.input_folders) == []:
            self.error('Error: No folders found.' % (args.input_folders))
            sys.exit(2)

        #Check tilt order information
        tilt_info_required = True
        dose_info_required = True
        if args.use_tilt_order_files:
            folder_list = sorted(glob.glob(args.input_folders)) if args.input_folders != None else (['.'] if '/'.join(args.input_files.split('/')[0:-1]) == '' else ['/'.join(args.input_files.split('/')[0:-1])])
            order_files = [folder + '/' + tilt_order_filename for folder in folder_list]
            if any([not os.path.isfile(order_file) for order_file in order_files]):
                tilt_info_required = True
                dose_info_required = True
            else:
                tilt_info_required = False
                tilt_order_info_list = [read_tilt_order(order_file) for order_file in order_files]
                if not any([doses_list == None for order_list, doses_list in tilt_order_info_list]):
                    dose_info_required = False


        tilt_info_args = ('tilt_scheme', 'min_angle', 'angle_step', 'starting_angle')

        if tilt_info_required:
            if args.custom_tilt_order != None:
                if any([getattr(args, arg) != None for arg in tilt_info_args]):
                    self.error('Tilt info arguments cannot be used with the --custom_tilt_order option. (ie %s)' % (', '.join(tilt_info_args)))
            else:
                if any([getattr(args, arg) == None for arg in tilt_info_args]):
                    self.error('Tilt info arguments are required. (ie. %s)' % (
                        ', '.join(tilt_info_args)))
                if args.tilt_scheme not in accepted_tilt_schemes:
                    self.error('Tilt scheme (%s) not supported.\nAccepted tilt schemes are %s' % (args.tilt_scheme, ', '.join(accepted_tilt_schemes)))
                if args.tilt_scheme not in starting_angle_tilt_schemes and args.starting_angle != default_starting_tilt_angle:
                    self.error('Starting angle not required with %s tilt scheme' % (args.tilt_scheme))

        do_any_doseweighting = args.do_motioncor2_doseweighting or args.do_custom_doseweighting


        if do_any_doseweighting:
            if args.pixel_size == None:
                self.error('Missing argument: --pixel_size')
            if dose_info_required:
                if args.dose_per_movie == None:
                    self.error('Missing argument: -dose or --dose_per_movie')
            if args.do_motioncor2_doseweighting and (args.frames == None or args.dose_per_movie == None):
                self.error('Missing argument: --frames and --dose_per_movie/-dose is required when doing motioncor2 dose weighting')


        #Find the software dependencies
        if not (spawn.find_executable("bimg")):
            self.error("Bsoft (bimg) not found.",
                       "Make sure Bsoft programs are in $PATH.")
        if not (spawn.find_executable(args.motioncor2)):
            self.error("%s not found." % (args.motioncor2),
                       "Make sure motioncor2 is in $PATH.")
        if args.do_custom_doseweighting:
            if not (spawn.find_executable(custom_doseweight_script)):
                self.error("%s not found." % (custom_doseweight_script),
                           "Make sure dose_filter is in $PATH.")














## Nitpicky details
default_motioncor2_exe = 'MotionCor2'
motioncor_file_suffix = 'motioncor_tilt_'
dose_weight_suffix = '_DW'  # set by motioncor
motioncor_file_name = 'tomo_motioncor2'
all_frames_suffix = '_all_frames'
default_starting_tilt_angle = 0  # only needed for the bidirectional tilt schemes
batch_file_name = 'batch_tomo_motioncor2.sh'
#motioncor options
default_patch = 4
default_iterations = 3
default_kv = 300
kv = default_kv

custom_doseweight_script = 'tomo_dose_filter'
tilt_order_filename = 'tilt.order'
overwrite_existing_tilt_order_files = False


#####



accepted_tilt_schemes = ['continuous_positive', 'continuous_negative', 'bidirectional_positive',
                         'bidirectional_negative', 'dose_symmetric_positive', 'dose_symmetric_negative']
starting_angle_tilt_schemes = ['bidirectional_positive', 'bidirectional_negative', 'dose_symmetric_positive', 'dose_symmetric_negative']


def tilt_order(tilt_scheme, tilt_num, total_tilts, zero_tilt_index):
    centre = zero_tilt_index + 1
    if tilt_scheme == 'continuous_positive':
        tilt_order = tilt_num
    elif tilt_scheme == 'continuous_negative':
        tilt_order = total_tilts - tilt_num + 1
    elif tilt_scheme == 'bidirectional_positive':
        if tilt_num >= centre:
            tilt_order = tilt_num - centre + 1
        else:
            tilt_order = total_tilts - tilt_num + 1
    elif tilt_scheme == 'bidirectional_negative':
        if tilt_num <= centre:
            tilt_order = centre - tilt_num + 1
        else:
            tilt_order = tilt_num

    elif tilt_scheme == 'dose_symmetric_positive':
        distance_from_minus_end = centre
        distance_from_plus_end = total_tilts - centre + 1
        distance_from_edge = min(distance_from_minus_end, distance_from_plus_end)
        subseries_size = (2*distance_from_edge)-1
        subseries_offset = centre-distance_from_edge
        sub_tilt_num = tilt_num - subseries_offset
        sub_centre = centre - subseries_offset
        if sub_tilt_num <= sub_centre:
            tilt_order = (subseries_size + 1 - (sub_tilt_num*2-1))
        else:
            tilt_order = (sub_tilt_num - sub_centre)*2
        if tilt_num <= subseries_offset:
            tilt_order = total_tilts - tilt_num + 1
        if tilt_num >= centre+distance_from_edge:
            tilt_order = tilt_num

    elif tilt_scheme == 'dose_symmetric_negative':
        distance_from_minus_end = centre
        distance_from_plus_end = total_tilts - centre + 1
        distance_from_edge = min(distance_from_minus_end, distance_from_plus_end)
        subseries_size = (2 * distance_from_edge) - 1
        subseries_offset = centre - distance_from_edge
        sub_tilt_num = tilt_num - subseries_offset
        sub_centre = centre - subseries_offset
        if sub_tilt_num < sub_centre:
            tilt_order = (subseries_size + 1 - (sub_tilt_num*2))
        else:
            tilt_order = (sub_tilt_num - sub_centre )*2 + 1
        if tilt_num <= subseries_offset:
            tilt_order = total_tilts - tilt_num + 1
        if tilt_num >= centre+distance_from_edge:
            tilt_order = tilt_num
    else:
        print('Tilt scheme not supported')
        raise ValueError
    return tilt_order


def is_unfinished(filename, dose_weight_append, do_motioncor2_doseweighting):
    fileroot = '.'.join(filename.split('.')[0:-1])
    if not os.path.isfile(filename) or (not os.path.isfile(fileroot + dose_weight_append + '.mrc') and do_motioncor2_doseweighting):
        return True
    else:
        return False


def make_executable(filename):
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)


def tilt_order_from_tilt_scheme(tilt_scheme, min_angle, angle_step, total_tilts, starting_tilt_angle):
    max_angle = min_angle + (angle_step * (total_tilts - 1))
    tilt_list = [min_angle + (angle_step * i) for i in range(0, total_tilts)]
    zero_tilt_index = tilt_list.index(min(tilt_list, key=lambda x: abs(x - starting_tilt_angle)))
    order_list = [tilt_order(tilt_scheme, i + 1, total_tilts, zero_tilt_index) for i in range(0, total_tilts)]
    starting_angle_string = ' starting from %d degrees' % starting_tilt_angle if tilt_scheme in starting_angle_tilt_schemes else ''
    print('Tiltseries from %d to %d degrees in steps of %d degrees%s (%d tilts in total) using a %s tilt scheme.' % (
        min_angle, max_angle, angle_step, starting_angle_string, total_tilts, tilt_scheme))
    return order_list


def display_number_list(number_list):
    return ', '.join((str(e) for e in number_list))

def validate_tilt_order(tilt_order, total_tilts, tomo_name):
    n = len(tilt_order)
    if n > total_tilts:
        tilt_order = tilt_order[0:total_tilts]
        print('Tilt order supplied (for %s) is too long. Truncating to %s' % (tomo_name, display_number_list(tilt_order)))
    elif n < total_tilts:
        print('Tilt order supplied (for %s) is too short. Quitting...' % tomo_name)
        sys.exit()
    sorted_list = sorted(tilt_order)
    if sorted_list != range(1,total_tilts+1): #check if it's a list of integers.
        print('Tilt order supplied (for %s) must be a continuous set of integers from 1. Quitting...' % tomo_name)
        sys.exit()
    return tilt_order


def validate_doses(doses_list, total_tilts, tomo_name):
    n = len(doses_list)
    if n > total_tilts:
        doses_list = doses_list[0:doses_list]
        print('Dose list supplied (for %s) is too long. Truncating to %s' % (tomo_name, display_number_list(doses_list)))
    elif n < total_tilts:
        print('Dose list supplied (for %s) is too short. Quitting...' % tomo_name)
        sys.exit()
    doses_list = [float(dose) for dose in doses_list]
    return doses_list


def read_tilt_order(tilt_order_path):
        f = open(tilt_order_path, 'r')
        lines = f.readlines()
        tilt_order = []
        tilt_doses = []
        for line in lines:
            value = line.split()
            if len(value) > 0:
                tilt_order.append(value[0])
            if len(value) > 1:
                tilt_doses.append(value[1])
        try:
            tilt_order = [int(item) for item in tilt_order]
        except ValueError:
            print('Tilt order file must be a set of integers')
            sys.exit()
        if len(tilt_doses) == len(tilt_order):
            try:
                tilt_doses= [float(item) for item in tilt_doses]
            except ValueError:
                print('Tilt order file doses must be a set of floats')
                sys.exit()
        elif len(tilt_doses) > 0:
            print('Tilt order file is formatted incorrectly')
        else:
            tilt_doses = None
        return tilt_order, tilt_doses

def read_custom_tilt_order(custom_tilt_order):
    try:
        order_list = [int(item) for item in custom_tilt_order.split(',')]
    except ValueError:
        print('Tilt order must be a set of integers')
        sys.exit()
    return order_list


def write_tilt_order(order_list, doses_list, folder):
    filename = folder +'/' + tilt_order_filename
    f = open(filename, 'w')
    if doses_list == None:
        outlist = ['%d\n' % (order) for order in order_list]
    else:
        outlist = ['%d\t%f\n' % (order, dose) for order, dose in zip(order_list, doses_list)]
    for tilt in outlist:
        f.write(tilt)
    f.close()


def tomogram_motioncor2(motioncor2, frames, input_files, folder, tomo_name, tilt_scheme, dose_per_movie, min_angle, angle_step,
         pixel_size, binning, throw, trunc, gpu, only_do_unfinished, do_motioncor2_doseweighting,
         do_custom_doseweighting, pre_dose, use_tilt_order_files, custom_tilt_order, write_tilt_order_files, patch, iterations, starting_tilt_angle):
    #Read files
    file_list = sorted(glob.glob(folder + '/' + input_files))
    total_tilts = len(file_list)

    #Parse tilt information
    print('\n')
    print('Tilt Series; %s' % tomo_name)
    do_any_doseweighting = do_motioncor2_doseweighting or do_custom_doseweighting
    doses_list = None
    tilt_order_path = folder + '/' + tilt_order_filename
    if use_tilt_order_files and os.path.isfile(tilt_order_path):
        order_list, doses_list = read_tilt_order(tilt_order_path)
        order_list = validate_tilt_order(order_list, total_tilts, tomo_name)
        print('Using tilt order from file; %s' % display_number_list(order_list))
    elif custom_tilt_order != None:
        order_list = read_custom_tilt_order(custom_tilt_order)
        order_list = validate_tilt_order(order_list, total_tilts, tomo_name)
        print('Using custom tilt order; %s' % display_number_list(order_list))
    else:
        order_list = tilt_order_from_tilt_scheme(tilt_scheme, min_angle, angle_step, total_tilts, starting_tilt_angle)
        order_list = validate_tilt_order(order_list, total_tilts, tomo_name)

    #Parse doses
    if do_any_doseweighting and doses_list == None:
        doses_list = [(item*dose_per_movie)+pre_dose for item in order_list]
        doses_list = validate_doses(doses_list, total_tilts, tomo_name)

    if do_any_doseweighting:
        print('Applying dose weighting; %s' % display_number_list(doses_list))
        if do_motioncor2_doseweighting:
            print("Motioncor2 doseweighting on. Note that motioncor2 doseweighting with -InitDose currently doesn't work! Also doses given will be minus the dose_per_movie to give initial doses")

    if write_tilt_order_files and (not os.path.isfile(tilt_order_path) or overwrite_existing_tilt_order_files):
        write_tilt_order(order_list, doses_list, folder)


    tomo_root = folder + '/' + tomo_name
    motioncor_file = folder + '/' + motioncor_file_name + '.sh'

    f = open(motioncor_file, 'w')
    f.write('#!/bin/sh\n')
    f.write('echo "Starting motion correction ... "\n')


    motioncor_image_paths = [None for i in range(0, total_tilts)]
    dw_motioncor_image_paths = [None for i in range(0, total_tilts)]
    for i, movie in enumerate(file_list):
        order_string = str(order_list.index(i + 1) + 1).zfill(2)
        motioncor_image_path = os.path.splitext(movie)[0] + '_' + motioncor_file_suffix + order_string + '.mrc'
        motioncor_image_paths[order_list.index(i + 1)] = motioncor_image_path
        dw_motioncor_image_paths[order_list.index(i + 1)] = os.path.splitext(motioncor_image_path)[0] + dose_weight_suffix + '.mrc'
        motioncor_line = '%s -InMrc %s -OutMrc %s -FtBin %d -Patch %d %d -Gpu %d -Iter %d -Throw %d -Trunc %d' % (
            motioncor2, movie, motioncor_image_path, binning, patch, patch, gpu, iterations, throw, trunc)
        if do_motioncor2_doseweighting:
            dose_per_frame = dose_per_movie / frames
            initdose = doses_list[order_list.index(i + 1)] - dose_per_movie #it's the dose at the start that is needed.
            if initdose < 0:
                print('Warning: initial dose can not be less than zero. Check tilt order file (for %s).' % (tomo_name))
                initdose = 0
            motioncor_line = '%s -kV %d -FmDose %f -InitDose %f -PixSize %f' % (motioncor_line, kv, dose_per_frame, initdose, pixel_size)
        if only_do_unfinished and not is_unfinished(motioncor_image_path, dose_weight_suffix, do_motioncor2_doseweighting):
            motioncor_line = '# ' + motioncor_line
        f.write(motioncor_line + '\n')

    bcat_line = 'bcat -output %s:mrc %s' % (tomo_root + '.st', ' '.join(motioncor_image_paths))
    f.write('echo "' + bcat_line + '"\n')
    f.write(bcat_line + '\n')
    if do_motioncor2_doseweighting:
        bcat_line = 'bcat -output %s:mrc %s' % (
            tomo_root + dose_weight_suffix + '.st', ' '.join(dw_motioncor_image_paths))
        f.write('echo "' + bcat_line + '"\n')
        f.write(bcat_line + '\n')
    if do_custom_doseweighting:
        custom_doseweight_line = '%s --tilt_series %s --pixel_size %f' % (custom_doseweight_script, tomo_root + '.st', pixel_size * binning)
        if (use_tilt_order_files and os.path.isfile(tilt_order_path)) or custom_tilt_order != None:
            custom_doseweight_line = '%s --custom_dose_series "%s"' % (custom_doseweight_line, display_number_list(doses_list))
        else:
            custom_doseweight_line = '%s --tilt_scheme %s --dose_per_tilt %f --min_angle %f --angle_step %f --starting_angle %f --pre_dose %f' % (
                custom_doseweight_line, tilt_scheme, dose_per_movie, min_angle, angle_step, starting_tilt_angle, pre_dose)
        f.write('echo "' + custom_doseweight_line + '"\n')
        f.write(custom_doseweight_line + '\n')
    f.close()
    make_executable(motioncor_file)
    return motioncor_file



def main(motioncor2, frames, input_files, input_folders, tomo_name, tilt_scheme, dose_per_movie, min_angle, angle_step,
         pixel_size, binning, throw, trunc, gpu, only_make_batch_file, only_do_unfinished, do_motioncor2_doseweighting,
         do_custom_doseweighting, pre_dose, use_tilt_order_files, custom_tilt_order, write_tilt_order_files, patch, iterations,
         starting_tilt_angle):
    if input_folders != None:
        folder_list = sorted(glob.glob(input_folders))
        batch_f = open(batch_file_name, 'w')
    else:
        folder = '/'.join(input_files.split('/')[0:-1])
        folder_list = ['.'] if folder == '' else [folder]
        input_files = input_files.split('/')[-1]

    for folder in folder_list:
        if input_folders != None:
            tomo_name = folder.split('/')[-1]
        #Main script for individual tilt series
        motioncor_file = tomogram_motioncor2(motioncor2, frames, input_files, folder, tomo_name, tilt_scheme, dose_per_movie, min_angle, angle_step,
         pixel_size, binning, throw, trunc, gpu, only_do_unfinished, do_motioncor2_doseweighting,
         do_custom_doseweighting, pre_dose, use_tilt_order_files, custom_tilt_order, write_tilt_order_files, patch, iterations, starting_tilt_angle)

        if input_folders != None:
            batch_f.write(motioncor_file + '\n')

    if input_folders != None:
        batch_f.close()
        make_executable(batch_file_name)

    if only_make_batch_file != True:
        if input_folders == None:
            os.system('./' + motioncor_file)
        else:
            os.system('./' + batch_file_name)






if __name__ == "__main__":
    argparser = ArgumentParser()
    args = argparser.parser.parse_args()
    argparser.validate(args)

    # main(motioncor2,frames,input_files,input_folders,tomo_name,tilt_scheme,dose_per_movie,min_angle,angle_step,pixel_size,binning,throw,trunc,gpu, only_make_batch_file, only_do_unfinished, skip_building_full_stack)
    main(args.motioncor2, args.frames, args.input_files, args.input_folders, args.tomo_name, args.tilt_scheme,
         args.dose_per_movie, args.min_angle, args.angle_step, args.pixel_size, args.binning, args.throw, args.trunc,
         args.gpu, args.only_make_batch_file, args.only_do_unfinished, args.do_motioncor2_doseweighting,
         args.do_custom_doseweighting, args.pre_dose, args.use_tilt_order_files, args.custom_tilt_order,
         args.write_tilt_order_files, args.patch, args.iterations, args.starting_angle)










