#!/usr/bin/env python

import mrcfile
import numpy as np
import sys
import os
import glob
import argparse


use_pfftw = True # This provides a faster FFT than the standard numpy version. It is optional.
if use_pfftw:
    try:
        import pyfftw
    except ImportError:
        print('Cannot find pyfftw module. Install this to get faster FFTs!')
        use_pfftw = False

class ArgumentParser():
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                              description="Applies sharpening to account for the amplification of low resolution information induced by dose weighting a tilt series."
                                              "This can be applied to 2d image stacks, 3d subvolumes or final maps")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        addr('-i', '--input_map', required=True, help='A wild card expression (IN QUOTES!) to the mrc maps that you want to sharpen. (also works on 2d images/stacks).')
        addr('--number_of_tilts', type=int, help='The number of tilts included in the reconstruction. This can be less than the total number collected so long as those removed were the last to be collected.')
        addr('-dose', '--dose_per_tilt', type=float, help='Dose applied per tilt image. (in e-/A^2)')
        add('--all_inputs_equal', action='store_true', help='Set this option if all input maps are equal (in size and pixel size). Arrays are precalculated to speed things up!')

        add('-apix', '--pixel_size', type=float, help='The pixel size of the images in angstrom. If not supplied, apix is read from the header.')
        add('--pre_dose', default=0, type=float, help='Initial dose before tilt series collected.')
        add('--in_place', action='store_true', help='Make the changes to the existing file.')
        add('--file_append', default='_dw_sharpened', type=str, help='String to append to the end of the file name.')
        add('--do_not_do_dose_weighting', action='store_true', help='Set this to just check the files and not actually apply dose weighting.')
        add('--interpret_as_slices', action='store_true', help='Force interpreting a 3d volume as a 2d image stack')
        add('--interpret_as_images', action='store_true', help='Force interpreting a stack of 2d images as a 3d volume')
        add('--verbosity', type=int, default=default_verbosity_level, help='verbosity level')

        if len(sys.argv) == 1:  # if no args print usage.
            self.usage()
            sys.exit()

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):

        if sys.version_info < (2, 7):
            self.error("Python version 2.7 or later is required.")


# Nitpicky details
starting_tilt_angle = 0  # only needed for the bidirectional tilt schemes
plot_filters = [] #[0] #[0]  # [0,10,29] #[0,1,2]
keep_header_apix = True #use the original pixel size in the header of the output file. (apix is still used for the dose weighting). This avoids mismatches in pixel size between the input and output stacks.
default_starting_tilt_angle = 0
default_verbosity_level = 3
####
if plot_filters != []:  # only use if matplotlib available
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D


def verbosity_print(verbosity, desired_verbosity, string):
    if verbosity >= desired_verbosity:
        print(string)


class DoseWeightSharpen:
    def __init__(self, input_file, dose_per_tilt, pre_dose, number_of_tilts, apix, interpret_as_slices, interpret_as_images, file_append, in_place, copy_file_init_from=None, plot_filters=[]):
        self.a = 0.245
        self.b = -1.665
        self.c = 2.81
        self.file_path = input_file
        self.dose_per_tilt = dose_per_tilt
        self.pre_dose = pre_dose
        self.number_of_tilts = number_of_tilts
        self.apix = apix
        self.interpret_as_slices = interpret_as_slices
        self.interpret_as_images = interpret_as_images
        self.in_place = in_place
        self.file_append = file_append
        self.plot_filters = plot_filters # list of indices to plot from the images list
        #self.plot_filters = [x for x in self.plot_filters if x >= 0 and x < self.number_of_files]  # remove nonsense values
        if copy_file_init_from==None:
            self.mrc, self.filter_array, self.apix, self.is_single_image, self.is_image_stack, self.is_single_volume, self.is_volume_stack = self.init_dw_sharpen(self.file_path, self.interpret_as_slices, self.interpret_as_images)
            self.copied_file_init = False
        else:
            self.copy_file_init_from_other(copy_file_init_from)
            self.copied_file_init = True
            self.mrc = None

    # self.dose_weight()

    def init_dw_sharpen(self, input_file, interpret_as_slices, interpret_as_images):
        mrcfile_open_args = [input_file]
        mrcfile_open_kwargs = {'mode': 'r'}
        mrc = mrcfile.open(*mrcfile_open_args, **mrcfile_open_kwargs)
        is_single_image, is_image_stack, is_single_volume, is_volume_stack = self.is_stack_or_volume(mrc, interpret_as_slices, interpret_as_images)
        is_stack = is_image_stack or is_volume_stack
        apix = self.multidim_apix_to_single_value(mrc.voxel_size) if self.apix == None else self.apix
        xsize = int(mrc.header.nx)
        ysize = int(mrc.header.ny)
        zsize = int(mrc.header.nz)
        zyx_shape = (zsize, ysize, xsize)
        if is_stack or is_single_image:
            freq_array_shape = zyx_shape[1:]
        else:
            freq_array_shape = zyx_shape
        verbosity_print(verbosity, 2, 'Calculating frequency array...')
        freq_array = self.create_frequency_array(freq_array_shape, apix)
        verbosity_print(verbosity, 2, 'Calculating filter array...')
        filter_array = self.create_filter_array(freq_array)
        del freq_array
        return mrc, filter_array, apix, is_single_image, is_image_stack, is_single_volume, is_volume_stack

    def copy_file_init_from_other(self, dw):
        copy_attributes = ['filter_array', 'apix', 'is_single_image', 'is_image_stack', 'is_single_volume', 'is_volume_stack']
        [setattr(self, attr, getattr(dw, attr)) for attr in copy_attributes]


    def is_stack_or_volume(self, mrc, interpret_as_slices, interpret_as_images, print_messages=True):
        is_single_image = mrc.is_single_image()
        spacegroup_is_image_stack = mrc.is_image_stack()
        spacegroup_is_single_volume = mrc.is_volume()
        is_volume_stack = mrc.is_volume_stack()
        if spacegroup_is_image_stack or spacegroup_is_single_volume:
            mz = int(mrc.header.mz)
            if mz == 1:
                is_image_stack = True
                is_single_volume = False
            else:
                is_image_stack = False
                is_single_volume = True
        else:
            is_image_stack = spacegroup_is_image_stack
            is_single_volume = spacegroup_is_single_volume
        messages = []
        img_type_str = 'UNKNOWN'
        img_type_str = 'a single image' if is_single_image else img_type_str
        img_type_str = 'an image stack' if is_image_stack else img_type_str
        img_type_str = 'a single volume' if is_single_volume else img_type_str
        img_type_str = 'a volume stack' if is_volume_stack else img_type_str
        msg = 'Input image is %s.' % img_type_str
        messages.append(msg)
        if interpret_as_images and is_single_volume:
            msg = 'Interpreting input image as an image stack'
            messages.append(msg)
            is_single_volume = False
            is_image_stack = True
        if interpret_as_slices and is_image_stack:
            msg = 'Interpreting input image as a single volume'
            messages.append(msg)
            is_single_volume = True
            is_image_stack = False
        if print_messages:
            for msg in messages:
                verbosity_print(verbosity, 2, msg)
        return is_single_image, is_image_stack, is_single_volume, is_volume_stack



    def dose_weight_sharpen(self):

        is_stack = self.is_image_stack or self.is_volume_stack
        #str = 'is stack' if is_stack else 'is volume'
        #print(str)
        read_mode = 'r+' if self.in_place else 'r'

        in_mrc = self.mrc if self.mrc != None else mrcfile.open(self.file_path, mode=read_mode)
        self.img = in_mrc.data
        if not self.in_place:
            split_path = os.path.splitext(self.file_path)
            output_file_path = split_path[0]+self.file_append+split_path[1]
            out_mrc = mrcfile.new(output_file_path, data=np.zeros(self.img.shape, dtype='float32'), overwrite=True)
            out_mrc.voxel_size = in_mrc.voxel_size
            out_mrc.set_image_stack() if is_stack else None
        else:
            out_mrc = in_mrc
            out_mrc.data.setflags(write=True)


        #check input matches expectation
        header_apix = in_mrc.voxel_size
        if self.apix != self.multidim_apix_to_single_value(header_apix):
            warning_msg = 'Warning apix from file header (%s) does not match (%s)' % (str(header_apix), self.apix)
            verbosity_print(verbosity, 2, warning_msg)

        inmrc_is_single_image, inmrc_is_image_stack, inmrc_is_single_volume, inmrc_is_volume_stack = self.is_stack_or_volume(in_mrc, self.interpret_as_slices, self.interpret_as_images, print_messages=False)
        if self.is_volume_stack and (self.is_volume_stack != inmrc_is_volume_stack):
            warning_msg = 'Expecting input mrc to be volume stack but it is not'
            verbosity_print(verbosity, 1, warning_msg)
        if self.is_single_image and (self.is_single_image != inmrc_is_single_image):
            warning_msg = 'Expecting input mrc to be single image but it is not'
            verbosity_print(verbosity, 1, warning_msg)
        if self.is_single_volume and (self.is_single_volume != inmrc_is_single_volume):
            warning_msg = 'Expecting input mrc to be single volume but it is not'
            verbosity_print(verbosity, 1, warning_msg)
        if self.is_image_stack and (self.is_image_stack != inmrc_is_image_stack):
            warning_msg = 'Expecting input mrc to be image stack but it is not'
            verbosity_print(verbosity, 1, warning_msg)


        expected_shape = self.filter_array.shape
        if not is_stack:
            image_shape = self.img.shape
            self.img = [self.img]
        else:
            number_of_stacked_images = self.img.shape[0]
            image_shape = self.img.shape[1:]

        if image_shape != expected_shape:
            error_message = 'Image %s is not the expected shape %s. Skipping...' % (self.file_path, str(expected_shape))
            verbosity_print(verbosity, 1, error_message)
            return

        if self.plot_filters != []:
            fig = plt.figure()

        for i, image in enumerate(self.img):
            if is_stack:
                verbosity_print(verbosity, 3, 'Reading image %d of %d...' % (i + 1, number_of_stacked_images))
            if i in self.plot_filters:
                surf = self.prepare_surf_plot(fig, self.filter_array)
            verbosity_print(verbosity, 3, 'Applying filter...')
            filter_array = self.fft_shift_filter(self.filter_array)
            filtered_image = self.overlay_filter(image, filter_array).astype('float32')
            verbosity_print(verbosity, 3, 'Image filtered.')
            if is_stack:
                out_mrc.data[i,...] = filtered_image
        if not is_stack:
            out_mrc.set_data(filtered_image)
        verbosity_print(verbosity, 2, 'Saving stack ...')

        out_mrc.flush()
        out_mrc.close()
        if not self.in_place:
            in_mrc.close()
        verbosity_print(verbosity, 2, 'Image saved.')
        if self.plot_filters != []:
            plt.show()

    def multidim_apix_to_single_value(self, multidim_apix):
        if type(multidim_apix) == np.recarray:
            apix = float(multidim_apix['x'])
        return apix

    def prepare_surf_plot(self, fig, filter_array):
        # plt.imshow(filter_array, cmap='gray')
        # cbar = plt.colorbar()
        scale_factor = 1
        scale_slice = slice(None, None, scale_factor)
        scale_slices = tuple([scale_slice for j in range(0, filter_array.ndim)])
        binned_filter_array = self.filter_array[scale_slices]  # a crude resampling for the plot.
        is_3d = True if filter_array.ndim > 2 else False
        xindex = 2 if is_3d else 1
        yindex = 1 if is_3d else 0
        zindex = 0 if is_3d else None
        x = np.arange(0, binned_filter_array.shape[xindex])
        y = np.arange(0, binned_filter_array.shape[yindex])
        X, Y = np.meshgrid(x, y)
        ax = fig.gca(projection='3d')
        if is_3d:
            binned_filter_array = binned_filter_array[int(binned_filter_array.shape[0] / 2.)]
        surf = ax.plot_surface(X, Y, binned_filter_array)
        return surf


    def create_frequency_array(self, shape, apix):
        is_3d = True if len(shape)>2 else False
        xindex = 2 if is_3d else 1
        yindex = 1 if is_3d else 0
        zindex = 0 if is_3d else None
        xsize = shape[xindex]
        ysize = shape[yindex]
        zsize = shape[zindex] if is_3d else None
        xcen = (xsize / 2)  # Brigg's center for array is half the image size +1 pix. I removed the 1 as didn't match with fft.
        ycen = (ysize / 2)
        zcen = (zsize / 2) if is_3d else None
        xrstep = 1 / (xsize * apix)  # reciprocal pixel size
        yrstep = 1 / (ysize * apix)
        zrstep = 1 / (zsize * apix) if is_3d else None
        # iterate over numpy array
        if is_3d:
            x,y,z = np.mgrid[0:zsize, 0:ysize, 0:xsize]
            freq_array = np.sqrt(((xrstep * (x - xcen)) ** 2) + ((yrstep * (y - ycen)) ** 2) + ((zrstep * (z - zcen)) ** 2))
        else:
            x,y = np.mgrid[0:ysize, 0:xsize]
            freq_array = np.sqrt(((xrstep * (x - xcen)) ** 2) + ((yrstep * (y - ycen)) ** 2))
        return freq_array

    def return_central_slice(self, image):
        slices = [slice(int(dim_size / 2.), None, int(dim_size / 2.)) if i != len(image.shape) - 1 else slice(int(dim_size / 2.),None, None) for i, dim_size in enumerate(image.shape)]
        return image[slices]

    def create_filter_array(self, freq_array):
        np.seterr(divide='ignore', invalid='ignore')  # avoids a zero divide error printed in output.
        t = np.divide(-1, np.multiply(np.add(np.multiply(np.power(freq_array, self.b), self.a), self.c), 2))
        q = np.divide(np.multiply(np.exp(np.multiply(t, self.dose_per_tilt+self.pre_dose)), np.subtract(1,  np.exp(np.multiply(t, self.dose_per_tilt*self.number_of_tilts)))),          np.multiply(self.number_of_tilts, np.subtract(1, np.exp(np.multiply(t,self.dose_per_tilt)))))
        q[np.where(np.isnan(q))] = 1
        q = (1 / q)
        return q



    def fft_shift_filter(self, filter_array):
        #axes = (0,1,2) if filter_array.ndim == 3 else (0,1)
        filter_array = np.fft.ifftshift(filter_array)
        return filter_array


    def overlay_filter(self, img, filter_array):
        if img.ndim == 3:
            is_3d = True
            axes = (0,1,2)
        else:
            is_3d = False
            axes = (0, 1)
        if use_pfftw:
            if is_3d:
                fft_func = pyfftw.interfaces.numpy_fft.fftn
                ifft_func = pyfftw.interfaces.numpy_fft.ifftn
            else:
                fft_func = pyfftw.interfaces.numpy_fft.fft2
                ifft_func = pyfftw.interfaces.numpy_fft.ifft2
            pyfftw.interfaces.cache.enable()
        else:
            if is_3d:
                fft_func = np.fft.fftn
                ifft_func = np.fft.ifftn
            else:
                fft_func = np.fft.fft2
                ifft_func = np.fft.ifft2
        filtered_img = np.real(ifft_func(np.multiply(fft_func(img, axes=axes), filter_array), axes=axes)).astype('float32')
        if use_pfftw:
            pyfftw.interfaces.cache.disable()
        return filtered_img



def dose_weight_sharpen(input_files, number_of_tilts, dose_per_tilt, pre_dose, all_inputs_equal, in_place, file_append, pixel_size, interpret_as_slices, interpret_as_images, plot_filters, do_not_do_dose_weighting, verbosity_level=default_verbosity_level):
    global verbosity
    verbosity = verbosity_level
    verbosity_print(verbosity, 3, 'Checking input files...')
    for input_file in input_files:
        if not os.path.isfile(input_file):
            raise IOError('Input file %s does not exist' % input_file)


    precalculate_arrays = True if all_inputs_equal and len(input_files) > 1 else False
    if precalculate_arrays:
        verbosity_print(verbosity, 1, 'Precalculating arrays...')
        precalculated_dw = DoseWeightSharpen(input_files[0], dose_per_tilt, pre_dose, number_of_tilts, pixel_size, interpret_as_slices, interpret_as_images, file_append, in_place, plot_filters=plot_filters)
    else:
        precalculated_dw = None
    number_of_files = len(input_files)
    verbosity_print(verbosity, 1, 'Starting dw sharpening%s...' % ' (%d files)' % number_of_files if number_of_files > 1 else '')
    for i, input_file in enumerate(input_files):
        verbosity_print(verbosity, 2, 'File %d of %d' % (i+1, number_of_files)) if number_of_files > 1 else None
        dw = DoseWeightSharpen(input_file, dose_per_tilt, pre_dose, number_of_tilts, pixel_size, interpret_as_slices, interpret_as_images, file_append, in_place, plot_filters=plot_filters, copy_file_init_from=precalculated_dw)
        dw.mrc = precalculated_dw.mrc if precalculate_arrays and i == 0 else None
        if do_not_do_dose_weighting == False:
            dw.dose_weight_sharpen()
        else:
            verbosity_print(verbosity, 2, 'Skipping actually doing the dose weighting as --do_not_do_dose_weighting set')
    verbosity_print(verbosity, 1, 'Done.')

def main(input_map,
        number_of_tilts,
        dose_per_tilt,
        pixel_size,
        pre_dose,
        in_place,
        file_append,
        do_not_do_dose_weighting,
        all_inputs_equal,
        interpret_as_slices,
        interpret_as_images,
        verbosity
        ):
    input_maps = sorted(glob.glob(input_map))
    dose_weight_sharpen(input_maps, number_of_tilts, dose_per_tilt, pre_dose, all_inputs_equal, in_place, file_append, pixel_size, interpret_as_slices, interpret_as_images, plot_filters, do_not_do_dose_weighting, verbosity)


if __name__ == "__main__":
    argparser = ArgumentParser()
    args = argparser.parser.parse_args()
    argparser.validate(args)

    main(args.input_map,
        args.number_of_tilts,
        args.dose_per_tilt,
        args.pixel_size,
        args.pre_dose,
        args.in_place,
        args.file_append,
        args.do_not_do_dose_weighting,
        args.all_inputs_equal,
        args.interpret_as_slices,
        args.interpret_as_images,
        args.verbosity
        )




