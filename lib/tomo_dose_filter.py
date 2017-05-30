#!/usr/bin/env python

import mrcfile
import numpy as np
import math
import sys
import os
import glob
import argparse


class ArgumentParser():
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                              description="Applies dose weighting to a tilt series stack.")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        addr('-i', '--tilt_series', required=True,
             help='A wild card expression (IN QUOTES!) to the tilt series stacks. (Images should be ordered correctly in the stack -ve to +ve)')
        addr('--tilt_scheme',
             help='A tilt scheme used from the following list; continuous_positive, continuous_negative, bidirectional_positive, bidirectional_negative')
        addr('-dose', '--dose_per_tilt', type=float, help='Dose applied per tilt image. (in e-/A^2)')
        addr('--min_angle', type=float, help='The minimum (most negative) angle used in the tilt series')
        addr('--angle_step', type=float, help='The angular step between tilt images')
        addr('-apix', '--pixel_size', required=True, type=float, help='The pixel size of the images in angstrom')
        add('--pre_dose', default=0, type=float, help='Initial dose before tilt series collected.')
        add('--file_append', default='dw', type=str, help='String to append to the end of the file.')
        add('--do_not_do_dose_weighting', action='store_true',
            help='Set this to just check the files and not actually apply dose weighting.')
        add('--custom_dose_series', default=None, type=str,
            help='A custom comma delimited list of the doses to apply. This overwites the --dose_per_tilt value given above (must be in the same order as the images. eg 2,4,6,8,10,12,14,16,18)')

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
        if args.tilt_series == None and glob.glob(args.tilt_series) == []:
            self.error('Error: No files found.' % (args.tilt_series))
            sys.exit(2)

        if args.tilt_scheme not in accepted_tilt_schemes and args.tilt_scheme != None:
            self.error('Error: Tilt scheme not supported.' % (args.tilt_scheme))
            sys.exit(2)

        if args.custom_dose_series == None:
            required_args = ('tilt_scheme', 'dose_per_tilt', 'min_angle', 'angle_step')
            error_msgs = []
            for arg in required_args:
                if getattr(args, arg) == None:
                    error_msgs.append('required argument --%s not given' % (arg))
            print(error_msgs)
            if error_msgs != []:
                self.error(*error_msgs)


# Nitpicky details
starting_tilt_angle = 0  # only needed for the bidirectional tilt schemes
plot_filters = []  # [0,10,29] #[0,1,2]
keep_header_apix = True #use the original pixel size in the header of the output file. (apix is still used for the dose weighting). This avoids mismatches in pixel size between the input and output stacks.

####
if plot_filters != []:  # only use if matplotlib available
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

accepted_tilt_schemes = ['continuous_positive', 'continuous_negative', 'bidirectional_positive',
                         'bidirectional_negative']


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
    else:
        print('Tilt scheme not supported')
        raise ValueError
    return tilt_order


class DoseWeight:
    def __init__(self, images, doses, apix, file_append, plot_filters=[0]):
        self.a = 0.245
        self.b = -1.665
        self.c = 2.81
        self.files = images  # list of 2D images or a single mrc stack of images
        self.is_stack = True if type(self.files) != list else False
        if self.is_stack:
            self.images, self.header_apix = self.read_image(self.files)
            self.number_of_files = self.images.shape[0]
            self.filtered_images = np.empty(self.images.shape, dtype='float32')
        else:
            self.images = self.files
            self.number_of_files = len(self.images)
        self.doses = doses
        self.apix = apix
        self.file_append = file_append
        self.plot_filters = plot_filters  # list of indices to plot from the images list
        self.plot_filters = [x for x in self.plot_filters if
                             x >= 0 and x < self.number_of_files]  # remove nonsense values

    # self.dose_weight()

    def dose_weight(self):
        np.seterr(divide='ignore')  # avoids a zero divide error printed in output.
        print('Pre calculating frequency array ...')
        if not self.is_stack:
            img, self.header_apix = self.read_image(self.images[0])
        else:
            img = self.images[0]
        shape = img.shape
        if len(shape) != 2:
            print('Images must be 2D. Quitting...')
            sys.exit(2)
        freq_array = self.create_frequency_array(shape, self.apix)
        print('Frequency array created.')
        if self.plot_filters != []:
            fig = plt.figure()
        for i, (image, dose) in enumerate(zip(self.images, self.doses)):
            print('Reading image %d of %d...' % (i + 1, self.number_of_files))
            if not self.is_stack:
                img, header_apix = self.read_image(image)
            else:
                img = image
            if img.shape != shape:
                print('Image %d is not the expected size. Skipping...' % (i + 1))
                continue
            print('Creating filter array...')
            filter_array = self.create_filter_array(dose, freq_array, self.a, self.b, self.c)
            print('Filter array created.')
            if i in self.plot_filters:
                # plt.imshow(filter_array, cmap='gray')
                # cbar = plt.colorbar()
                scale_factor = 16
                binned_filter_array = filter_array[::scale_factor, ::scale_factor]  # a crude resampling for the plot.
                x = np.arange(0, binned_filter_array.shape[1])
                y = np.arange(0, binned_filter_array.shape[0])
                X, Y = np.meshgrid(x, y)
                ax = fig.gca(projection='3d')
                surf = ax.plot_surface(X, Y, binned_filter_array)
            print('Overlaying filter...')
            filtered_image = self.overlay_filter(img, filter_array)
            print('Image filtered.')
            if not self.is_stack:
                print('Saving image ...')
                filename, file_extension = os.path.splitext(image)
                outfile = filename + '_' + self.file_append + file_extension
                out_apix = self.header_apix if keep_header_apix else self.apix
                self.write_image(filtered_image, outfile, out_apix)
                print('Image saved.')
            else:
                self.filtered_images[i] = filtered_image
        if self.is_stack:
            print('Saving stack ...')
            filename, file_extension = os.path.splitext(self.files)
            outfile = filename + '_' + self.file_append + file_extension
            out_apix = self.header_apix if keep_header_apix else self.apix
            self.write_image(self.filtered_images, outfile, out_apix)
            print('Stack saved.')
        if self.plot_filters != []:
            plt.show()

    def read_image(self, image):
        with mrcfile.open(image) as mrc:
            return mrc.data, mrc.voxel_size

    def write_image(self, image, path, apix=1):
        with mrcfile.new(path, overwrite=True) as mrc:
            mrc.set_data(image)
            mrc.set_image_stack()
            mrc.voxel_size = apix
            mrc.close()

    def create_frequency_array(self, shape, apix):
        freq_array = np.zeros(shape)
        xsize = shape[1]
        ysize = shape[0]
        xcen = (
        xsize / 2)  # Brigg's center for array is half the image size +1 pix. I removed the 1 as didn't match with fft.
        ycen = (ysize / 2)
        xrstep = 1 / (xsize * apix)  # reciprocal pixel size
        yrstep = 1 / (ysize * apix)
        # iterate over numpy array
        it = np.nditer(freq_array, flags=['multi_index'], op_flags=['readwrite'])
        while not it.finished:
            x = it.multi_index[1]
            y = it.multi_index[0]
            d = math.sqrt(((xrstep * (x - xcen)) ** 2) + ((yrstep * (y - ycen)) ** 2))
            it[0] = d
            it.iternext()
        return freq_array

    def create_filter_array(self, dose, freq_array, a, b, c):
        # q = exp((-dose)./(2.*((a.*(freq_array.^b))+c)));
        q = np.exp(np.divide(-dose, np.multiply(np.add(np.multiply(np.power(freq_array, b), a), c), 2)))
        return q

    def overlay_filter(self, img, filter_array):
        fft = np.fft.fftshift(np.fft.fft2(img))
        filtered_fft = np.multiply(fft, filter_array)
        filtered_img = np.absolute(np.fft.ifft2(np.fft.ifftshift(filtered_fft)), dtype='float32')
        return filtered_img


def tilt_series_dose_weight(tilt_series, dose_per_tilt, file_append, apix, plot_filters, tilt_scheme, min_angle,
                            angle_step, do_not_do_dose_weighting, custom_dose_series, pre_dose):
    dw = DoseWeight(tilt_series, [], apix, file_append, plot_filters)
    if custom_dose_series == None:
        total_tilts = dw.number_of_files
        max_angle = min_angle + (angle_step * (total_tilts - 1))
        tilt_list = [min_angle + (angle_step * i) for i in range(0, total_tilts)]
        zero_tilt_index = tilt_list.index(min(tilt_list, key=lambda x: abs(x - starting_tilt_angle)))
        order_list = [tilt_order(tilt_scheme, i + 1, total_tilts, zero_tilt_index) for i in range(0, total_tilts)]
        doses = [(order * dose_per_tilt) + pre_dose for order in order_list]
        print('Tiltseries from %d to %d degrees in steps of %d degrees (%d tilts in total) using a %s tilt scheme.' % (
        min_angle, max_angle, angle_step, total_tilts, tilt_scheme))
    else:
        doses = custom_dose_series.split(',')
        doses = [float(dose) for dose in doses]
        if len(doses) != dw.number_of_files:
            print('Not the correct number of entries in the dose list. Skipping...')
            return
    print('The following doses are used for dose weighting each tilt image: %s' % (str(doses)))
    dw.doses = doses
    if do_not_do_dose_weighting == False:
        dw.dose_weight()
    else:
        print('Skipping actually doing the dose weighting as --do_not_do_dose_weighting set')
        return


def main(tilt_series, dose_per_tilt, file_append, apix, plot_filters, tilt_scheme, min_angle, angle_step,
         do_not_do_dose_weighting, custom_dose_series, pre_dose):
    tilt_series = sorted(glob.glob(tilt_series))
    for stack in tilt_series:
        tilt_series_dose_weight(stack, dose_per_tilt, file_append, apix, plot_filters, tilt_scheme, min_angle,
                                angle_step, do_not_do_dose_weighting, custom_dose_series, pre_dose)


if __name__ == "__main__":
    argparser = ArgumentParser()
    args = argparser.parser.parse_args()
    argparser.validate(args)

    main(args.tilt_series, args.dose_per_tilt, args.file_append, args.pixel_size, plot_filters, args.tilt_scheme,
         args.min_angle, args.angle_step, args.do_not_do_dose_weighting, args.custom_dose_series, args.pre_dose)



