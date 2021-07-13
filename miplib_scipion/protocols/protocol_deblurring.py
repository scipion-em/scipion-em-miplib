# **************************************************************************
# *
# * Authors:     Ana Cayuela (acayuela@cnb.csic.es)
# *
# * Natl. Center of Biotechnology (CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
from os.path import exists
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwem.emlib.image import ImageHandler
from pwem.emlib.lib import Image
from pwem.protocols.protocol import EMProtocol
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram

from miplib.psf import psfgen
from miplib.processing.deconvolution import deconvolve
import miplib.ui.cli.miplib_entry_point_options as options
from miplib.data.containers.image import Image as MIPImage
import miplib.analysis.resolution.fourier_ring_correlation as frc
from miplib.data.containers.fourier_correlation_data import FourierCorrelationDataCollection

class MIPLIB_Deblurring(EMProtocol, ProtTomoBase):
    """
       This protocol deblurs a tomogram using the FRC
       """
    _label = 'Tomogram deblurring'

    MODE_DEBLURRING = 0
    MODE_DENOISING = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTomogram', params.PointerParam, pointerClass='Tomogram', label='Tomogram',
                      allowsNull=False)
        form.addParam('Niter', params.IntParam, label='Number of iterations', default=10,
                      help='Number of iterations of the Richardson-Lucy deconvolution')
        form.addParam('mode', params.EnumParam, choices=['Deblurring', 'Denoising'], default=0, label='Operation')
        form.addParam('lambda', params.FloatParam, label='Total Variation weight', default=0.0,
                      expertLevel=LEVEL_ADVANCED, condition='mode==MODE_DEBLURRING',
                      help='Weight of the penalization by Total Variation')
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('deblurStep', self.inputTomogram.getObjId())
        self._insertFunctionStep('createOutputStep')

    def deblurStep(self, objId):
        # Grab volume from input
        tomogram = self.inputTomogram.get()

        if self.mode.get()==self.MODE_DEBLURRING:
            args_list = ("image psf"
                         "moving_image"
                         " --max-nof-iterations={}  --first-estimate=image "
                         " --blocks=1 --pad=30 --resolution-threshold-criterion=fixed "
                         " --tv-lambda=0 --bin-delta=1  --frc-curve-fit-type=smooth-spline").format(self.Niter.get()).split()
        elif self.mode.get()==self.MODE_DENOISING:
            args_list = []
        args = options.get_deconvolve_script_options(args_list)
        Xdim, Ydim, Zdim = tomogram.getDim()
        ih = ImageHandler()
        tomogramImage = ih.read(tomogram)
        tomogramData = tomogramImage.getData()
        spacing = tomogram.getSamplingRate()
        # writer = imwrap.TiffImageWriter(self._getTmpPath())

        for z in range(Zdim):
            sliceZ = tomogramData[z,:,:]
            m = sliceZ.min()
            M = sliceZ.max()
            sliceZ = 16536*(sliceZ-m)/(M-m)
            sliceZMIP = MIPImage(sliceZ, [spacing, spacing])
            # print(sliceZMIP.shape)
            # print("Max: ", np.max(sliceZMIP))
            # print("Min: ", np.min(sliceZMIP))
            # save = Image()
            # save.setData(sliceZMIP)
            # save.write(self._getExtraPath("image%d.xmp"%z))
            frc_results = FourierCorrelationDataCollection()
            try:
                frc_results[0] = frc.calculate_single_image_frc(sliceZMIP, args)
                fwhm = [frc_results[0].resolution['resolution'], ] * 2
                print("Slice %d, FWHM=%s"%(z,str(fwhm)))
                psf_generator = psfgen.PsfFromFwhm(fwhm)
                psf = psf_generator.xy()
                task = deconvolve.DeconvolutionRL(sliceZMIP, psf, None, args)
                task.execute()
                sliceZMIP = task.get_result()
                # print("Max: ", np.max(sliceZMIP))
                # print("Min: ", np.min(sliceZMIP))
                mp=sliceZMIP.min()
                Mp=sliceZMIP.max()
                tomogramData[z,:,:] = m+(M-m)*(sliceZMIP-mp)/(Mp-mp)
            except:
                print("Skipping slice %d"%z)

        tomogramImage.setData(tomogramData)
        tomogramImage.write(self._getPath("enhancedTomogram.mrc"))

    def createOutputStep(self):
        fnVol = self._getPath('enhancedTomogram.mrc')
        if exists(fnVol):
            tomogram = Tomogram()
            tomogram.setFileName(fnVol)
            tomogram.setSamplingRate(self.inputTomogram.get().getSamplingRate())
            self._defineOutputs(outputTomogram=tomogram)
            self._defineSourceRelation(self.inputTomogram.get(), tomogram)

    # --------------------------- INFO functions -----------------------------------
    def _methods(self):
        methods = []

        if self.isFinished():
            methodStr = "The tomogram %s was deblurred with its own FRC [[Koho2019]] with %d Richardson-Lucy iterations " %\
                        (self.getObjectTag('inputVolume'),self.Niter)
            methods.append(methodStr)
            methodStr += " The tomogram %s was created." % (self.getObjectTag('outputTomogram'))
        return methods

    def _citations(self):
        return ['Koho2019']