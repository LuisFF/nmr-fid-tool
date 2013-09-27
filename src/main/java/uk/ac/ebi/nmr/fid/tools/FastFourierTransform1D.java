/*
 * Copyright (c) 2013 EMBL, European Bioinformatics Institute.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package uk.ac.ebi.nmr.fid.tools;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import uk.ac.ebi.nmr.fid.Spectrum;

/**
 *
 * @author  Luis F. de Figueiredo
 * User: ldpf
 * Date: 24/05/2013
 * Time: 16:48
 * To change this template use File | Settings | File Templates.
 */
public class FastFourierTransform1D implements FastFourierTransform {
    
    Spectrum spectrum;

    public FastFourierTransform1D(Spectrum spectrum) {
        this.spectrum = spectrum;

    }

    @Override
    public Spectrum computeFFT(int offset) {
        double [] fid = spectrum.getFid();
        double [] realPart;
        double [] imaginaryPart;

        // run the FFT
        DoubleFFT_1D fftd = new DoubleFFT_1D(fid.length/2-offset);
        fftd.complexForward(fid,offset);

        double [] realChannel = new double[fid.length/2];
        double [] imagChannel = new double[fid.length/2];
        double maxFFTdata=0;
        double minFFTdata=0;

        // calculate the maximum and minimum to normalize the intentisity
        for(int i=0; i<fid.length; i++){
            minFFTdata=Math.min(minFFTdata,fid[i]);
            maxFFTdata=Math.max(maxFFTdata,fid[i]);
        }

        // extract the real an imaginary parts and normalize
        for(int i =0; i< fid.length; i+=2){
            realChannel[i/2]=fid[i]/(maxFFTdata-minFFTdata);
            imagChannel[i/2]=fid[i+1]/(maxFFTdata-minFFTdata);
        }
        // flip each half...
        double [] tmpReal = new double[realChannel.length];
        double [] tmpImag = new double[imagChannel.length];
        for(int i =0; i< realChannel.length/2; i++){
            tmpReal[realChannel.length/2-i]=realChannel[i];
            tmpImag[realChannel.length/2-i]=imagChannel[i];
            tmpReal[realChannel.length-i-1]=realChannel[realChannel.length/2+i-1];
            tmpImag[imagChannel.length-i-1]=imagChannel[realChannel.length/2+i-1];
        }
        spectrum.setRealChannelData(tmpReal);
        spectrum.setImaginaryChannelData(tmpImag);
        return spectrum;
    }

    @Override
    public Spectrum computeFFT() {
        return computeFFT(0);
    }
}
