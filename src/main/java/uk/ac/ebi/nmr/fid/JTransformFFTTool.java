/*
 * Copyright (c) 2013. EMBL, European Bioinformatics Institute
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

package uk.ac.ebi.nmr.fid;


import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import org.apache.log4j.Logger;

/**
 * @name    JTransformFFTTool
 * @date    2013.01.31
 * @version $Rev$ : Last Changed $Date$
 *
 * @author  Luis F. de Figueiredo
 * @author  pmoreno
 * @author  $Author$ (this version)
 * @brief   ...class description...
 *
 */
public class JTransformFFTTool extends AbstractFastFourierTransformTool implements FastFourierTransformTool{

    private static final Logger LOGGER = Logger.getLogger( JTransformFFTTool.class );

    public JTransformFFTTool(Spectrum spectrum) {
        this.fid=spectrum;

    }

    @Override
    double[] implementedFFT(double[] apodizedData) {
        DoubleFFT_1D dfftd = new DoubleFFT_1D(apodizedData.length/2);
        dfftd.complexForward(apodizedData);
        return getRealPart(apodizedData);
    }

    public void fft(boolean isForward) {
        throw new UnsupportedOperationException("Not supported yet.");
    }


    @Override
    public Proc getProcessing() {
        return fid.getProc();  //To change body of implemented methods use File | Settings | File Templates.
    }
}
