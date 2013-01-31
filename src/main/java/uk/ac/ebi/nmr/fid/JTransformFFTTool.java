/**
 * JTransformFFTTool.java
 *
 * 2013.01.31
 *
 * This file is part of the CheMet library
 * 
 * The CheMet library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * CheMet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with CheMet.  If not, see <http://www.gnu.org/licenses/>.
 */

package uk.ac.ebi.nmr.fid;


import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import org.apache.log4j.Logger;

/**
 * @name    JTransformFFTTool
 * @date    2013.01.31
 * @version $Rev$ : Last Changed $Date$
 * @author  pmoreno
 * @author  $Author$ (this version)
 * @brief   ...class description...
 *
 */
public class JTransformFFTTool extends AbstractFastFourierTransformTool implements FastFourierTransformTool{

    private static final Logger LOGGER = Logger.getLogger( JTransformFFTTool.class );
    
    public JTransformFFTTool(Fid fid, Acqu acquisition) throws Exception{                
        this.processing = new Proc(acquisition);     
        this.acquisition = acquisition;
        this.fid = fid;
    }

    public JTransformFFTTool(Fid fid, Acqu acquisition, Proc processing) {
        this.fid=fid;
        this.acquisition=acquisition;
        this.processing=processing;
    }

    @Override
    double[] implementedFFT(double[] apodizedData) {
        DoubleFFT_1D dfftd = new DoubleFFT_1D(apodizedData.length/2);
        dfftd.complexForward(apodizedData);
        return apodizedData;
    }

    public void fft(boolean isForward) {
        throw new UnsupportedOperationException("Not supported yet.");
    }


}
