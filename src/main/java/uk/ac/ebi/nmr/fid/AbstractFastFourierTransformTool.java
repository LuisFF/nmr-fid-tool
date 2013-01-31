/**
 * AbstractFastFourierTransformTool.java
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


import org.apache.log4j.Logger;
import uk.ac.ebi.nmr.fid.tools.ApodizationTool;

/**
 * @name    AbstractFastFourierTransformTool
 * @date    2013.01.31
 * @version $Rev$ : Last Changed $Date$
 * @author  pmoreno
 * @author  $Author$ (this version)
 * @brief   Abstract Fast Fourier Transform Tool class to be extended by all FFT implementations.
 *
 */
public abstract class AbstractFastFourierTransformTool {

    Acqu acquisition;
    double[] data; //proc_buffer   apodization::transform::do_fft
    //proc_buffer   apodization::transform::do_fft
    Fid fid;
    Proc processing;

    public AbstractFastFourierTransformTool() {
    }

    void applyRedfieldOnSequentialData(Acqu acquisition, Proc processing) {
        // for sequential data apply Redfield trick: multiply data by 1 -1 -1 1
        if (acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i = 0; i < processing.getTdEffective(); i += 4) {
                data[i + 1] = -(data[i + 1]);
                data[i + 2] = -(data[i + 2]);
            }
        }
    }

    public double[] computeFFT() {
        initDataFormat(acquisition, processing);
        shiftData(processing, acquisition, fid);
        applyRedfieldOnSequentialData(acquisition, processing);
        tdRelatedNonUnderstoodRearrangementForSequential(processing, acquisition);
        int signals;
        processing.setLineBroadning(0.3);
        /// applyWindowFunctions //window function need to be applied before FT
        ApodizationTool apodization = new ApodizationTool(data, acquisition.getAcquisitionMode(), processing);
        double[] apodizedData = null;
        apodization.exponential();
        apodizedData = apodization.getSpectrum();
        return implementedFFT(apodizedData);
    }

    public Acqu getAcquisition() {
        return acquisition;
    }

    public double[] getData() {
        return data;
    }

    public Fid getFid() {
        return fid;
    }

    public Proc getProcessing() {
        return processing;
    }

    /**
     * This is where the implementation of each fft package go.
     * @param apodizedData
     * @return
     */
    abstract double[] implementedFFT(double[] apodizedData);

    void initDataFormat(Acqu acquisition, Proc processing) {
        // instanciating the array where the fourier transformed spectra will be stored....
        switch (acquisition.getAcquisitionMode()) {
            case DISP:
            case SIMULTANIOUS:
                data = new double[2 * processing.getTransformSize()]; // allocate space for processing
                // allocate space for processing
                break;
            case SEQUENTIAL:
                data = new double[4 * processing.getTransformSize()]; // allocate space for processing
                // allocate space for processing
                break;
            default:
                break;
        }
    }

    public void setData(double[] data) {
        this.data = data;
    }

    public void setData(int index, double point) {
        this.data[index] = point;
    }

    void shiftData(Proc processing, Acqu acquisition, Fid fid) {
        // perform a left or right shift of the data (ignoring the corresponding portion of the data)
        // the code from cuteNMR was simplified
        for (int i = 0; i < (processing.getTdEffective() - Math.abs(processing.getShift())); i++) {
            int dataIndex = (processing.getShift() >= 0) ? i : i - processing.getShift(); // or shift the placement of the data to the right
            // or shift the placement of the data to the right
            int fidIndex = (processing.getShift() >= 0) ? i + processing.getShift() : i; // start in the correct order
            // start in the correct order
            switch (acquisition.getFidType()) {
                case INT32:
                    data[dataIndex] = (double) fid.getData()[fidIndex];
                    break;
                case DOUBLE:
                    data[dataIndex] = fid.getData()[fidIndex];
                    break;
                case FLOAT:
                    break;
                case INT16:
                    break;
                default:
                    break;
            }
        }
    }

    void tdRelatedNonUnderstoodRearrangementForSequential(Proc processing, Acqu acquisition) {
        //// try to understand this bit of code!!!!
        //nonsense case if SI is set to be less than TD/2
        int td = (processing.getTdEffective() > 2 * processing.getTransformSize()) ? 2 * processing.getTransformSize() : processing.getTdEffective();
        // move the data from position i to position 2*i, why?
        if (acquisition.getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i = td - 1; i > 0; i--) {
                data[2 * i] = data[i];
                data[i] = 0;
            }
        }
        /////////////////////////////////////////////
    }


}
