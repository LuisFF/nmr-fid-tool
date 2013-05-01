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


/**
 * @name    AbstractFastFourierTransformTool
 * @date    2013.01.31
 * @version $Rev$ : Last Changed $Date$
 * @author  Luis F. de Figueiredo
 * @author  pmoreno
 * @author  $Author$ (this version)
 * @brief   Abstract Fast Fourier Transform Tool class to be extended by all FFT implementations.
 *
 */
public abstract class AbstractFastFourierTransformTool {


    double[] data; //proc_buffer   apodization::transform::do_fft
    //proc_buffer   apodization::transform::do_fft
    Spectrum fid;


    public AbstractFastFourierTransformTool() {
    }

    void applyRedfieldOnSequentialData() {
        // for sequential data apply Redfield trick: multiply data by 1 -1 -1 1
        if (fid.getAcqu().getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i = 0; i < fid.getProc().getTdEffective(); i += 4) {
                data[i + 1] = -(data[i + 1]);
                data[i + 2] = -(data[i + 2]);
            }
        }
    }

    public double[] computeFFT() {
        initDataFormat();
        shiftData();
        applyRedfieldOnSequentialData();
        tdRelatedNonUnderstoodRearrangementForSequential();
        int signals;
        fid.getProc().setLineBroadening(0.3);
        /// applyWindowFunctions //window function need to be applied before FT
        //TODO adapt the FFT to the new object Spectrum
//        Apodizator apodization = new ExponentialApodizator(data, acquisition.getAcquisitionMode(), processing);
        double[] apodizedData = null;

        try {
//            apodizedData = apodization.calculate();
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return implementedFFT(apodizedData);
    }
    
    double[] getRealPart(double[] data) {
        double[] realPart = new double[data.length/2];
        int realIndex = 0;
        for (int i = 0; i < data.length; i+=2) {
            realPart[realIndex++] = data[i];
        }
        return realPart;
    }



    public double[] getData() {
        return data;
    }

    public double[] getFid() {
        return fid.getFid();
    }



    /**
     * This is where the implementation of each fft package go.
     * @param apodizedData
     * @return
     */
    abstract double[] implementedFFT(double[] apodizedData);

    void initDataFormat() {
        // instanciating the array where the fourier transformed spectra will be stored....
        switch (fid.getAcqu().getAcquisitionMode()) {
            case DISP:
            case SIMULTANIOUS:
                data = new double[2 * fid.getProc().getTransformSize()]; // allocate space for processing
                // allocate space for processing
                break;
            case SEQUENTIAL:
                data = new double[4 * fid.getProc().getTransformSize()]; // allocate space for processing
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

    void shiftData() {
        // perform a left or right shift of the data (ignoring the corresponding portion of the data)
        // the code from cuteNMR was simplified
        for (int i = 0; i < (fid.getProc().getTdEffective() - Math.abs(fid.getProc().getShift())); i++) {
            int dataIndex = (fid.getProc().getShift() >= 0) ? i : i - fid.getProc().getShift(); // or shift the placement of the data to the right
            // or shift the placement of the data to the right
            int fidIndex = (fid.getProc().getShift() >= 0) ? i + fid.getProc().getShift() : i; // start in the correct order
            // start in the correct order
            switch (fid.getAcqu().getFidType()) {
                case INT32:
                    data[dataIndex] = (double) fid.getFid()[fidIndex];
                    break;
                case DOUBLE:
                    data[dataIndex] = fid.getFid()[fidIndex];
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

    void tdRelatedNonUnderstoodRearrangementForSequential() {
        //// try to understand this bit of code!!!!
        //nonsense case if SI is set to be less than TD/2
        int td = (fid.getProc().getTdEffective() > 2 * fid.getProc().getTransformSize()) ?
                2 * fid.getProc().getTransformSize() : fid.getProc().getTdEffective();
        // move the data from position i to position 2*i, why?
        if (fid.getAcqu().getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i = td - 1; i > 0; i--) {
                data[2 * i] = data[i];
                data[i] = 0;
            }
        }
        /////////////////////////////////////////////
    }

    public Proc getProcessing(){
        return fid.getProc();
    }

    public Acqu getAcquisition(){
        return fid.getAcqu();
    }
}
