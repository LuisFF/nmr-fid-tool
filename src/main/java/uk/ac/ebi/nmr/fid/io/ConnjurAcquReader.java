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

package uk.ac.ebi.nmr.fid.io;

import edu.uchc.connjur.core.DimInfo;
import edu.uchc.connjur.spectrumtranslator.DataControl;
import edu.uchc.connjur.spectrumtranslator.NumericStorageType;
import edu.uchc.connjur.spectrumtranslator.SwDataSet;
import edu.uchc.connjur.spectrumtranslator.bruker.BrukerDataSetReader;
import uk.ac.ebi.nmr.fid.Acqu;

import java.io.File;
import java.io.IOException;

/**
 * Connjur integration to read acquisition parameters
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:25
 *
 */
public class ConnjurAcquReader implements AcquReader {


    private BrukerDataSetReader dataSetReader;
    private Acqu.Spectrometer spectrometer;

    public ConnjurAcquReader(File file,Acqu.Spectrometer spectrometer) {
        this.dataSetReader = new BrukerDataSetReader();
        this.spectrometer=spectrometer;
        try {
            dataSetReader.setDataSource(file);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }
    @Override
    public Acqu read() throws Exception {
        Acqu acqu = new Acqu(spectrometer);

        DataControl dataControl = new DataControl();
        SwDataSet dataSet = dataSetReader.read(null,dataControl);

        DimInfo dimInfo = dataSet.getDataByDim(0);
        acqu.setAquiredPoints(dimInfo.getTotalPoints());
        // somehow the conversion from float to double adds some decimal places..
        acqu.setSpectralFrequency(Double.parseDouble(Float.toString(dimInfo.getSpectralFrequency())));
        acqu.setSpectralWidth(Double.parseDouble(Float.toString(dimInfo.getSweepWidthPPM())));
        // spectral width in Hz
//        System.out.println("sweep width "+dimInfo.getSweepWidth());

        acqu.set32Bit(dataSet.getSourceDataType().equals(NumericStorageType.INT32) ||
                dataSet.getSourceDataType().equals(NumericStorageType.IEEE_FLOAT32));

        acqu.setDspGroupDelay(Double.parseDouble(dataSet.supplementalMetadata().get("GroupDelay")));




        // Additional data but not from acqu file...
        // from the dataSet
        for (String key : dataSet.supplementalMetadata().keySet()){
            System.out.println(key+" "+ dataSet.supplementalMetadata().get(key));
        }
        System.out.println(dataSet.sizeAndDescription().getDescription());
        // from the dimInfo
        System.out.println("Analysis domain " + dimInfo.getAnalysisDomain().name());
        System.out.println("First Order Phase " + dimInfo.getFirstOrderPhase());
        System.out.println("Zero Order Phase " + dimInfo.getZeroOrderPhase());
        System.out.println("Carrier PPM " + dimInfo.getCarrierPPM());
        System.out.println("Data mode - Numbers per data point " + dimInfo.getDataMode().NUMBERS_PER_DATA_POINT);
        System.out.println("Data mode - code " + dimInfo.getDataMode().CODE);
        System.out.println("first scale point " + dimInfo.getFirstScalePoint());

        return acqu;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
