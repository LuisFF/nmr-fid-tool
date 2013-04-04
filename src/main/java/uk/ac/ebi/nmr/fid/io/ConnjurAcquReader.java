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
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:25
 * To change this template use File | Settings | File Templates.
 */
public class ConnjurAcquReader implements AcquReader {


    private BrukerDataSetReader dataSetReader;

    public ConnjurAcquReader(File file) {
        dataSetReader = new BrukerDataSetReader();
        try {
            dataSetReader.setDataSource(file);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }
    @Override
    public Acqu read() throws Exception {
        Acqu acqu = new Acqu();

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
