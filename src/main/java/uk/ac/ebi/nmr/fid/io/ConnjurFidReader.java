package uk.ac.ebi.nmr.fid.io;

import edu.uchc.connjur.spectrumtranslator.DataControl;
import edu.uchc.connjur.spectrumtranslator.SwDataSet;
import edu.uchc.connjur.spectrumtranslator.bruker.BrukerDataSetReader;
import uk.ac.ebi.nmr.fid.Fid;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 21/03/2013
 * Time: 12:00
 * To change this template use File | Settings | File Templates.
 */
public class ConnjurFidReader implements FidReader {


    private BrukerDataSetReader dataSetReader;

    public ConnjurFidReader(File file) {
        dataSetReader = new BrukerDataSetReader();
        try {
            dataSetReader.setDataSource(file);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    @Override
    public Fid read() throws Exception {
        
        Fid fid = null;

        DataControl dataControl = new DataControl();
        SwDataSet dataSet = dataSetReader.read(null,dataControl);
        // extract the real and imaginary values of the FID
        Iterator<edu.uchc.connjur.core.Fid> fidIterator = dataSet.fidIterator();
        edu.uchc.connjur.core.Fid fidReal = fidIterator.next();
        edu.uchc.connjur.core.Fid fidIm = fidIterator.next();

        //put them together because the FID only accepts a list of alternating real + imaginary numbers
        int[] fidCombined = new int[fidReal.getSize()+fidIm.getSize()];
        for(int i = 0; i < fidReal.getSize();i++){
            fidCombined[i*2]=(int) Math.round(fidReal.getValue(i));
            fidCombined[i*2+1]=(int) Math.round(fidIm.getValue(i));
        }

        return new Fid(fidCombined);  //To change body of implemented methods use File | Settings | File Templates.
    }
}
