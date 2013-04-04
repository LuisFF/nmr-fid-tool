package uk.ac.ebi.nmr.fid.io;

import edu.uchc.connjur.spectrumtranslator.bruker.BrukerDataSetReader;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

import java.io.File;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:20
 * To change this template use File | Settings | File Templates.
 */
public class ConnjurProcReader implements ProcReader {

    private BrukerDataSetReader dataSetReader;
    private Acqu acqu;

    public ConnjurProcReader (File file, Acqu acqu){
        dataSetReader = new BrukerDataSetReader();
        this.acqu=acqu;
        try {
            dataSetReader.setDataSource(file);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
    @Override
    public Proc read() throws IOException {
//        Proc proc = new Proc();
        return null;  //TODO change body of implemented methods
    }
}
