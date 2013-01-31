package uk.ac.ebi.nmr.fid.io;


import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Fid;
import uk.ac.ebi.nmr.fid.LdpfsFourierTransformTool;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 21/09/2012
 * Time: 14:46
 * To change this template use File | Settings | File Templates.
 */
public class BrukerReader {

    private String fidFileName = null;
    private File acquFile = null;
    private File procFile = null;


    public BrukerReader() {
    }

    public BrukerReader(String filename) throws FileNotFoundException {

        this.fidFileName = filename;
        File fidFile = new File(this.fidFileName);
        String workingDIR = fidFile.getParent();
        this.acquFile = new File(workingDIR + "/acqu");
        // TODO make other PROC files available
        this.procFile = new File(workingDIR + "/pdata/1/proc");

    }


    public LdpfsFourierTransformTool read() throws Exception {

        Acqu acquisition = null;
        Fid fid = null;

        acquisition = new AcquReader(acquFile).read();
        // for the moment set only one processing though it is irrelevant for now
        readPROC(1);

        return new LdpfsFourierTransformTool(new FidReader(new FileInputStream(fidFileName), acquisition).read(), acquisition);


    }


    /**
     * read proc file and extract the paramenters. Note that can be more than one processing.
     *
     * @param processingNb
     */
    private void readPROC(int processingNb) {

    }


}