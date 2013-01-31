package uk.ac.ebi.nmr.fid.io;

import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Fid;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 18:41
 * To change this template use File | Settings | File Templates.
 */
public class FidReaderTest {

    @Test
    public void testReader(){
        File fidFile = new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                "resources/examples/file_formats/bmse000109/1H/fid");

        try {
            Acqu acquisition = new AcquReader("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                    "resources/examples/file_formats/bmse000109/1H/acqu").read();

            //Fid fid = new FidReader(fidFile, acquisition).read();

        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }
}
